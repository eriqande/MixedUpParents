#include <Rcpp.h>
using namespace Rcpp;

// This file has all the stuff in it for calculating the pairwise genotype
// probabilities.


// Enumeration for column indices
enum IGXCol {
  __kIdx = 0,
  __pIdx = 1,
  __cIdx = 2,
  __anck = 3,
  __ancp = 4,
  __start = 5,
  __stop = 6,
  __mIdx = 7,
  __pos = 8,
  __gk = 9,
  __gp = 10
};



// Constants (for now)
const int MaxAnc = 3;   // max number of ancestries
const double epsilon = 0.01;  // genotyping error rate
const int WriteRcouts = 0;




//' Helper function to turn trits to a vector of one or two ancestries
//'
//' As a reminder:
//' - 2 is two copies of the 1-ancestry
//' - 4 is one copy of the 2-ancestry and one of the 1-ancestry
//' - 6 is two copies of the 2-ancestry
//' - 12 is one copy of the 2-ancestry and one copy of the 3-ancestry
//' - 10 is one copy of 3 and one copy of 1
//' - 18 is two copies of the 3-ancestry
//' This is currently configured to deal with up to 3 ancestries (MaxAnc = 3)
//' But this can easily be changed by setting MaxAnc = 4 or 5, etc.
//'
//' It is important to understand that this returns the ancestries in base-0.  So,
//' ancestry 1 is called 0, ancestry 2 is called 1, etc.
//' @export
//' @examples
//' trits <- c(2, 4, 6, 10, 12, 18)
//' names(trits) <- trits
//' lapply(trits, trit2vec)
// [[Rcpp::export]]
 IntegerVector trit2vec(int t) {
    int MaxAnc = 3;
    int b,d;
    IntegerVector ret;

    for(b=MaxAnc-1; b>=0; b--) {
      d = t / pow(3, b);
      if(d>0) ret.push_back(b);
      t -= d * pow(3, b);
    }

    return(ret);
}



// Helper function to calculate SOS probs in the unrelated case
double calculatePAkUnrelated(const IntegerVector& kHap, const NumericMatrix& AD, int k) {
  double PAk_un = 1.0;
  int kd = kHap.length();

  for (int i = 0; i < kd; ++i) {
    PAk_un *= pow(AD(k, kHap[i]), (3 - kd));
  }
  PAk_un *= (1.0 + (kd == 2));

  return PAk_un;
}




// Helper function to calculate genotype probabilities for unrelated case.  We include
// the parent genotype in here, because if it is missing, we don't want to include
// the locus in the kid genotype calculation.
double calculatePGkUnrelated(const IntegerVector& kHap, int kd, int g, int gp, int midx, const IntegerVector& isDiag, const NumericMatrix& AF) {
  double geno_prob = 1.0;
  if(g == -1 || (gp == -1 && isDiag(midx) == 0)) return(1.0);  // missing data case.  Never include if kid is missing.  If parent genotype is missing and it is not a Diagnostic SNP, we also do not include it.
  if (kd == 1) {
    int a = kHap(0);
    double fa = AF(midx, a);
    geno_prob = pow(fa, g) * pow(1 - fa, 2 - g) * (1 + (g == 1));
  } else if (kd == 2) {
    int a = kHap(0);
    int b = kHap(1);
    double fa = AF(midx, a);
    double fb = AF(midx, b);
    geno_prob = (g == 1) ? (fa * (1 - fb) + fb * (1 - fa)) :
      (g == 0) ? ((1 - fa) * (1 - fb)) : (fa * fb);
  } else {
    Rcpp::stop("Invalid kd value in genotype calculation.");
  }
  return geno_prob;
}

// Helper function to calculate the prob of the kid's genotype (gk)
// given a parental relationship and
// conditional on the ancestry of the segment segregated from the parent (As),
// the ancestry of the segment gained from the population (Ap), and the allele
// frequencies, FOR DIAGNOSTIC Markers.  This cycles over all the markers from
// lo to hi and takes the product over all diagnostic markers in that interval.
// Because we have inferred the ancestry of the segments using markers in the interval,
// whether they are missing or not in the parent, we don't actually have to explicitly
// account for the ancestry of the parent at each marker.  And we do not omit
// markers that are not missing in the kid but are missing in the parent, because we
// then it is commensurate to the PGkUnrelated case.
double logPGkParentalDiag(
    int As, // ancestry segregated from the parent
    int Ap, // ancestry received from the population
    int lo, // lo index of markers in the segment
    int hi, // hi index of markers in the segment
    const IntegerVector& isDiag, // vector of 1 of marker is a "diagnostic" marker, and 0 otherwise
    const NumericMatrix& AF, // matrix of all allele frequencies
    const IntegerMatrix& IXG // holds the genotypes
) {
  double geno_prob, ret = 0.0;

  for(int m=lo; m<=hi; m++) {
    int midx = IXG(m, __mIdx);
    if(isDiag(midx) == 1) {  // only accumulate a product for diagnostic markers
      int gk = IXG(m, __gk);
      if(gk == -1) geno_prob = 1.0;
      else
        geno_prob = (gk == 0) ?  (1 - AF(midx, As)) * (1 - AF(midx, Ap)) :
          (gk == 2) ? AF(midx, As) * AF(midx, Ap) :
          (1 - AF(midx, As)) * AF(midx, Ap) + (1 - AF(midx, Ap)) * AF(midx, As);

      if(WriteRcouts) Rcout << "PGk_par:"<< m << " " << midx << " " << isDiag(midx) << " " << AF(midx, As) << " " << (1 - AF(midx, As)) << " " << AF(midx, Ap) << " " << (1 - AF(midx, Ap)) << " "<< gk << " " << "NA" << " " << geno_prob << std::endl;
      // Rcout << "Diag: " << midx << " " <<
      ret += log(geno_prob);
    }
  }
  return(ret);
}


//' Helper function to compute the sum of the log-probs of the offspring genotypes
//' at the VARIABLE markers given the parental genotypes
 double logPGkParentalVar(
     int As, // ancestry segregated from the parent
     int An, // ancestry of the segment *not* segregated by the parent
     int Ap, // ancestry received from the population
     int lo, // lo index of markers in the segment
     int hi, // hi index of markers in the segment
     const IntegerVector& isDiag, // vector of 1 of marker is a "diagnostic" marker, and 0 otherwise
     const NumericMatrix& AF, // matrix of all allele frequencies
     const IntegerMatrix& IXG // holds the genotypes
 ) {
   double geno_prob = 1.0;  // if isDiag == 1 this is what you get
   double ret = 0.0;        // initialize to accumulate a sum of log probs
   double corr = 0.002;     // a correction factor so that allele frequencies of 0 do not occur, so that the log does not blow up

   for(int m=lo; m<=hi; m++) {
     int midx = IXG(m, __mIdx);
     if(isDiag(midx) == 0) {  // only accumulate a product/sum for the variable markers
       int gk = IXG(m, __gk);
       int gp = IXG(m, __gp);
       double f0s = (1.0 - AF(midx, As)) == 0.0 ? corr : (1.0 - AF(midx, As)); // freq of the 0 allele on the segregated ancestry
       double f1s = AF(midx, As) == 0.0 ? corr : AF(midx, As); // freq of the 1 allele on the segregated ancestry
       double f0p = (1.0 - AF(midx, Ap)) == 0.0 ? corr : (1.0 - AF(midx, Ap)); // freq of the 0 allele on the ancestry from the population
       double f1p = AF(midx, Ap) == 0.0 ? corr : AF(midx, Ap); // freq of the 1 allele on the ancestry from the population

       if(gk == -1 || gp == -1) {
         geno_prob = 1.0;  // omit if missing in either the kid or parent
       }
       else {

         //Rcout << "f0s, f1s, f0p, f1p = " << f0s << " " << f1s << " " << f0p << " " << f1p << std::endl;
         // we just break this down into cases based on gp then gk.
         if(gp == 0) {
           if(gk == 2) geno_prob = epsilon;  // Got a 0 from parent so there must have been a genotyping error
           else if(gk == 1) geno_prob = f1p;  // Got the 0 from the parent and the 1 from the population
           else if(gk == 0) geno_prob = f0p;
         }
         else if(gp == 2) {
           if(gk == 0) geno_prob = epsilon;
           else if(gk == 1) geno_prob = f0p;  // Got the 1 from the parent and the 0 from the population
           else if(gk == 2) geno_prob = f1p;
         }
         else if(gp == 1) { // in these cases we must calculate the probability that the 1 allele was segregated from the parent, or not
           double f0n = (1.0 - AF(midx, An)) == 0.0 ? corr : (1.0 - AF(midx, An)); // freq of the 0 allele on the non-segregated ancestry
           double f1n = AF(midx, An) == 0.0 ? corr : AF(midx, An); // freq of the 1 allele on the non-segregated ancestry

           double Pseg1 = (f1s * f0n) / ((f1s * f0n) + (f0s * f1n)); // the probability that the segregated haplotype from the parent carries a 1 allele
           double Pseg0 = (f0s * f1n) / ((f1s * f0n) + (f0s * f1n)); // the probability that the segregated haplotype from the parent carries a 0 allele

           //Rcout << "f0n, f1n, Pseg0, Pseq1 = " << f0n << " " << f1n << " " << Pseg0 << " " << Pseg1 << std::endl;

           if(gk == 0) geno_prob = Pseg0 * f0p;  // the parent segregates a 0 and another 0 is received from the population
           else if(gk == 2) geno_prob = Pseg1 * f1p; // the parent segregates a 1 and another 1 is received from the population
           else if(gk == 1) geno_prob = (Pseg0 * f1p) + (Pseg1 * f0p);  // the parent segs a 0 and the pop a 1, or the parent segs a 1 and the pop a 0
         }
       }
       if(WriteRcouts) Rcout << "PGk_par:"<< m << " " << midx << " " << isDiag(midx)  << " " << f1s << " " << f0s << " " << f1p << " " << f0p << " " << gk << " " << gp << " " << geno_prob << std::endl;
       // Rcout << "geno_prob = " << geno_prob << std::endl;
       ret += log(geno_prob);
     }


   }  // close loop over m
   return(ret);
 }



//' Low level function to compute the pairwise genotype probabilities
//'
//' This should never really be used directly by users.  Rather, use the
//' `pairwise_genotype_probs()` function that calls this.
//' @param IXG an integer matrix of the intersected intervals with columns
//' - kIdx: 0-based index of the kid
//' - pIdx: 0-based index of candidate parent
//' - cIdx: 0-based index of chromosome
//' - anck: ternary trit giving ancestry of kid at this segment
//' - ancp: ternary trit giving ancestry of candidate parent at this segment
//' - start: integer base pair of the start of the segment. Will probably be used just
//'   to know when we are reaching a new segment.
//' - stop: integer base pair of the stop of the segment (probably won't be used)
//' - mIdx: 0-based index of the marker
//' - pos:  base-pair position of the marker (likely will not be used)
//' - gk:   count of 1 alleles in the kid genotype (0, 1, or 2, or -1 for missing data)
//' - gp:   count of 1 alleles in the parent genotype (0, 1, or 2, or -1 for missing data)
//' @param AF a numeric matrix with number of columns equal to the number of ancestries and
//' ordered like A, B, C,... and number of rows equal to the number of variable and diagnostic
//' markers combined.
//' @param isD an integer vector of 0s and 1s with length equal to the number of
//' variable and diagnostic markers combined. A 1 means the marker is a species-diagnostic
//' marker and a 0 means otherwise.
//' @param AD a numeric matrix with number of columns equal to the number of ancestries
//' and number of rows equal to the total number of individuals that got put into the
//' integer representation. Each row sums to one.  These are the admixture fractions.
//' @param debug an IntegerVector, but we only use the first element.  debug(0) == 0 means no debug information.
//' Otherwise the function returns a list with a lot of extra information that can be used
//' to verify that the calculations are being done correctly. If debug(0) == 1, then information
//' is returned at the rate of one row per segment, and no genotype information gets
//' returned.  If debug(0) == 2, then information is returned at one row per marker, and it
//' includes the genotype probabilities for each marker in the unrelated case.
//' If debug(0) == 3 then information is returned with one row for each marker and each
//' distinct value of the ancestry of the segregated haplotype.
//' @export
// [[Rcpp::export]]
List pgp_rcpp(
  IntegerMatrix IXG,
  NumericMatrix AF,
  IntegerVector isD,
  NumericMatrix AD,
  IntegerVector debug
 ) {
  int N = AD.nrow();  // number of individuals
  int M = AF.nrow();  // number of markers
  int R = IXG.nrow(); // total number of rows of genotypes within the IXG matrix
  int lo,hi;  // lo and hi are for the first and last markers in a segment
  int kd, pd; // The number of distinct ancestries in the kid or parent
  double PGk_un;  // For storing the prob of the genotypes given the ancestry (and kid unrelated)
  List ret;  // for returning the values
  int NumParents = unique(IXG(_, __pIdx)).length();
  int theKidx = IXG(0, __kIdx);
  int Prow = 0;  // tells us which row of the return vectors we are on
  double logPunrel = 0.0, logPparental = 0.0;  // for accumulating sums over segments within individuals
  // preallocate to the return vectors
  IntegerVector ret_kIdx(NumParents, theKidx), ret_pIdx(NumParents, -1);
  NumericVector ret_logPunrel(NumParents, 0.0), ret_logPparental(NumParents, 0.0);


  // For storing debug mode stuff and other things
  std::vector<int> rLo, rHi, Gk_list, Gp_list, mIdx_list, isD_list, start_store, stop_store;
  std::vector<double> PAk_unList, PGk_unList, PAk_parList, logProbDiag, logProbVar, geno_prob_vec, fa_vec, fb_vec;

  ////// for storing stuff in debug mode //////
  List Haplist_k, Haplist_p,  ADk_list, popvec;
  IntegerVector ancvec_k, ancvec_p, Asvec, Anvec, Apvec;

  std::vector<int> p2k1_vec, p2k2_vec, kid_hap_from_pop_vec;
  std::vector<double> pseg_prob_vec, pop2kid_prob_vec;

  // initialize things at the first row
  lo=0;

  // keep doing the following until lo == R
  while(lo < R) {
    // If here, we have entered a new segment, so we record things about it and then
    // determine the index of the last marker within it.  We find the end of the
    // segment by figuring out where the kIdx, pIdx, cIdx, or start value changes,
    // or where we are at the very final row (at R - 1)

    // for storing current kIdx (k), pIdx (p), cIdx (c) and start (s)
    int k = IXG(lo, __kIdx), p = IXG(lo, __pIdx), c = IXG(lo, __cIdx), s = IXG(lo, __start);
    for (hi = lo + 1; hi < R; hi++) {
      if (k != IXG(hi, __kIdx) || p != IXG(hi, __pIdx) || c != IXG(hi, __cIdx) || s != IXG(hi, __start)) break;
    }
    hi--; // decrement it back by one to be the actual final row index

    // at this point, lo and hi are the inclusive lo and hi row indexes of the intersected segment

    // vectors for storing the ancestry of the haplotypes in kid and parent
    // If it is length 1, there are two copies of one ancestry.
    IntegerVector kHap = trit2vec(IXG(lo, __anck));
    IntegerVector pHap = trit2vec(IXG(lo, __ancp));

    kd = kHap.length();
    pd = pHap.length();


    // write kidx, parentidx, cidx, lo, and hi
    if(WriteRcouts) Rcout << "KPC:"<< k << " " << p << " " << c << " " << IXG(lo, __anck) << " " << IXG(lo, __ancp) << " " << AD(k, 0) << " " << AD(k, 1) << " " << AD(k, 2) << " " << AD(p, 0) << " " << AD(p, 1) << " " << AD(p, 2) << " " << lo << " " << hi << std::endl;

    if(debug(0) == 1) {  // Just checking to make sure lo and hi are correct
      rLo.push_back(lo);
      rHi.push_back(hi);
      Haplist_k.push_back(kHap);
      Haplist_p.push_back(pHap);
      ancvec_k.push_back(IXG(lo, __anck));
      ancvec_p.push_back(IXG(lo, __ancp));
    }

    // Here, first we can calculate the probability of the kid having
    // the ancestries in kHap from the sack-o-segs model (with no parentage)
    double PAk_un = calculatePAkUnrelated(kHap, AD, k);
    if(WriteRcouts) Rcout << "PAk_un: "<< PAk_un << std::endl;

    // Then after we set the sack-o-segs contribution of the ancestry of the segments
    // we can cycle over the markers from lo to hi and accumulate the product of
    // the genotype probability of each, conditional on the kid's ancestry. (Still in the
    // Unrelated case).
    PGk_un = 1.0;  // to accumulate the product over loci
    double logPGk_un = 0.0; // to accumulate the product as a sum of logs
    for(int m = lo; m <= hi; m++) {
      int g = IXG(m, __gk);
      int gpar = IXG(m, __gp);
      int midx = IXG(m, __mIdx);

      double geno_prob = calculatePGkUnrelated(kHap, kd, g, gpar, midx, isD, AF);
      PGk_un*= geno_prob;  // accumulate the product
      logPGk_un += log(geno_prob);

      if(WriteRcouts) Rcout << "PGk_un:"<< m << " " << midx << " " << isD(midx)  << " "<< g << " " << gpar << " " << geno_prob << std::endl;

      if(debug(0) == 2) {
        rLo.push_back(lo);
        rHi.push_back(hi);
        Haplist_k.push_back(kHap);
        Haplist_p.push_back(pHap);
        ancvec_k.push_back(IXG(lo, __anck));
        ADk_list.push_back(AD(IXG(m, __kIdx),_));
        PAk_unList.push_back(PAk_un);
        mIdx_list.push_back(midx);
        isD_list.push_back(isD(midx));
        Gk_list.push_back(g);
        Gp_list.push_back(gpar);
        geno_prob_vec.push_back(geno_prob);
        PGk_unList.push_back(PGk_un);
      }
    } // close loop over markers
    logPunrel += log(PAk_un) + logPGk_un;

    /////////////////// Done with the Unrelated Case Probability calcs ///////////


    ////////////////////////////////////////////////////////////////
    // Now, in this block we need to cycle over the different possible ancestries
    // that might have been segregated from the parent.  For each one, we record the
    // segregation probs given the sack-o-segs model, and then for each marker in this
    // segment (i.e., cycling from lo to hi, once for each different segregation possibility)
    // we add in the log-prob of the different genotypes found, given that segregation
    // pattern and genotyping error.
    // cycling over the ancestry of the haplotype segregated from the parent
    double sumOverSeggedHaps = 0.0;  // initialize to collect a sum of (unlogged) probabilities.
    for(int segpidx = 0; segpidx < pd; segpidx++)  {
      double pseg_prob = 1.0 / pd; // prob of segregating the ancestry (it is 1/2 if there are two ancestries)
      int As = pHap(segpidx); // the ancestry of what got segregated from the parent.
      int An = (pd == 1) ? As : ((segpidx == 0) ? pHap(1) : pHap(0)) ; // Ancestry of segment that was not segregated from the parent (necessary when parent's genotype is 1)

      IntegerVector fromPop;  // vector to hold the possible haplotypes from the population in the kid

      // if the kid has two doses of a single ancestry, then we know that one of them
      // must have come from the population, regardless of what came from the parent.
      if(kd == 1) {
        fromPop.push_back(kHap(0));
      } else { // otherwise fill fromPop with the ancestries in kHap that are not the same as As
        for(int h=0; h<kHap.length(); h++) {
          if(kHap(h) != As) fromPop.push_back(kHap(h));
        }
      }

      // Rcout << "fromPop.length() = " << fromPop.length() << std::endl;

      // now, cycle over the possible ancestries received from the population
      for(int popsegidx = 0; popsegidx < fromPop.length(); popsegidx++) {
        int Ap = fromPop(popsegidx);
        double pop_seg_factor = 1.0 / fromPop.length();  // if there are two possible ancestries gained from the population we take the average of them by weighting each by 1/2


        // The prob of the kid having gotten those ancestries under the parental hypothesis is the prob of segregating
        // it from the parent (1 or 1/2), times the admixture fraction of the ancestry that came from the population (SOS model)
        // times a factor of 1/2 or 1 which depends on whether the kid has the ancestry of the parent or not (pop_seg_factor).
        double PAk_par = pseg_prob * pop_seg_factor * AD(k, Ap);

        if(WriteRcouts) Rcout << "PAk_par:" << PAk_par << std::endl;

        if(WriteRcouts) Rcout << "AsApAn:" << As << " " << Ap << " " << An << std::endl;


        // get the sum of the log-probabilities of all the genotypes at the diagnostic markers, conditional
        // on the segregated and population-received ancestries
        double logPGk_diag = logPGkParentalDiag(As, Ap, lo, hi, isD, AF, IXG);

        double logPGk_var = logPGkParentalVar(As, An, Ap, lo, hi, isD, AF, IXG);



        // now, we need to multiply the probability of the ancestries segregated to the probs
        // of the variants and sum that over the ancestries segregated.  Currently a simple
        // hack to avoid underflow.
        double TMP = log(PAk_par) + logPGk_diag + logPGk_var;
        sumOverSeggedHaps += (TMP < -600 ? 0 : PAk_par * exp(logPGk_diag) * exp(logPGk_var));

        //Rcout << "logPGk_diag = " << logPGk_diag << "   logPGk_var = " << logPGk_var << "  " << log(PAk_par) + logPGk_diag + logPGk_var << "  " << PAk_par * exp(logPGk_diag) * exp(logPGk_var) << "  " << sumOverSeggedHaps << "  " << logPparental << std::endl;

        // debug stuff
        if(debug(0) == 3) {
          Haplist_k.push_back(kHap);
          Haplist_p.push_back(pHap);
          rLo.push_back(lo);
          rHi.push_back(hi);
          start_store.push_back(IXG(lo, __start));
          stop_store.push_back(IXG(lo, __stop));
          ancvec_k.push_back(IXG(lo, __anck));
          ancvec_p.push_back(IXG(lo, __ancp));
          Asvec.push_back(As);
          Anvec.push_back(An);
          Apvec.push_back(Ap);
          popvec.push_back(fromPop);
          PAk_parList.push_back(PAk_par);
          logProbDiag.push_back(logPGk_diag);
          logProbVar.push_back(logPGk_var);
        }

      } // end loop over popsegidx

    }  // end the for loop over segpidx

    logPparental += (sumOverSeggedHaps == 0 ? 0.0 : log(sumOverSeggedHaps));
    //////////////////////////////////////////////////////////////
    // move onto the next segment
    lo = hi + 1;

    // Accumulate the non-debug results---basically we want one row per candidate
    // parent, so.
    if(p != IXG(lo, __pIdx) || lo >= R) {  // Before advancing to a new parent or finishing up completely, record the results and reset the sums
      ret_pIdx(Prow) = p;
      ret_logPunrel(Prow) = logPunrel;
      ret_logPparental(Prow) = logPparental;

      logPunrel = 0.0;     // reset these to accumulate a sum
      logPparental = 0.0;
      Prow++; // increment to the next return row
    }

  } // end the while loop over the segments

  if(debug(0) == 1) {
    ret = List::create(
      _["Lo"] = wrap(rLo),
      _["Hi"] = wrap(rHi),
      _["anc_k"] = ancvec_k,
      _["haplist_k"] = Haplist_k,
      _["anc_p"] = ancvec_p,
      _["haplist_p"] = Haplist_p
    );
  }
  else if(debug(0) == 2) {
    ret = List::create(
      _["Lo"] = wrap(rLo),
      _["Hi"] = wrap(rHi),
      _["anc_k"] = ancvec_k,
      _["haplist_k"] = Haplist_k,
      _["admix_fracts_k"] = ADk_list,
      _["prob_k_SOS"] = PAk_unList,
      _["mIdx"] = mIdx_list,
      _["isDiag"] = isD_list,
      _["kid_geno"] = Gk_list,
      _["par_geno"] = Gp_list,
      _["geno_prob"] = geno_prob_vec,
      _["prob_kid_geno"] = PGk_unList
    );
  }
  else if(debug(0) == 3) {
    ret = List::create(
      _["Lo"] = wrap(rLo),
      _["Hi"] = wrap(rHi),
      _["start"] = wrap(start_store),
      _["stop"] = wrap(stop_store),
      _["anc_k"] = ancvec_k,
      _["haplist_k"] = Haplist_k,
      _["anc_p"] = ancvec_p,
      _["haplist_p"] = Haplist_p,
      _["Asvec"] = Asvec,
      _["Anvec"] = Anvec,
      _["popvec"] = popvec,
      _["Apvec"] = Apvec,
      _["PAk_par"] = PAk_parList,
      _["LogProbDiag"] = logProbDiag,
      _["LogProbVar"] = logProbVar
    );
  }
  else {
    ret = List::create(
      _["kIdx"] = ret_kIdx,
      _["pIdx"] = ret_pIdx,
      _["probKidUnrel"] = ret_logPunrel,
      _["probKidParental"] = ret_logPparental
    );
  }


  return(ret);

}


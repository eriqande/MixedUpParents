#include <Rcpp.h>
using namespace Rcpp;

// This file has all the stuff in it for calculating the pairwise genotype
// probabilities.


// These are just indexes to use to extract values from columns in IGXmat
#define kIdx  0
#define pIdx  1
#define cIdx  2
#define anck  3
#define ancp  4
#define start 5
#define mIdx  7
#define pos 8
#define gk  9
#define gp  10


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
  int k,p,c,s;  // for storing current kIdx, pIdx, cIdx and start
  int kd, pd; // The number of distinct ancestries in the kid or parent
  int p2k1, p2k2;  // holds the possible ancestries that came from the population (whatever
                      // did not come from the parent).  If there is only one possible ancestry
                      // then that is stored in p2k1, and p2k2 is set to -1. It could be two separate ancestries
                      // if for example, the parent has two doses of 0 and the kid has one
                      // dose each of 1 and 2.  In which case p2k1 = 1 and p2k2 = 2.
  double PAk_un;  // For storing the prob of the ancestry given the kid is unrelated
  double PGk_un;  // For storing the prob of the genotypes given the ancestry (and kid unrelated)
  List ret;  // for returning the values

  ////// for storing stuff in debug mode //////
  std::vector<int> rLo, rHi;  // for testing
  List Haplist_k, Haplist_p,  ADk_list;
  IntegerVector ancvec_k, ancvec_p, segvec_p;

  IntegerVector kHap, pHap;  // vectors for storing the ancestry of the haplotypes in kid and parent
                             // If it is length 1, there are two copies of one ancestry.
  std::vector<int> p2k1_vec, p2k2_vec;
  std::vector<double> PAk_unList, PGk_unList, fa_vec, fb_vec, geno_prob_vec;
  std::vector<int> Gk_list;  // for genotype of the kid
  std::vector<int> mIdx_list;


  // initialize things at the first row
  lo=0;

  // keep doing the following until lo == R
  while(lo < R) {
  // If here, we have entered a new segment, so we record things about it and then
  // determine the index of the last marker within it.  We find the end of the
  // segment by figuring out where the kIdx, pIdx, cIdx, or start value changes,
  // or where we are at the very final row (at R - 1)

  k = IXG(lo, kIdx);
  p = IXG(lo, pIdx);
  c = IXG(lo, cIdx);
  s = IXG(lo, start);

  for(hi=lo+1;hi<R;hi++) {
    if(k != IXG(hi, kIdx) || p != IXG(hi, pIdx) || c != IXG(hi, cIdx) || s != IXG(hi, start)) {
      break;
    }
  }
  hi--;  // decrement it back by one to be the actual final row index

  // at this point, lo and hi are the inclusive lo and hi row indexes of the intersected segment


  ////////////////////////////////////////////////////////////////
  // Now, in this block we need to cycle over the different possible ancestries
  // that might have been segregated from the parent.  For each one, we record the
  // segregation probs given the sack-o-segs model, and then for each marker in this
  // segment (i.e., cycling from lo to hi, once for each different segregation possibility)
  // we add in the log-prob of the different genotypes found, given that segregation
  // pattern and genotyping error.
  kHap = trit2vec(IXG(lo, anck));
  pHap = trit2vec(IXG(lo, ancp));

  kd = kHap.length();
  pd = pHap.length();


  if(debug(0) == 1) {  // all this needs to get imbedded into a loop over markers.
    rLo.push_back(lo);
    rHi.push_back(hi);
    Haplist_k.push_back(kHap);
    Haplist_p.push_back(pHap);
    ancvec_k.push_back(IXG(lo, anck));
    ancvec_p.push_back(IXG(lo, ancp));
  }

  // Here, first we can calculate the probability of the kid having
  // the ancestries in kHap from the sack-o-segs model (with no parentage)
  PAk_un = 1.0;
  for(int i = 0; i< kd; i++) {
    PAk_un *= pow(AD(k, kHap(i)), (3 - kd));  // if kd = 1 then we square the frequency (like a homozygote) otherwise kd==2 and don't
  }
  PAk_un *= (1.0 + kd==2); // If kd == 2 then we add the 2 (like we would for a heterozygote)

  // Then after we set the sack-o-segs contribution of the ancestry of the segments
  // we can cycle over the markers from lo to hi and accumulate the product of
  // the genotype probabilty of each, conditional on the kid's ancestry. (Still in the
  // Unrelated case).
  PGk_un = 1.0;  // to accumulate the product over loci
  for(int m = lo; m <= hi; m++) {
    int a, b;  // for storing ancestries (to make notation easier)
    int g = IXG(m, gk); // for storing the genotype (to make notation easier)
    int midx = IXG(m, mIdx); // for storing the marker index
    double fa=-1.0, fb=-1.0; // to store the freq of the 1 allele
    double geno_prob = 0.0;


    if(kd == 1) {  // simple case---just like a normal population
      a = kHap(0);  // get the index of the ancestry
      fa = AF(midx, a);
      geno_prob =  pow(fa, g) * pow(1 - fa, 2 - g) * (1 + (g == 1)); // standard Binomial with two draws
      PGk_un *= geno_prob;
    } else if(kd == 2) {
      a = kHap(0);  // get the index of the ancestry
      fa = AF(midx, a);
      b = kHap(1);  // get the index of the ancestry
      fb = AF(midx, b);

      if(g == 1) {
        geno_prob *= (fa * (fb - 1)) + (fb * (fa - 1));
      } else if(g == 0) {
        geno_prob *= (1 - fa) * (1 - fb);
      } else {
        geno_prob *= fa * fb;
      }
      PGk_un*= geno_prob;

    } else {
      Rcpp::stop("Neither 1 nor 2 copies of ancestry");
    }

    if(debug(0) == 2) {
      rLo.push_back(lo);
      rHi.push_back(hi);
      Haplist_k.push_back(kHap);
      Haplist_p.push_back(pHap);
      ancvec_k.push_back(IXG(lo, anck));
      ADk_list.push_back(AD(IXG(m, kIdx),_));
      PAk_unList.push_back(PAk_un);
      mIdx_list.push_back(midx);
      fa_vec.push_back(fa);
      fb_vec.push_back(fb);
      Gk_list.push_back(g);
      geno_prob_vec.push_back(geno_prob);
      PGk_unList.push_back(PGk_un);
    }
  } // close loop over markers

  // cycling over the ancestry of the haplotype segregated from the parent
  for(IntegerVector::iterator segp = pHap.begin(); segp != pHap.end(); ++segp) {

    // now figure out what might have been inherited from the population into the kid.
    p2k1 = p2k2 = -1;

    // if the kid has two doses of a single ancestry, then we know that one of them
    // must have come from the population, regardless of what came from the parent.
    if(kd == 1) {
      p2k1 = kHap(0);
    } else {  // otherwise, we set p2k1 and p2k2 to whatever did not get inherited from
              // the parent.  Note that if p2k2 is not -1 then we need to account for
              // the two ancestries in the kid, one of which must be wrong (because the
              // other must have come from the parent, which matches neither of them).
      if(kHap(0) != *segp) {
        p2k1 = kHap(0);
        if(kHap(1) != *segp)
          p2k2 = kHap(2);
      } else{
        p2k1 = kHap(1);
      }
    }

    // Now, segp is the ancestry of the segment passed from the parent and
    // p2k1 is what came from the population.  If p2k2 is not -1, then with prob
    // 1/2 the ancestry from the population is p2k1 and with prob 1/2 it is p2k2
    // and we will have to sum over those cases.

    // Now we can take the product over all the markers on this segment of the
    // probability of the kids genotype given

    // debug stuff
    if(debug(0) == 3) {  // all this needs to get imbedded into a loop over markers.
      Haplist_k.push_back(kHap);
      Haplist_p.push_back(pHap);
      rLo.push_back(lo);
      rHi.push_back(hi);
      ancvec_k.push_back(IXG(lo, anck));
      ancvec_p.push_back(IXG(lo, ancp));
      segvec_p.push_back(*segp);
      p2k1_vec.push_back(p2k1);
      p2k2_vec.push_back(p2k2);
    }

  }  // end the for loop over segp

  //////////////////////////////////////////////////////////////
  // move onto the next segment
  lo = hi + 1;

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
      _["afreq_a"] = fa_vec,
      _["afreq_b"] = fb_vec,
      _["kid_geno"] = Gk_list,
      _["geno_prob"] = geno_prob_vec,
      _["prob_kid_geno"] = PGk_unList
    );
  }
  else if(debug(0) == 3) {
    ret = List::create(
      _["Lo"] = wrap(rLo),
      _["Hi"] = wrap(rHi),
      _["anc_k"] = ancvec_k,
      _["haplist_k"] = Haplist_k,
      _["anc_p"] = ancvec_p,
      _["haplist_p"] = Haplist_p,
      _["segvec_p"] = segvec_p,
      _["p2k1_vec"] = wrap(p2k1_vec),
      _["p2k2_vec"] = wrap(p2k2_vec)
    );
  }


  return(ret);

}

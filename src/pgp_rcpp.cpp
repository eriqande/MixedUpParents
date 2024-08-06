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
#define gk  8
#define gp  9


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
//' @export
//' @examples
//' trits <- c(2, 4, 6, 10, 12, 18)
//' names(trits) <- trits
//' lapply(trits, trit2vec)
// [[Rcpp::export]]
 IntegerVector trit2vec(IntegerVector x) {
    int MaxAnc = 3;
    int t = x[0];
    int b,d;
    IntegerVector ret;

    for(b=MaxAnc-1; b>=0; b--) {
      d = t / pow(3, b);
      if(d>0) ret.push_back(b + 1);
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
//' @export
// [[Rcpp::export]]
DataFrame pgp_rcpp(
  IntegerMatrix IXG,
  NumericMatrix AF,
  IntegerVector isD,
  NumericMatrix AD
 ) {
  int N = AD.nrow();  // number of individuals
  int M = AF.nrow();  // number of markers
  int R = IXG.nrow(); // total number of rows of genotypes within the IXG matrix
  int i,lo,hi;  // lo and hi are for the first and last markers in a segment
  int k,p,c,s;  // for storing current kIdx, pIdx, cIdx and start
  int NewSeg = 1;
  IntegerVector rLo, rHi;

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


  // For testing here (and it checks out!)
  // rLo.push_back(lo);
  // rHi.push_back(hi);

  ////////////////////////////////////////////////////////////////
  // Now, in this block we need to cycle over the different possible ancestries
  // that might have been segregated from the parent.  For each one, we record the
  // segregation probs given the sack-o-segs model, and then for each markers in this
  // segment (i.e., cycling from lo to hi, once for each different segregation possibility)
  // we add in the log-prob of the different genotypes found, given that segregation
  // pattern and genotyping error.





  //////////////////////////////////////////////////////////////
  // move onto the next segment
  lo = hi + 1;

  }

  return(DataFrame::create( _["Lo"] = rLo, _["Hi"] = rHi));
}


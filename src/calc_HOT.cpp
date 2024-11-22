// this is here simply as a convenient place for calculating the HOT
// statistic for comparisons of our method.
#include <Rcpp.h>
using namespace Rcpp;


//' Count the number of homozygous opposites between a kid and candidate parents
//'
//' This is included as a place to have this calculation for comparison, but is
//' not actually used in the MixedUpParents methodology.
//' @param G a integer matrix with loci in rows and samples in columns.  0,1,2 for homozygous
//' reference, heterozygote, and homozygous alternate.  -1 denotes missing data.
//' @param kid  the 0-based subscript of the focal kid
//' @param par integer vector of the 0-based subscripts of the candidate parents you
//' want to compare the kid to.
//' @return Retuns a list that can be easily converted to a tibble.  It has four vector
//' elements: kIdx (0-based index of kid) pIdx (0-based index of candidate parent), num_non_missing,
//' and num_hot.
//' @export
// [[Rcpp::export]]
List calc_HOT(IntegerMatrix G, int kid, IntegerVector par) {
  int n_rows = G.nrow();
  int n_par = par.size();

  // Initialize result vectors
  IntegerVector num_non_missing(n_par);
  IntegerVector num_hot(n_par);

  // Vector to store kIdx, containing the value `kid` repeated `n_par` times
  IntegerVector kIdx(n_par, kid);

  // Clone par to create an independent pIdx
  IntegerVector pIdx = clone(par);

  // Loop over each column index in par
  for (int i = 0; i < n_par; ++i) {
    int p = par[i];
    int non_missing_count = 0;
    int hot_count = 0;

    // Loop over each row to compare kid and p columns
    for (int j = 0; j < n_rows; ++j) {
      int val_kid = G(j, kid);
      int val_p = G(j, p);

      // Check if neither value is -1 (non-missing data)
      if (val_kid != -1 && val_p != -1) {
        non_missing_count++;  // Count valid (non-missing) comparisons

        // Count "HOT" cases: (kid = 0 and p = 2) or (kid = 2 and p = 0)
        if ((val_kid == 0 && val_p == 2) || (val_kid == 2 && val_p == 0)) {
          hot_count++;
        }
      }
    }

    // Store results in vectors
    num_non_missing[i] = non_missing_count;
    num_hot[i] = hot_count;
  }

  // Create a list to return with kIdx and pIdx as the first two elements
  List result = List::create(
    Named("kIdx") = kIdx,
    Named("pIdx") = pIdx,
    Named("num_non_missing") = num_non_missing,
    Named("num_hot") = num_hot
  );

  return result;
}



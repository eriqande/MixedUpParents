#include <Rcpp.h>
using namespace Rcpp;

//' Intersect the ancestry intervals of individuals
//'
//' This is a low level function that should probably not be directly
//' accessed by users.  It is designed to be relatively general and I think
//' it will be better to wrap it within an R function for the user, later.
//' @param grp1 The 0-based indices of the individuals in the first group of indviduals to intersect
//' @param grp2 The 0-based indices of the individuals in the second group
//' @param nC The number of chromosomes
//' @param V1 The Interval Matrix for the group 1 individuals
//' @param V2 The Interval Matrix for the group 2 individuals
//' @param X1 The index for V1
//' @param X2 The index for V2
//' @param nv1 vector of number of intervals for each individual in V1
//' @param nv2 vector of number of intervals for each individual in V1
//' @return  This returns a matrix that is somewhat like the input V1 or V2.
//' It has columns of Idx1 (0), Idx2 (1), cIdx (2), anc1 (3), anc2 (4), start (5), stop (6).  This holds
//' the intersected intervals (or -1 if none exist on the chromosome).
//' @export
// [[Rcpp::export]]
IntegerMatrix intersect_ancestry_intervals_rcpp(
  IntegerVector grp1,
  IntegerVector grp2,
  int nC,
  IntegerMatrix V1,
  IntegerMatrix V2,
  IntegerMatrix X1,
  IntegerMatrix X2,
  IntegerVector nv1,
  IntegerVector nv2
) {
  int n1 = grp1.length();
  int n2 = grp2.length();
  int r = 0;
  int first1, first2, last1, last2; // for the indexes of individuals
  int i1, i2, i, j;
  int OverMax=0;

  // first thing off the bat, we estimate well over the maximum number of rows that
  // we might be returning here to preallocate some memory.
  for(int i=0;i<n1;i++) {
    for(int j=0;j<n2;j++) {
      OverMax += nv1(grp1(i)) + nv2(grp2(j));
    }
  }

  IntegerMatrix ret(OverMax, 7);

  for(int g1=0; g1<n1; g1++) {  // first indiv
    for(int g2=0; g2<n2; g2++) {  // second indiv
      for(int c=0; c<nC; c++) {  // chromosomes
        i1 = grp1(g1);
        i2 = grp2(g2);
        int nada = 1;  // By default we assume there are no intersections, even if there are
                       // intervals on the chromosome.  If we actually do get an intersection,
                       // this gets set to 0.  If it is still 1 by the end, we load up the row
                       // with -1s.

        first1 = X1(i1 * nC + c, 0);
        last1 = X1(i1 * nC + c, 1);

        first2 = X2(i2 * nC + c, 0);
        last2 = X2(i2 * nC + c, 1);

        i = first1;
        j = first2;

        while(i<=last1 && j<=last2) {
          int startA = V1(i, 0);
          int endA = V1(i, 1);
          int startB = V2(j,0);
          int endB = V2(j, 1);
          int anc1 = V1(i, 2);
          int anc2 = V2(j, 2);

          // Return a -1 if either of the chromosomes has no intervals
          if(startA == -1 || endA == -1 || startB == -1 || endB == -1) {
            ret(r, 0) = i1;
            ret(r, 1) = i2;
            ret(r, 2) = c;
            ret(r, 3) = -1;
            ret(r, 4) = -1;
            ret(r, 5) = -1;
            ret(r, 6) = -1;
            r++;  // increment it for the next row
            // set i and j to ensure that we exit the while loop
            i = last1 + 1;
            j = last2 + 1;
          }
          else {
            // Check if intervals intersect
            if (endA >= startB && endB >= startA) {
              // Compute intersection
              double intersectionStart = std::max(startA, startB);
              double intersectionEnd = std::min(endA, endB);

              // Add intersection to result matrix
              ret(r, 0) = i1;
              ret(r, 1) = i2;
              ret(r, 2) = c;
              ret(r, 3) = anc1;
              ret(r, 4) = anc2;
              ret(r, 5) = intersectionStart;
              ret(r, 6) = intersectionEnd;
              r++;
              nada = 0;

              // Move to the next interval in A or B
              if (endA <= endB)
                i++;
              else
                j++;
            }
            else {
              // Move to the next interval in A or B based on the order of their start points
              if (startA < startB)
                i++;
              else
                j++;
            }
          }
        }  // close the while loop

        if(nada == 1) {
          ret(r, 0) = i1;
          ret(r, 1) = i2;
          ret(r, 2) = c;
          ret(r, 3) = -1;
          ret(r, 4) = -1;
          ret(r, 5) = -1;
          ret(r, 6) = -1;
          r++;  // increment it for the next row
        }

      }  // close loop over chromosomes
    } // close loop over g2
  } // close loop over g1

  // resize the matrix to include only the rows used
  ret = ret(Range(0, r - 1), _);

  return(ret);
}


# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Intersect the ancestry intervals of individuals
#'
#' This is a low level function that should probably not be directly
#' accessed by users.  It is designed to be relatively general and I think
#' it will be better to wrap it within an R function for the user, later.
#'
#' The way I am going to use this, mostly, moving forward is to have a single
#' individual in grp 1, and multiple individuals in grp 2 (though we could) have
#' multiple individuals in each.  The important thing here is that the grp1 individual
#' is going to be considered the candidate offspring, and the grp2 individuals are going
#' to be considered as the candidate parents.  However, that asymmetry does not really
#' affect how this function works since I will implement the asymmetrical parts
#' in the wrapper function.
#' @param grp1 The 0-based indices of the individuals in the first group of indviduals to intersect
#' @param grp2 The 0-based indices of the individuals in the second group
#' @param nC The number of chromosomes
#' @param V1 The Interval Matrix for the group 1 individuals
#' @param V2 The Interval Matrix for the group 2 individuals
#' @param X1 The index for V1
#' @param X2 The index for V2
#' @param nv1 vector of number of intervals for each individual in V1
#' @param nv2 vector of number of intervals for each individual in V1
#' @return  This returns a matrix that is somewhat like the input V1 or V2.
#' It has columns of Idx1 (0), Idx2 (1), cIdx (2), anc1 (3), anc2 (4), start (5), stop (6).  This holds
#' the intersected intervals (or -1 if none exist on the chromosome).
#' @export
intersect_ancestry_intervals_rcpp <- function(grp1, grp2, nC, V1, V2, X1, X2, nv1, nv2) {
    .Call('_MixedUpParents_intersect_ancestry_intervals_rcpp', PACKAGE = 'MixedUpParents', grp1, grp2, nC, V1, V2, X1, X2, nv1, nv2)
}

rcpp_hello <- function() {
    .Call('_MixedUpParents_rcpp_hello', PACKAGE = 'MixedUpParents')
}


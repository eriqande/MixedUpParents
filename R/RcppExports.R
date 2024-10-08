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

#' Helper function to turn trits to a vector of one or two ancestries
#'
#' As a reminder:
#' - 2 is two copies of the 1-ancestry
#' - 4 is one copy of the 2-ancestry and one of the 1-ancestry
#' - 6 is two copies of the 2-ancestry
#' - 12 is one copy of the 2-ancestry and one copy of the 3-ancestry
#' - 10 is one copy of 3 and one copy of 1
#' - 18 is two copies of the 3-ancestry
#' This is currently configured to deal with up to 3 ancestries (MaxAnc = 3)
#' But this can easily be changed by setting MaxAnc = 4 or 5, etc.
#'
#' It is important to understand that this returns the ancestries in base-0.  So,
#' ancestry 1 is called 0, ancestry 2 is called 1, etc.
#'
#' @param t the trit to convert to a vector
#' @export
#' @examples
#' trits <- c(2, 4, 6, 10, 12, 18)
#' names(trits) <- trits
#' lapply(trits, trit2vec)
trit2vec <- function(t) {
    .Call('_MixedUpParents_trit2vec', PACKAGE = 'MixedUpParents', t)
}

#' Low level function to compute the pairwise genotype probabilities
#'
#' This should never really be used directly by users.  Rather, use the
#' `pairwise_genotype_probs()` function that calls this.
#' @param IXG an integer matrix of the intersected intervals with columns
#' - kIdx: 0-based index of the kid
#' - pIdx: 0-based index of candidate parent
#' - cIdx: 0-based index of chromosome
#' - anck: ternary trit giving ancestry of kid at this segment
#' - ancp: ternary trit giving ancestry of candidate parent at this segment
#' - start: integer base pair of the start of the segment. Will probably be used just
#'   to know when we are reaching a new segment.
#' - stop: integer base pair of the stop of the segment (probably won't be used)
#' - mIdx: 0-based index of the marker
#' - pos:  base-pair position of the marker (likely will not be used)
#' - gk:   count of 1 alleles in the kid genotype (0, 1, or 2, or -1 for missing data)
#' - gp:   count of 1 alleles in the parent genotype (0, 1, or 2, or -1 for missing data)
#' @param AF a numeric matrix with number of columns equal to the number of ancestries and
#' ordered like A, B, C,... and number of rows equal to the number of variable and diagnostic
#' markers combined.
#' @param isD an integer vector of 0s and 1s with length equal to the number of
#' variable and diagnostic markers combined. A 1 means the marker is a species-diagnostic
#' marker and a 0 means otherwise.
#' @param AD a numeric matrix with number of columns equal to the number of ancestries
#' and number of rows equal to the total number of individuals that got put into the
#' integer representation. Each row sums to one.  These are the admixture fractions.
#' @param debug an IntegerVector, but we only use the first element.  debug(0) == 0 means no debug information.
#' Otherwise the function returns a list with a lot of extra information that can be used
#' to verify that the calculations are being done correctly. If debug(0) == 1, then information
#' is returned at the rate of one row per segment, and no genotype information gets
#' returned.  If debug(0) == 2, then information is returned at one row per marker, and it
#' includes the genotype probabilities for each marker in the unrelated case.
#' If debug(0) == 3 then information is returned with one row for each marker and each
#' distinct value of the ancestry of the segregated haplotype.
#' @export
pgp_rcpp <- function(IXG, AF, isD, AD, debug) {
    .Call('_MixedUpParents_pgp_rcpp', PACKAGE = 'MixedUpParents', IXG, AF, isD, AD, debug)
}

rcpp_hello <- function() {
    .Call('_MixedUpParents_rcpp_hello', PACKAGE = 'MixedUpParents')
}


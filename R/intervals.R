
# Collection of functions for interfacing with the intervals package


#' turn a tibble with start and stop columns into an Intervals object
#' @param tib a tibble with the columns start and stop
#' @keywords internal
#' @export
tib2interval <- function(tib) {
  intervals::Intervals(
    as.matrix(tib),
    closed = c(TRUE, TRUE),
    type = "Z"
  )
}


#' interval union with a shorter name and default no check validity
#'
#' @param a an Interval object
#' @param b another Interval object
#' @keywords internal
#' @export
iu <- function(a, b) {
  intervals::interval_union(a, b, check_valid = FALSE)
}


#' interval intersection with a shorter name and default no check validity
#'
#' @param a an Interval object
#' @param b another Interval object
#' @keywords internal
#' @export
ii <- function(a, b) {
  intervals::interval_intersection(a, b, check_valid = FALSE)
}


#' interval difference with a shorter name and default no check validity
#'
#' @param a an Interval object
#' @param b another Interval object
#' @keywords internal
#' @export
id <- function(a, b) {
  intervals::interval_difference(a, b, check_valid = FALSE)
}

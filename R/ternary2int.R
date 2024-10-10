#' Convert vector representation of a ternary number to a decimal integer
#'
#' Probably won't use this, but I just did it cuz it is a good way
#' to store ancestry types.
#' @param Te A list of vectors, each one a vector of 0's, 1's, or 2's.  The first position is the
#' 0,1,2 position, the second is the 0,3,6 position, the third is the
#' 0,9,18 position, etc.
#' @return An integer vector
#' @details This assumes list input and uses lapply to return a vector
#' @export
#' @examples
#' TeEx <- list(
#'  c(2, 0, 0),
#'  c(0, 2, 0),
#'  c(0, 0, 2),
#'  c(1, 1, 0),
#'  c(1, 0, 1),
#'  c(0, 1, 1)
#' )
#' ternary2int(TeEx)
ternary2int <- function(Te) {
  if(!is.list(Te)) stop("Te must be a list of vectors")
  ret <- lapply(Te, function(x) {
    if(any(x > 2 | x < 0)) {
      warning("values <0 or >2 in ternary ancestry vector")
      return(NA)
    }
    L <- length(x)
    return(sum(x * 3^((1:L) - 1) ) )
  })

  as.integer(unlist(ret))
}

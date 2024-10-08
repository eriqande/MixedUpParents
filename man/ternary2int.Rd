% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ternary2int.R
\name{ternary2int}
\alias{ternary2int}
\title{Convert vector representation of a ternary number to a decimal integer}
\usage{
ternary2int(Te)
}
\arguments{
\item{Te}{A list of vectors, each one a vector of 0's, 1's, or 2's.  The first position is the
0,1,2 position, the second is the 0,3,6 position, the third is the
0,9,18 position, etc.}
}
\value{
An integer vector
}
\description{
Probably won't use this, but I just did it cuz it is a good way
to store ancestry types.
}
\details{
This assumes list input and uses lapply to return a vector
}
\examples{
TeEx <- list(
 c(2, 0, 0),
 c(0, 2, 0),
 c(0, 0, 2),
 c(1, 1, 0),
 c(1, 0, 1),
 c(0, 1, 1)
)
ternary2int(TeEx)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ancestral_n_segments.R
\name{ancestral_n_segments}
\alias{ancestral_n_segments}
\title{Find the segments flanked by common number of copies of ancestry-specific alleles}
\usage{
ancestral_n_segments(D, chrom_levels)
}
\arguments{
\item{D}{a tibble with the columns:
\itemize{
\item \code{diag_spp}: the name (preferably a short code, like three or four letters) of the species
for which this particular position carries diagnostic or species-specific alleles.
\item \code{indiv}: the name of the individual or sample
\item \code{chrom}: chromsome
\item \code{pos}: position on the chromosome
\item \code{n}: number of copies of the species-specific allele (0, 1, 2, or NA)
}}

\item{chrom_levels}{the different chromosome values in order. (This function
will create a factor chrom_f)}
}
\description{
Find the segments flanked by common number of copies of ancestry-specific alleles
}
\examples{
ans <- ancestral_n_segments(diag_markers_10_fish, unique(diag_markers_10_fish$chrom))
}

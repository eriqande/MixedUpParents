


#### Import the pipe operator from magrittr ####
#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @noRd
NULL


#' @importFrom dplyr filter group_by lag mutate n select summarise ungroup
#' @importFrom ggplot2 aes facet_wrap geom_rect geom_segment ggplot scale_colour_manual scale_fill_manual theme_bw
#' @importFrom purrr map
#' @importFrom Rcpp evalCpp
#' @importFrom tibble tibble
#' @importFrom tidyr expand_grid nest unnest
#' @useDynLib MixedUpParents
NULL

# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if (getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "chrom",
      "chrom_f",
      "chrom_top",
      "cu",
      "data",
      "diag_spp",
      "indiv",
      "marker_bottom",
      "marker_top",
      "n_f",
      "pos",
      "seg_bottom",
      "seg_colour_name",
      "seg_tib",
      "seg_top",
      "spp_f",
      "start"
    )
  )
}

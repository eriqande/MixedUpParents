


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


#' @importFrom dplyr arrange count filter group_by lag lead mutate n pull recode select summarise ungroup
#' @importFrom ggplot2 aes facet_wrap geom_rect geom_segment ggplot scale_colour_manual scale_fill_manual theme_bw
#' @importFrom grDevices rainbow
#' @importFrom purrr map map_int
#' @importFrom Rcpp evalCpp
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr complete expand_grid nest pivot_longer pivot_wider replace_na unnest
#' @useDynLib MixedUpParents
NULL




# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if (getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "A0",
      "A0B0C2",
      "A1",
      "A1B1C0",
      "A2",
      "B0",
      "B1",
      "B2",
      "C0",
      "C1",
      "C2",
      "IntervalsEmpty",
      "anc_int",
      "ancestry",
      "block",
      "chrom",
      "chrom_f",
      "chrom_top",
      "cIdx",
      "copy_num",
      "cu",
      "data",
      "diag_spp",
      "drop_it",
      "first",
      "iIdx",
      "indiv",
      "inttib",
      "intv",
      "last",
      "marker_bottom",
      "marker_top",
      "mod_start",
      "n_f",
      "pos",
      "rle_group",
      "row_num",
      "seg_bottom",
      "seg_colour_name",
      "seg_tib",
      "seg_top",
      "spp_f",
      "start"
    )
  )
}





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


#' @importFrom dplyr arrange bind_rows case_when count distinct filter group_by join_by lag lead left_join mutate n pull recode rename select starts_with summarise ungroup
#' @importFrom ggplot2 aes coord_cartesian element_blank facet_wrap geom_col geom_hline geom_line geom_rect geom_segment ggplot guide_legend guides scale_colour_manual scale_fill_brewer scale_fill_manual scale_x_continuous scale_y_continuous theme theme_bw
#' @importFrom grDevices rainbow
#' @importFrom grid unit
#' @importFrom purrr map map2 map_int pmap
#' @importFrom Rcpp evalCpp
#' @importFrom readr read_tsv write_tsv
#' @importFrom stringr str_c str_detect str_split str_sub
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr complete expand_grid extract nest pivot_longer pivot_wider replace_na unnest
#' @useDynLib MixedUpParents
NULL




# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if (getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(".", ".data", "A0", "A0B0C2", "A1", "A1B1C0", "A2", "admix_fract",
      "anc_int", "ancestry", "Ancestry", "Ancestry Copy Numbers", "B0", "B1", "B2", "block",
      "C0", "C1", "C2", "chrom", "chrom_f", "chrom_top", "cIdx", "cmax",
      "cmin", "cn_vec", "copy_num", "cstart", "cu", "data",
      "diag_spp", "dose", "drop_it", "fc0", "fc1", "first", "fract",
      "freq", "geno", "gk", "gp", "group", "idx", "iIdx", "ind_chrom",
      "ind_group", "ind_group_f", "indiv", "indiv_f", "IntervalsEmpty",
      "inttib", "intv", "isCertain", "isDiag", "last", "len", "marker_bottom",
      "marker_top", "markerMash", "mIdx", "mod_start", "n_chr", "n_f",
      "newPs", "pnew1", "pos", "rle_group", "row_num", "seg_bottom",
      "seg_colour_name", "seg_tib", "seg_top", "sp", "spp_f", "spp_idx",
      "start", "tmp", "tot_length", "X1", "xend", "xlab", "xpos", "ylab",
      "yMb_start", "yMb_stop", "ypos")
  )
}


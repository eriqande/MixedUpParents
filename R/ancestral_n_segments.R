

#' Find the segments flanked by common number of copies of ancestry-specific alleles
#'
#' @param D a tibble with the columns:
#' - `diag_spp`: the name (preferably a short code, like three or four letters) of the species
#'   for which this particular position carries diagnostic or species-specific alleles.
#' - `indiv`: the name of the individual or sample
#' - `chrom`: chromsome
#' - `pos`: position on the chromosome
#' - `n`: number of copies of the species-specific allele (0, 1, 2, or NA)
#' @inheritParams chromo_segs
#' @export
#' @examples
#' ans <- ancestral_n_segments(diag_markers_10_fish, unique(diag_markers_10_fish$chrom))
ancestral_n_segments <- function(D, chrom_levels) {
  Dn <- D %>%
    group_by(diag_spp, indiv) %>%
    nest()

  Dn %>%
    mutate(seg_tib = map(data, .f = chromo_segs, chrom_levels = chrom_levels)) %>%
    select(-data) %>%
    unnest(seg_tib) %>%
    ungroup()

}


# let's do it as a function.  This is designed to operate on a sinble
# individual (like nice_fish)
#' Collapse runs of adjacent same-numbers of ancestral alleles into segments
#'
#' This is designed to operate on a single individual and with reference to
#' a single diagnostic species.
#' @param D a tibble with columns:
#' - `chrom`: the name of the chromosome
#' - `pos`: the position along that chromosome
#' - `n`: the number of copies (0, 1, or 2, or NA) of the species-specific alleles
#'    carried by this individual at this position.
#' @param chrom_levels the different chromosome values in order. (This function
#' will create a factor chrom_f)
#' @return This returns a tibble that has the following columns:
#' - `chrom`: the chromosome
#' - `chrom_f`: the chromosome as a factor with levels from `chrom_levels`
#' - `n`: the number of copies of the species-diagnostic allele
#' - `n_f`: n as a factor rather than an integer (useful for plotting with categorical colors)
#' - `start`: the position of the marker at which the segment begins.
#' - `stop`: the position of the marker at which the segement ends
#' @export
chromo_segs <- function(D, chrom_levels) {

  Df <- D %>%
    mutate(chrom_f = factor(chrom, levels = chrom_levels))

  # get a version with no missing data
  Dn <- Df %>%
    filter(!is.na(n))

  # now, we collapse runs of 0, 1, or 2, into segments.  Singleton
  # values, and values that change between points will be "empties"
  Segs <- Dn %>%
    group_by(chrom, chrom_f) %>%
    mutate(
      diff = n != lag(n, default = TRUE),
      cu = cumsum(diff)
    ) %>%
    group_by(chrom, chrom_f, cu) %>%
    summarise(
      n = n[1],
      n_f = factor(n),
      start = pos[1],
      stop = pos[n()],
      .groups = "drop"
    ) %>%
    select(-cu)

  Segs

}

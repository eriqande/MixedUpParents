#' Calculate admixture fractions from the extended segments
#'
#' This just sums up the length of genome from each of the ancestries.
#' @param E the extended segments tibble like what comes out of
#' [extend_ancestral_segments_3()]. Basically a tibble with the columns
#' `indiv`, `chrom_f`, `copy_num`, `start`, `stop`.  The column `copy_num`
#' holds the ancestral copy number indicator which will have values like
#' "A0B1C1" or "A2B0C0", etc.
#' @export
admix_fract_from_extended_segs <- function(E) {

  # we can do this within the tidyverse by breaking up the copy num specifier
  # and unnesting.  It is actually pretty cool...
  E2 <- E %>%
    mutate(
      tmp = strsplit(copy_num,""),
      Ancestry = map(tmp, function(x) x[seq(1, length(x), by = 2)]),
      dose = map(tmp, function(x) as.integer(x[seq(2, length(x), by = 2)])),
    ) %>%
    select(-tmp) %>%
    unnest(cols = c(Ancestry, dose))

  # Now, we just have to do a simple summarise followed by a mutate
  E2 %>%
    group_by(indiv, Ancestry) %>%
    summarise(tot_length = sum(dose * (stop - start))) %>%
    mutate(admix_fract = tot_length / sum(tot_length)) %>%
    ungroup() %>%
    select(-tot_length)

}

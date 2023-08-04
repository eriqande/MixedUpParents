

#' Estimate allele frequencies in the different backgrounds by EM
#'
#' Takes the ancestral copy-number segs tibble and the variable
#' markers data and joins and filters them and then counts up
#' alleles in different ancestry copy number segments and then
#' simply uses R to implement an EM algorithm.
#' @param E the ancestral copy number segs tibble like that produced
#' by extend_ancestral_segments_3().
#' @param V the variable-markers genotypes. Must have indiv IDs and
#' chromosome names that are consistent with those in E.


estimate_allele_freqs_EM <- function(E, V) {

  # get the positions intersected into the segments
  EV <- V %>%
    mutate(chrom_f = factor(chrom, levels = levels(E$chrom_f))) %>%
    left_join(E, ., by = c("chrom_f", "indiv"), relationship = "many-to-many") %>%
    filter(pos >=start & pos <= stop)

  # then count up the number of each of the genotypes (0, 1, 2) in each
  # background
  EV2 <- EV %>%
    rename(geno = n) %>%
    count(chrom_f, pos, copy_num, geno) %>%
    filter(!is.na(geno)) %>%
    mutate(isCertain = str_detect(copy_num, "2"))
}

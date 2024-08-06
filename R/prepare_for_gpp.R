#' Collect inputs for calculating genotype pair probabilities
#'
#' The purpose of this is to collect the necessary inputs and prepare
#' everything to be passed into another function that will do the
#' intersecting of intervals and then the calculation of genotype
#' pair probabilities.
#'
#' The things done by this function are:
#' - Estimating allele frequencies of the variable markers in the
#'   different backgrounds, using an EM algorithm.
#' - Estimating admixture fractions of the individuals using the
#'   extended ancestral segments in each of them.
#' - Creating an integer representation of the extended segments that
#'   can be used for efficiently intersecting segments.
#'
#' The `V` and `M` inputs are largely the "raw" long format inputs,
#' and `E` is what comes out of `extend_ancestral_segments_3()`.
#' @param V Long tibble of variable-marker genotypes with columns
#' `chrom`, `pos`, `indiv`, `n`.
#' @param M long tibble of diagnostic markers with columns `diag_spp`,
#' `chrom`, `pos`, `indiv`, `n`
#' @param E The tibble of ancestral segments as it comes out of
#' `extend_ancestral_segments_3()`.
#' @param diag_err_rate the rate at which diagnostic alleles might not be
#' so diagnostic for a given species.  Defaults to 1 in 200.
#' @return This function retuns a list with a number of components.  This
#' list serves as one of the inputs to `pairwise_genotype_probs()`. The components
#' of the return list are:
#' - `mKey`: a tibble showing the integer equivalent of chrom-key + "-" position
#' for indexing the markers.
#' - `AF`: A tibble of allele frequencies.
#' - ``
#' @export
prepare_for_gpp <- function(
    V,
    M,
    E,
    diag_err_rate = 0.005
) {

  #### Estimate allele frequencies and put them in a wide format ####
  # this gets a wide tibble of the estimated frequency of the 1 allele in the
  # different backgrounds.
  VAF <- estimate_allele_freqs_EM(E = E, V = V)$alle_freqs %>%
    select(chrom_f, pos, Ancestry, pnew1) %>%
    pivot_wider(names_from = Ancestry, values_from = pnew1) %>%
    mutate(chrom = as.character(chrom_f), .before = chrom_f) %>%
    select(-chrom_f)

  #### Prepare "frequencies" of the diagnostic SNPs in the different backgrounds ####
  # This includes the chance for some errors (low frequency of diagnostic alleles on
  # the wrong background, etc.)
  DAF <- M %>%
    distinct(diag_spp, chrom, pos) %>%
    mutate(freq = 1 - diag_err_rate) %>%
    pivot_wider(names_from = diag_spp, values_from = freq, values_fill = diag_err_rate) %>%
    mutate(isDiag = 1L, .after = pos)


  #### Make an integer representation of the extended segments ####
  # This makes list output that has some important things in it, like
  # a key to the integers assigned to individuals and chromosomes that we
  # will use later for creating new integer representations.
  ir_E <- integer_representation_EAS(E)


  #### Combine the VAF and DAF into a single data frame, AF, of allele frequencies ####
  # And we will use chromosome indexes instead of characters here to make sure
  # things will sort out correctly, and also record the 0-based row number of each
  # marker in there and create a marker key (mKey) to record those things.
  AF0 <- VAF %>%
    mutate(isDiag = 0L, .after = pos) %>%
    bind_rows(DAF) %>%
    mutate(cIdx = as.integer(factor(chrom, levels = ir_E$cKey$chrom)) - 1L, .before = chrom) %>%
    select(-chrom) %>%
    arrange(cIdx, pos) %>%
    mutate(
      markerMash = sprintf("%d-%d", cIdx, pos),
      mIdx = as.integer(factor(markerMash, levels = unique(markerMash))) - 1L,
      .after = pos
    )
  mKey <- AF0 %>%
    select(markerMash, mIdx)
  AF <- AF0 %>%
    select(-markerMash)


  #### Combine the variable and the diagnostic marker genotypes into a single tibble ####
  # use the cIdx instead of chrom, and also include an mIdx in there.
  Genos <- M %>%
    select(-diag_spp) %>%
    bind_rows(V) %>%
    mutate(
      iIdx = as.integer(factor(indiv, levels = ir_E$iKey$indiv)) - 1L,
      cIdx = as.integer(factor(chrom, levels = ir_E$cKey$chrom)) - 1L,
      .before = chrom
    ) %>%
    arrange(iIdx, cIdx, pos) %>%
    mutate(
      mIdx = as.integer(factor(sprintf("%d-%d", cIdx, pos), levels = mKey$markerMash)) - 1L,
      .after = cIdx
    ) %>%
    select(-chrom, -indiv)




  #### Estimate admixture fractions of the individuals and include the iIdx with them ####
  AD <- admix_fract_from_extended_segs(E) %>%
    pivot_wider(names_from = Ancestry, values_from = admix_fract) %>%
    mutate(
      iIdx = as.integer(factor(indiv, levels = ir_E$iKey$indiv)) - 1L,
      .after = indiv
    )


  #### Return what we need ####
  list(
    mKey = mKey,
    AFmat = AF %>% select(-cIdx, -pos, -mIdx, -isDiag) %>% as.matrix(),
    ADmat = AD %>% select(-indiv, -iIdx) %>% as.matrix(),
    isDiagVec = AF$isDiag,
    AF = AF,
    AD = AD,
    Genos = Genos,
    ir_E = ir_E
  )
}


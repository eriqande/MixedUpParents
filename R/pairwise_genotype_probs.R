#' Calculate the genotype probabilities for pairs of individuals
#'
#' More description later
#' @param L the list output from `prepare_for_gpp()`.
#' @param kid the ID of the individual (the actual ID, not the iIdx) to be treated
#' as the kid.
#' @param par a vector of the actual IDs of the individuals to be treated as parents
#' @param debug an integer flag  debug[1] == 0 means return no debug information.
#' Otherwise the function returns a list with a lot of extra information that can be used
#' to verify that the calculations are being done correctly. If debug[1] == 1, then information
#' is returned at the rate of one row per segment, and no genotype information gets
#' returned.  If debug[1] == 2, then information is returned at one row per marker, and it
#' includes the genotype probabilities for each marker in the unrelated case.
#' If debug[1] == 3 then information is returned with one row for each marker and each
#' distinct value of the ancestry of the haplotype segregated from the parent. (This
#' allows us to check the calculations for the probability of genotypes given the
#' ancestries and the segregated haplotype). **NOTE:** debug > 1 runs pretty slowly
#' because there are so many push_back calls.  If you are using these, it would be
#' best to pass it farly small data sets (i.e. one candidate kid compare to only a
#' handful of candidate parents).
pairwise_genotype_probs <- function(
  L,
  kid,
  par,
  debug = 0L
) {

  # figure out the iIdx's of the kid and par
  stopifnot(length(kid) == 1)
  tmp <- L$ir_E$iKey %>%
    filter(indiv == as.character(kid))
  if(nrow(tmp) == 0) stop("Individual with ped_id ", kid, " not found!")
  if(nrow(tmp) > 1) stop("Individual with ped_id ", kid, " found in more than one row")
  kIdx <- tmp$idx
  rm(tmp)


  # figure out the iIdx's of the par vector
  tmp <- L$ir_E$iKey %>%
    filter(indiv %in% as.character(par))
  if(nrow(tmp) < length(par)) {
    warning("Not all parent ids found.  Requested ", length(par), " in par and matched ", nrow(tmp), " in ir_E$iKey.")
  }
  pIdx <- tmp$idx
  rm(tmp)

  # intersect the intervals between kid and everyone in par
  IX <- intersect_ancestry_intervals_rcpp(
    grp1 = kIdx,
    grp2 = pIdx,
    nC = L$ir_E$nChrom,
    V1 = L$ir_E$IntervalsMatrix,
    V2 = L$ir_E$IntervalsMatrix,
    X1 = L$ir_E$IndexMatrix,
    X2 = L$ir_E$IndexMatrix,
    nv1 = L$ir_E$nIntvlsVec,
    nv2 = L$ir_E$nIntvlsVec
  )
  colnames(IX) <- c("kIdx", "pIdx", "cIdx", "anck", "ancp", "start", "stop")
  IXtib <- as_tibble(IX)


  # now, join the genotypes in each interval onto the Genos, and turn it
  # into an integer matrix that I can pass to RCpp.
  # This is probably not super efficient, but I think it should be fast
  # enough for what we have going on here, right now.
  IXGmat <- IXtib %>%
    left_join(
      L$Genos,
      by = join_by(kIdx == iIdx, cIdx),
      relationship = "many-to-many"
    ) %>%
    filter(pos > start, pos < stop) %>%
    rename(gk = n) %>%
    left_join(
      L$Genos %>% select(iIdx, mIdx, n) %>% rename(gp = n),
      by = join_by(pIdx == iIdx, mIdx),
      relationship = "many-to-many"
    ) %>%
    mutate(
      gk = ifelse(is.na(gk), -1L, gk),
      gp = ifelse(is.na(gp), -1L, gp)
    ) %>%
    as.matrix()

  # And now I am ready to pass IXGmat, L$AFmat, L$isDiagVec, and L$ADmat off to
  # a function written in RCpp.  I think that is all I will need to calculate the
  # genotype probs for all these pairs.

  # Here are some lines to see what those all look like:

  # head(IXGmat); dim(IXGmat)
  # head(L$AFmat); dim(L$AFmat)
  # head(L$isDiagVec); length(L$isDiagVec)
  # head(L$ADmat); dim(L$ADmat)

  pgp_result <- pgp_rcpp(IXGmat, L$AFmat, L$isDiagVec, L$ADmat, as.integer(debug)) %>% as_tibble()

  return(pgp_result)

}

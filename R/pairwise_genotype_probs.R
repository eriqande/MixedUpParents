#' Calculate the genotype probablities for a pair of individuals
#'
#' More description later
#' @param L the list output from `prepare_for_gpp()`.
#' @param kid the ID of the individual (the actual ID, not the iIdx) to be treated
#' as the kid.
#' @param par a vector of the actual IDs of the individuals to be treated as parents
#' @export
pairwise_genotype_probs <- function(
  L,
  kid,
  par
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

  ret <- pgp_rcpp(IXGmat, L$AFmat, L$isDiagVec, L$ADmat) %>% as_tibble()

  return(ret)

}

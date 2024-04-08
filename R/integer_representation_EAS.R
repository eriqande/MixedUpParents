#' Create an integer representation of everyone's extended ancestry segments
#'
#' This integer representation will be passed to RCpp for fast intersection
#' of segments between pairs of individuals.  We are basically giving individuals
#' and chromosomes base-0 integer indexes.
#' @param E The output of `extend_ancestral_segments_3()`. This table should
#' include all the individuals that you interested in doing pairwise comparisons
#' between.
#' @return Returns a list with six components:
#'   - `nIndiv`: the number of individuals
#'   - `nChrom`: the number of chromosomes
#'   - `iKey`: a tibble giving the base-0 indexes of the individuals (and their names)
#'   - `cKey`: a tibble giving the base-0 indexes of the chromosomes. This and the previous
#'   output are useful for reassociating the sample and chromosome names with later output.
#'   - `IntervalsMatrix`: an integer matrix with five columns: start, stop,
#'   anc_int, iIdx, and cIdx.  Each row corresponds to an interval (or a chromosome with no
#'   intervals) in and individual.  iIdx and cIdx are not really needed since we
#'   also have the index, but I am going to leave those on there for now just
#'   while developing. I put them on the end so that if I drop them it won't change
#'   the subscripting of the columns at all.
#'   -`IndexMatrix`: An integer matrix with columns first, last, iIdx, and cIdx.  Again,
#'   iIdx and cIdx are not needed and I will likely drop them once I get used to
#'   indexing arrays with these.
#'
#' As for indexing arrays, here is how it works.  Imagine you want to get the
#' intervals for sample i on chrom c (where both i and c are zero-based indexes).
#' Let the number of chromosomes be nC and let the IndexMatrix be called `X`, then the
#' lo and hi indexes (base-0) of the rows in IntervalMatrix that hold those intervals
#' are given (in base-0 subscripting, for C++) by `X(i * nC + c, 0)` and `X(i * nC + c, 1)`.
#' By this system we can quickly go through lots of chromosomes in different individuals
#' and intersect the intervals thereon.
#' @export
integer_representation_EAS <- function(E) {

  # get indiv levels
  ilev <- unique(E$indiv)

  # note that chrom_f is already a factor

  # make a tibble that is a key for the indivs
  iKey <- tibble(
    indiv = ilev
  ) %>%
    mutate(idx = 0:(n() - 1))

  # and a tibble that is a key for the chromosomes
  cKey <- tibble(
    chrom = levels(E$chrom_f)
  ) %>%
    mutate(idx = 0:(n() - 1))



  # now make the integer matrix.  We do so on the "completed" tibble
  # so that it will be easy to index into it over indivs and chroms.
  # At the same time we make integer replacements for the copy_num syntax.
  # Note that I am doing these as the integer values of ternary numbers, (i.e.,
  # we can describe the ancestry in base-3 and convert to an integer.)  In other
  # words, because each individuals is diploid, the A B and C values can only be
  # 0, 1, or 2, which is a ternary system.  A is the "0,1,2" digit (i.e., the 3^0 column), while
  # B is the 3^1 column.  0, 3, or 6.  and C in that case is the 3^2 column: 0, 9, or 18.
  # We could write function to return the integer of a ternary vector of 0, 1, and 2,
  # but I will just hard-code it here.
  E2 <- E %>%
    mutate(
      anc_int = recode(
        copy_num,
        A2B0C0 = 2L,
        A0B2C0 = 6L,
        A0B0C2 = 18L,
        A1B1C0 = 4L,
        A1B0C1 = 10L,
        A0B1C1 = 12L
      )
    ) %>%
    tidyr::complete(
      indiv,
      chrom_f,
      fill = list(anc_int = -1L, start = -1, stop = -1)
    ) %>%
    mutate(
      iIdx = as.integer(factor(indiv, levels = ilev)) - 1L,
      cIdx = as.integer(chrom_f) - 1L,
      .before = start
    ) %>%
    mutate(
      start = as.integer(start),
      stop = as.integer(stop)
    ) %>%
    select(iIdx, cIdx, start, stop, anc_int)


  # we are going to take a moment to count how many intervals
  # each individual has, since that can be used to predict
  # an upper bound on the number of intersected intervals between
  # sets of individuals
  numIntvs <- E2 %>%
    count(iIdx) %>%
    pull(n)

  # now, we actually make something to index into the rows here
  # for dealing with multiple intervals in a chromosome
  E_index <- E2 %>%
    mutate(row_num = 0:(n() - 1)) %>%
    group_by(iIdx, cIdx) %>%
    summarise(
      first = first(row_num),
      last = last(row_num),
      .groups = "drop"
    )


  # here is what we return
  list(
    nIndiv = length(ilev),
    nChrom = length(levels(E$chrom_f)),
    nIntvlsVec = numIntvs,
    iKey = iKey,
    cKey = cKey,
    IntervalsMatrix = as.matrix(E2 %>% select(start, stop, anc_int, iIdx, cIdx)),
    IndexMatrix = as.matrix(E_index %>% select(first, last, iIdx, cIdx))
  )
}

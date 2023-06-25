#' Extend segments found by `ancestral_n_segments` for scenarios with three ancestries
#'
#' This function uses the interval arithmetic to come up with segments of different
#' copy numbers of the three ancestries that we will call A, B, and C within the
#' function.  In the process, it stores all the intervals/segments as Intervals
#' objects from the 'intervals' package.
#'
#' This function relies heavily an nesting.
#' @param D a tibble like that which comes out of function ancestral_n_segments.
#' It is crucial that it has these columns:
#' - `diag_spp`: a character vector of the codes/names of which species the marker
#'    is diagnostic for.
#' - `indiv`: the sample ID for the individual.
#' - `chrom_f`: the chromosome as a factor with levels from `chrom_levels`
#' - `n`: the number of copies of the species-diagnostic allele
#' - `start`: the position of the marker at which the segment begins.
#' - `stop`: the position of the marker at which the segement ends
#' @param diag_spp_levels  a vector that will be used to turn the
#' `diag_spp` column into a factor, which will then define species A, B, and C.
#' For the example this package was made for it would be, for example,
#' `c("WCT", "RBT", "YCT")`.
extend_ancestral_segments_3 <- function(D, diag_spp_levels) {

  # make A, B, C for the diag_spp (call it ancestry), then nest
  # by indiv, chrom_f, ancestry, and n, and create Intervals
  D2 <- D %>%
    mutate(
      ancestry = c("A", "B", "C")[as.integer(factor(diag_spp, levels = diag_spp_levels))],
      .after = indiv
    ) %>%
    select(indiv, chrom_f, ancestry, n, start, stop) %>%
    group_by(indiv, chrom_f, ancestry, n) %>%
    nest() %>%
    ungroup() %>%
    arrange(indiv, chrom_f, ancestry, n) %>%
    mutate(intv = map(.x = data, .f = tib2interval)) %>%
    select(-data) %>%
    complete(indiv, chrom_f, ancestry, n) %>%
    mutate(
      intv = map(
        .x = intv,
        .f = function(x) {
          if (is.null(x)) {
            return(IntervalsEmpty)
          } else {
            return(x)
          }
        }
      )
    ) %>%
    pivot_wider(
      names_from = c(ancestry, n),
      names_sep = "",
      values_from = intv
    ) %>%
    select(indiv, chrom_f, A0, A1, A2, B0, B1, B2, C0, C1, C2)

  # and now we do the integer arithmetic to get the different
  # categories like A1B0C1...

  D3 <- D2 %>%
    mutate(
      A1B1C0 = furrr::future_pmap(
        .l = D2,
        .f = function(A0, A1, A2, B0, B1, B2, C0, C1, C2, ...) {
          ii(
            iu(A1, id(ii(B1, C0), iu(A0, A2))),
            iu(B1, id(ii(A1, C0), iu(B0, B2)))
          )
        }),
      A1B0C1 = furrr::future_pmap(
        .l = D2,
        .f = function(A0, A1, A2, B0, B1, B2, C0, C1, C2, ...) {
          ii(
            iu(A1, id(ii(B0, C1), iu(A0, A2))),
            iu(C1, id(ii(A1, B0), iu(C0, C2)))
          )
        }),
      A0B1C1 = furrr::future_pmap(
        .l = D2,
        .f = function(A0, A1, A2, B0, B1, B2, C0, C1, C2, ...) {
          ii(
            iu(B1, id(ii(A0, C1), iu(B0, B2))),
            iu(C1, id(ii(A0, B1), iu(C0, C2)))
          )
        }),
      A2B0C0 = furrr::future_pmap(
        .l = D2,
        .f = function(A0, A1, A2, B0, B1, B2, C0, C1, C2, ...) {
            iu(A2, id(ii(B0, C0), iu(A0, A1)))
        }),
      A0B2C0 = furrr::future_pmap(
        .l = D2,
        .f = function(A0, A1, A2, B0, B1, B2, C0, C1, C2, ...) {
          iu(B2, id(ii(A0, C0), iu(B0, B1)))
        }),
      A0B0C2 = furrr::future_pmap(
          .l = D2,
          .f = function(A0, A1, A2, B0, B1, B2, C0, C1, C2, ...) {
            iu(C2, id(ii(A0, B0), iu(C0, C1)))
          })
    )

  # and now pivot the full-copy-number columns back to long format
  # and remove any that are empty
  D4 <- D3 %>%
    select(-(A0:C2)) %>%
    pivot_longer(
      cols = A1B1C0:A0B0C2,
      names_to = "copy_num",
      values_to = "intv"
    ) %>%
    filter(!(map_int(intv, nrow) == 0))

}

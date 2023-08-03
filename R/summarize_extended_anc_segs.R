

#' Very simple summaries of extended ancestry segments
#'
#' Just a small function to make some summaries of
#' a tibble of extended segments across a lot of
#' individuals.  To be used to get a sense for how
#' much variation there is across individuals in
#' ancestry.
#' @param E the tibble of extended ancestry segments.
#' Has the format of the output of extend_ancestral_segments_3()
#' @export
summarize_extended_anc_segs <- function(E) {

  lengths <- E %>%
    group_by(indiv) %>%
    summarise(tot_length_MB = sum(stop - start) / 1e6)

  fracts <- E %>%
    group_by(indiv, copy_num) %>%
    summarise(len = sum(stop - start), .groups = "drop_last") %>%
    mutate(fract = len / sum(len)) %>%
    ungroup() %>%
    select(-len)

  fracts_wide <- fracts %>%
    pivot_wider(
      names_from = copy_num,
      values_from = fract,
      values_fill = 0.0
    )

  # now, we can do a sort order.  Let's sort on fraction of A first, sort of,
  # and then on B
  A_ord <- fracts_wide %>%
    arrange(desc(A2B0C0), desc(A1B1C0), desc(A1B0C1), desc(A0B2C0), desc(A0B1C1), desc(A0B0C2)) %>%
    pull(indiv)

  catOrd <- c("A2B0C0", "A1B1C0", "A1B0C1", "A0B2C0", "A0B1C1", "A0B0C2")

  # for fun, let's plot these fractions:
  fracts_ord <- fracts %>%
    mutate(
      indiv_f = factor(indiv, levels = A_ord),
      copy_num_f = factor(copy_num, levels = catOrd)
    )

  ggplot(fracts_ord, aes(x = indiv_f, y = fract, fill = copy_num_f)) +
    geom_col()
}

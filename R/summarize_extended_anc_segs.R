

#' Very simple summaries of extended ancestry segments
#'
#' Just a small function to make some summaries of
#' a tibble of extended segments across a lot of
#' individuals.  To be used to get a sense for how
#' much variation there is across individuals in
#' ancestry.
#' @param E the tibble of extended ancestry segments.
#' Has the format of the output of extend_ancestral_segments_3()
#' @param inds_per_row Number of indviuals on each row of the plot
#' @param sort_spec string that would be an input to `arrange()` in order
#' to sort individuals in certain ways.  The default is by the admixture
#' proportion of the "A" species, i.e. `QA`.  You could use `QB` or `QC`,
#' or, if you wanted to sort by the fractions of different full ancestry patterns,
#' you could use, for example:
#' `sort_spec <- c("desc(A2B0C0)", "desc(A1B1C0)", "desc(A1B0C1)", "desc(A0B2C0)", "desc(A0B1C1)", "desc(A0B0C2)")`.
#' But, it is probably best to just sort by admixture fractions. That way the admixture fraction
#' line gets plotted reasonable.
#' @param admixture_line string giving the Q value to plot on the line that overlaps the
#' columns.  Could be `A`, `B`, or `C`.  If it is NULL, then null line is plotted
#' showing the admixture fraction.
#' @param catOrd the order in which you want the colors of the different full
#' ancestry pattern categories to appear in the columns.
#' @export
barplot_ancestry_copy_num_fracts <- function(
  E,
  inds_per_row = 400,
  sort_spec = c("desc(QA)"),
  admixture_line = "A",
  catOrd = c("A2B0C0", "A1B1C0", "A1B0C1", "A0B1C1", "A0B0C2", "A0B2C0")
) {

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
    ) %>%
    mutate(
      QA = (2 * A2B0C0 + A1B1C0 + A1B0C1) / 2,
      QB = (2 * A0B2C0 + A1B1C0 + A0B1C1) / 2,
      QC = (2 * A0B0C2 + A1B0C1 + A0B1C1) / 2
    )

  # now, we can do a sort order.  Let's sort on fraction of A first, sort of,
  # and then on B


  #sort_spec <- c("desc(A2B0C0)", "desc(A1B1C0)", "desc(A1B0C1)", "desc(A0B2C0)", "desc(A0B1C1)", "desc(A0B0C2)")
  sort_spec_p <- rlang::parse_exprs(sort_spec)
  A_ord <- fracts_wide %>%
    arrange(!!! sort_spec_p) %>%
    #arrange(desc(QA)) %>%
    #arrange(desc(A2B0C0), desc(A1B1C0), desc(A1B0C1), desc(A0B2C0), desc(A0B1C1), desc(A0B0C2)) %>%
    pull(indiv)



  # for fun, let's plot these fractions:
  fracts_ord <- fracts %>%
    mutate(
      indiv_f = factor(indiv, levels = A_ord),
      `Ancestry Copy Numbers` = factor(copy_num, levels = catOrd),
      group = as.integer(floor((as.integer(indiv_f) - 1) / inds_per_row)) + 1L,
      xpos = ((as.integer(indiv_f) - 1L) %% inds_per_row) + 1L
    ) %>%
    arrange(indiv_f) %>%
    mutate(
      ind_group = paste0("Samples ", inds_per_row * (group - 1) + 1, " to ", pmin(inds_per_row * group, nlevels(indiv_f))),
      ind_group_f = factor(ind_group, levels = unique(ind_group))
    ) %>%
    select(-ind_group)

  # let's get these set up so that we can plot the admixture fractions as a line
  Q_fracts <- fracts_wide %>%
    select(indiv, starts_with("Q")) %>%
    left_join(fracts_ord %>% select(indiv, group, ind_group_f, xpos) %>% distinct) %>%
    arrange(group, xpos)

  line_var <- paste0("Q", admixture_line)

  tmp <- ggplot(fracts_ord, aes(x = xpos, y = fract)) +
    geom_col(aes(fill = `Ancestry Copy Numbers`)) +
    geom_line(data = Q_fracts, aes(x = xpos, y = .data[[line_var]]), linewidth = 0.4) +
    facet_wrap(~ ind_group_f, ncol = 1) +
    theme_bw() +
    theme(legend.position="top") +
    coord_cartesian(clip = "off", expand = FALSE) +
    ggprism::annotation_ticks(
      sides = "r",
      type = "major",
      outside = TRUE,
      tick.length = unit(3, "pt")
    ) +
    ylab(str_c("Fraction of each full ancestry pattern (colored bars) or admixture fraction of ", admixture_line, " (line)")) +
    xlab("Individual samples") +
    theme(plot.margin = unit(c(5.5, 7.5, 5.5, 5.5), "pt")) +
    scale_fill_brewer(type = "div", palette = "RdYlBu")

  ggsave(tmp,  filename = "grouped-ancestries.pdf", width = 10, height = 10)
}

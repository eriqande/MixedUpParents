

#' Very simple summaries of extended ancestry segments
#'
#' Just a small function to make some summaries of
#' a tibble of extended segments across a lot of
#' individuals.  To be used to get a sense for how
#' much variation there is across individuals in
#' ancestry.
#' @param E the tibble of extended ancestry segments.
#' Has the format of the output of extend_ancestral_segments_3()
#' @param plot_type string saying what type of plot to make: "fracts" just
#' plots the overall fractions of ancestry copy-numbers in each individual,
#' while "segments" actually places the ancestry copy-number segments along
#' the vertical axis according to where they are in the genome.
#' @param diagnostic_markers A tibble of species-diagnostic markers.  Must have
#' a column `diag_spp` which holds a value of A, B, or C, `chrom`, `pos`, `indiv`,
#' and `n`.  `n` is the number of copies of the species diagnostic allele at the
#' locus.  Obviously the chromosome names and the sample names must match those
#' in E.
#' @param variable_markers a tibble that holds the variable markers. Must
#' have a column `var_spp` which holds a value of A, B, or C, naming the
#' species in which the marker is known to be variable.
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
ancestry_barplot <- function(
  E,
  plot_type = "fracts",
  diagnostic_markers = NULL,
  variable_markers = NULL,
  inds_per_row = 400,
  sort_spec = c("desc(QA)"),
  admixture_line = "A",
  catOrd = c("A2B0C0", "A1B1C0", "A1B0C1", "A0B1C1", "A0B0C2", "A0B2C0"),
  fill_scale = scale_fill_brewer(type = "div", palette = "RdYlBu")
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
    left_join(fracts_ord %>% select(indiv, group, ind_group_f, xpos) %>% distinct, by = "indiv") %>%
    arrange(group, xpos)

  line_var <- paste0("Q", admixture_line)

  if(plot_type == "fracts") {
    ret <- ggplot(fracts_ord, aes(x = xpos, y = fract)) +
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
      fill_scale

  } else if(plot_type == "segments") {
    # in this case we apply the previous work to make bar segments
    # first, get the min and max points on each chromosome.  These are
    # where the first and last diagnostic markers are on each
    cEnds <- E %>%
      group_by(chrom_f) %>%
      summarise(
        cmin = min(start),
        cmax = max(stop)
      ) %>%
      ungroup() %>%
      mutate(
        length = cmax - cmin,
        cstart = lag(cumsum(length), default = 0) # this is the "zero" point on each chromosome
      )
    # Now, the yPos for any start or stop of a segment in Mb is found as, for example,
    # yMb = (cstart + start - cmin + 1) / 1e6

    # Note, if we had the length of each chromosome, we could do something a little different
    # as an option, but I think this is going to look best, anyway.

    # join those chrom start and stop values and calculate the yMb_start and yMb_stop
    # for each segment.  Also join on Q_fracts so we know how to order the indivs
    segs_tib <- E %>%
      left_join(
        cEnds,
        by = "chrom_f"
      ) %>%
      mutate(
        yMb_start = (cstart + start - cmin + 1) / 1e6,
        yMb_stop = (cstart + stop - cmin + 1) / 1e6
      ) %>%
      left_join(Q_fracts, by = "indiv") %>%
      arrange(ind_group_f, xpos, yMb_start) %>%
      mutate(`Ancestry Copy Numbers` = factor(copy_num, levels = catOrd))

    bar_halfwidth <- 0.45
    first_layer <- ggplot() +
      geom_rect(
        data = segs_tib,
        mapping = aes(
          xmin = xpos - bar_halfwidth,
          xmax = xpos + bar_halfwidth,
          ymin = yMb_start,
          ymax = yMb_stop,
          fill = `Ancestry Copy Numbers`
        )
      ) +
      geom_hline(data = cEnds, mapping = aes(yintercept = cstart / 1e6), linewidth = 0.05) +
      facet_wrap(~ ind_group_f, ncol = 1) +
      theme_bw() +
      theme(
        legend.position="top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      ) +
      fill_scale +
      scale_x_continuous(expand = c(0.01, 0.01)) +
      scale_y_continuous(expand = c(0, 0)) +
      ylab("Position in the portion of the genome (in megabases) subtended by species-diagnostic markers") +
      xlab("Individual samples")

    ret <- first_layer

    # we place the first segment at xpos - bar_halfwidth + smidge and it extends to
    # xpos - bar_halfwidth + smidge + marker_length
    smidge <- 0.02  # distance between the marker hash marks
    marker_length <- 0.13

    if(!is.null(diagnostic_markers)) {
      DM <- diagnostic_markers %>%
        mutate(
          chrom_f = factor(chrom, levels = levels(E$chrom_f)),
          n_chr = as.character(n)
        ) %>%
        filter(!is.na(n)) %>%
        left_join(cEnds, by = "chrom_f")  %>%
        mutate(ypos = (cstart + pos - cmin + 1) / 1e6) %>%
        left_join(distinct(segs_tib %>% select(indiv, chrom_f, ind_group_f, xpos)), by = join_by(indiv, chrom_f)) %>%
        filter(!is.na(xpos)) %>% # remove individuals that were not in the segs_tib
        arrange(xpos, chrom_f) %>%
        mutate(
          x = xpos - bar_halfwidth + (as.integer(factor(diag_spp, levels = sort(unique(diag_spp)))) - 1L) * (marker_length + smidge),
          xend = x + marker_length
        )

      ret <- first_layer +
        geom_segment(data = DM, aes(x = x, xend = xend, y = ypos, yend = ypos, colour = n_chr), linewidth = 0.01, lineend = "butt") +
        scale_colour_manual(values = c(`0` = "gray90", `1`="gray50", `2` = "black")) +
        guides(colour = guide_legend(title = "Number of species-\ndiagnostic alleles",override.aes = list(linewidth=3)))
    }


    if(!is.null(variable_markers)) {
      x <- 0
    }
    #ggsave(ret, filename = "stacked_segs.pdf", width = 10, height = 30)
  }



  ret
}

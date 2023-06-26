

#' Make a plot of the ancestral segments with 0, 1, or 2 copies of a species-specific allele
#'
#' @param D a tibble like the one that comes out of `ancestral_n_segments()`.
#' @param M a tibble of the markers like in \code{\link{diag_markers_10_fish}}
#' @param spp_levels a vector with the order of the species desired from top
#' to bottom for each chromosome. Defaults to `unique(D$diag_spp)`.
#' @param SF fraction of a chromosome's vertical space occupied by the species specific
#' markers and segments (as opposed to being a blank break between chromosomes)
#' @param MOP "marker overhang part", the fraction of the vertical space devoted to markers
#' and segments that is accounted for by the marker "whisker" that extends from either
#' the top of the bottom of the segment
#' @param E the output from `extend_ancestral_segments_3()` if you want to plot
#' those as well.  Must have been made from D.
#' @export
#' @examples
#' data(diag_markers_10_fish)
#' M <- diag_markers_10_fish
#' D <- ancestral_n_segments(M, unique(M$chrom))
#' g <- plot_ancestral_segs(D, M)
#'
#' # with many individuals, it is recommended to plot this large:
#' # like: ggsave(g, filename = "gplot.pdf", width = 30, height = 40)
#'
#' # If you want to add the extended segments to it:
#' # E <- extend_ancestral_segments_3(D, c("WCT", "RBT", "YCT"))
#' # g2 <- plot_ancestral_segs(D, M, E = E)
#'
#' #' # with many individuals, it is recommended to plot this large:
#' # like: ggsave(g2, filename = "E-gplot.pdf", width = 30, height = 40)
plot_ancestral_segs <- function(
    D,
    M,
    spp_levels = unique(D$diag_spp),
    SF = 0.66,
    MOP = 0.2,
    E = NULL
) {

  if (!setequal(unique(D$diag_spp), spp_levels)) {
    stop("spp_levels must include all diagnostic species values in D, and no more.")
  }
  NS <- length(spp_levels)         # num diagnostic species
  if (NS > 4) {
    stop("Only set up to make colors for <=4 species with species-specific alleles.")
  }
  NC <- length(levels(D$chrom_f))  # num chromosomes

  spp_height <- SF / NS            # total vertical space out of one unit taken up by non-blank space

  # get the segments figured out
  D2 <- D %>%
    mutate(
      spp_f = factor(diag_spp, levels = spp_levels),
      chrom_top = NC - as.integer(chrom_f) + 1,
      marker_top = chrom_top + (1 - as.integer(spp_f)) * spp_height,
      marker_bottom = chrom_top - (as.integer(spp_f)) * spp_height,
      seg_top = marker_top - MOP * spp_height,
      seg_bottom = marker_bottom + MOP * spp_height,
      seg_colour_name = paste(spp_f, n_f, sep = "-")
    )

  # get the markers figured out
  M2 <- M %>%
    mutate(
      chrom_f = factor(chrom, levels = levels(D$chrom_f)),
      spp_f = factor(diag_spp, levels = spp_levels),
      n_f = factor(n),
      chrom_top = NC - as.integer(chrom_f) + 1,
      marker_top = chrom_top + (1 - as.integer(spp_f)) * spp_height,
      marker_bottom = chrom_top - (as.integer(spp_f)) * spp_height,
      seg_colour_name = paste(spp_f, n_f, sep = "-")
    )

  # make a tibble for picking out the colors, Purples, Reds, Greens, Blues from ColorBrewer
  CB_vec <- c(
    RColorBrewer::brewer.pal(3, "Purples"),
    RColorBrewer::brewer.pal(3, "Reds"),
    RColorBrewer::brewer.pal(3, "Greens"),
    RColorBrewer::brewer.pal(3, "Blues")
  )
  segCats <- expand_grid(
    spp_f = factor(spp_levels, levels = spp_levels),
    n_f = factor(c(0, 1, 2))
  )
  seg_colors <- segCats %>%
    mutate(color = CB_vec[1:n()])
  seg_colour_values <- seg_colors$color
  names(seg_colour_values) <- paste(seg_colors$spp_f, seg_colors$n_f, sep = "-")

  g <- ggplot() +
    geom_rect(
      data = D2,
      mapping = aes(xmin = start, xmax = stop, ymin = seg_bottom, ymax = seg_top, fill = seg_colour_name)
    ) +
    geom_segment(
      data = M2,
      mapping = aes(x = pos, xend = pos, y = marker_bottom, yend = marker_top, colour = seg_colour_name)
    ) +
    facet_wrap(~indiv) +
    scale_fill_manual(values = seg_colour_values) +
    scale_colour_manual(values = seg_colour_values) +
    theme_bw()

  if(!is.null(E)) {
    E2 <- E %>%
      mutate(
        chrom_top = NC - as.integer(chrom_f) + 1
      )

    copy_num_color_names <- c("A1B1C0", "A1B0C1", "A0B1C1", "A2B0C0", "A0B2C0", "A0B0C2")
    new_cols <- rainbow(6)
    names(new_cols) <- copy_num_color_names
    all_cols <- c(seg_colour_values, new_cols)
    g2 <- g +
      geom_rect(
        data = E2,
        mapping = aes(xmin = start, xmax = stop, ymin = chrom_top + SF / 16, ymax = chrom_top + SF / 4, fill = copy_num),
        colour = "black",
        linewidth = 0.1
      ) +
      scale_fill_manual(values = all_cols)

    g <- g2

  }

  #ggsave(g, filename = "test.pdf", width = 30, height = 40)
  g

}

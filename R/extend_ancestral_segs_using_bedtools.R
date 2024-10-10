#' extend ancestral segments with three species using bedtools
#'
#' This is a function to extend the ancestral segments using
#' bedtools.  It is more general, and I hope, much faster than
#' the previous approach I used with the intervals package, which is
#' incredibly slow.
#'
#' This is more general because it can accommodate any number of
#' ancestries, though I have only really implemented it with three
#' ancestries.
#'
#' Running only on a single core this function is about 30 times faster
#' than `extend_ancestral_segmentes_3()`, even when the latter function is running on 8 cores.
#' I checked the output on some test cases and this gives nearly identical results to
#' `extend_ancestral_segmentes_3()`, with the differences being one or two small
#' segments, and start and stop positions moved by 1 in some cases.  But overall it is
#' pretty much the same, and WAY faster.
#' @param X a tibble like that which comes out of `ancestral_n_segments()`.
#' It has columns of:
#' - `diag_spp`: character name of the species the 1 allele at the
#'   marker is diagnostic for
#' - `indiv`: character individuals ID
#' - `chrom`: character chromosome name
#' - `chrom_f`: a factor of the chromosome names for proper sorting
#' - `n`: integer number of doses of the ancestry in the segment
#' - `n_f`: a factor of n.
#' - `start`: the starting position (base-1 indexed of the interval)
#' - `stop`: the ending position of the interval.
#' @param species_levels a character vector giving the names of the ancestries
#' in the order that you want them to be considered (essentially the levels of a
#' factor), like `c("WCT", "RBT", "YCT")`
#' @param TMP temp directory in which to run bedtools.  Defaults to a tempfile().
#' Setting it manually can be useful for development and debugging.
#' @export
extend_ancestral_segments_using_bedtools <- function(X, species_levels, TMP = tempfile()) {


  # set up a temp directory for these to work in
  message("Extending ancestral segments using bedtools in temp directory: ", TMP)
  dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

  # use diag_species to make a 1-based integer for each ancestry,
  # and also make the starting position 0-based for bedtools, and
  # make the chromosome names actually be the indiv + the chromosome
  # so that they get associated together properly.
  X2 <- X %>%
    mutate(
      spp_idx = as.integer(factor(diag_spp, levels = species_levels)),
      ind_chrom = paste0(indiv, "_-_", chrom),
      start = ifelse(start > 1, start - 1, start)
    ) %>%
    select(spp_idx, n, ind_chrom, start, stop)


  # nest those up so we have a list column with a tibble for every spp_idx by n
  # combination.
  X3 <- X2 %>%
    arrange(spp_idx, n, ind_chrom, start, stop) %>%
    group_by(spp_idx, n) %>%
    nest()

  # Create empty files for segment types that might not exist. We just create
  # all of them first and then overwrite the ones we really have values for
  # later.
  all_them <- expand_grid(sp = unique(X3$spp_idx), n = 0:2) %>%
    mutate(file = file.path(TMP, paste0("z", sp, "_", n, ".bed")))
  dump <- file.create(all_them$file, showWarnings = FALSE)


  # now we can write out a bed file for every segment of a particular
  # spp_idx by n combination.  we will name these, for example,
  # z1_1 for segments with spp_idx 1 in 1 copy, and z3_2 for those
  # segments that have spp_idx 3 in 2 doses.
  for(r in 1:nrow(X3)) {
    write_tsv(
      X3$data[[r]],
      file = file.path(
        TMP,
        paste0("z", X3$spp_idx[r], "_", X3$n[r], ".bed")
      ),
      col_names = FALSE,
      progress = FALSE
    )
  }

  # those files that we just wrote out are the raw material for doing the intersections
  # and unions that we need to define all these ancestral intervals.  We just need
  # to string together the appropriate intersections and unions via pipes for running
  # bedtools. Like:
  # bedtools intersect -a z2_1.bed -b z3_0.bed | bedtools subtract -a - -b z1_0.bed | bedtools subtract -a - -b z1_2.bed | cat - z1_1.bed | bedtools merge -d 8 | bedtools sort
  # for example.

  # let's start with the easy ones:  A2B0C0, A0B2C0, etc.
  two_doses <- list()
  NS <- length(species_levels)
  for(i in 1:NS) {
    other_spp <- setdiff(1:NS, i)
    # get the intersection of all the other ancestries being zero
    if(length(other_spp) == 1) {
      CALL <- paste0("cd ", TMP, "; cat ", other_spp[1], "_0.bed")
    }
    if(length(other_spp) > 1) {
      CALL <- paste0(
        "cd ", TMP, "; bedtools intersect -a z", other_spp[1], "_0.bed  -b z", other_spp[2], "_0.bed"
      )
      if(length(other_spp) > 2) {
        for(m in 3:length(other_spp)) {
          CALL <- paste0(CALL, " | bedtools intersect -a - -b ", other_spp[m], "_0.bed")
        }
      }
    }
    # then subtract from that the cases where ancestry i is either 0 or 1
    CALL <- paste0(CALL, " | bedtools subtract -a - -b z", i, "_0.bed | bedtools subtract -a - -b z", i, "_1.bed")

    # and finally add all the segments where ancestry i has two copies, merge intervals within 8 bp of one
    # another, and sort them
    CALL <- paste0(CALL, " | cat - z", i, "_2.bed | bedtools sort | bedtools merge -d 8 > double_", i, ".bed")

    system(CALL)

    two_doses[[i]] <- read_tsv(
      file.path(TMP, paste0("double_", i, ".bed")),
      col_names = FALSE,
      progress = FALSE,
      show_col_types = FALSE
    )
  }


  # Cool.  Now we just have to do it for the cases where this is one dose of each of
  # two different ancestries...
  hets <- list()
  for(i in 1:NS) {
    if(i<NS) for(j in (i+1):NS) {
      rem_spp <- setdiff(1:NS, c(i, j))

      #print(c(i, j))

      if(length(rem_spp)==0) {  # if there are only two species
        CALL_i <- paste0("cd ", TMP, "; cat z", i, "_0.bed z", i, "_2.bed | bedtools subtract -a z", j, "_1.bed -b - > i1_term.bed")
        #print(CALL_i)
        CALL_j <- paste0("cd ", TMP, "; cat z", j, "_0.bed z", j, "_2.bed | bedtools subtract -a z", i, "_1.bed -b - > j1_term.bed")
        #print(CALL_j)

        system(CALL_i)
        system(CALL_j)

        CALL_int <- paste0("cd ", TMP, "; bedtools intersect -a i1_term.bed -b j1_term.bed |  bedtools sort | bedtools merge -d 8 > het_result_", i, "_", j, ".bed")
        #print(CALL_int)
        system(CALL_int)

      } else if(length(rem_spp) == 1) {  # three-species case
        CALL_i_ixn <- paste0("cd ", TMP, "; bedtools intersect -a z", j, "_1.bed -b z", rem_spp, "_0.bed > i_rem_spp_intxn.bed")
        CALL_j_ixn <- paste0("cd ", TMP, "; bedtools intersect -a z", i, "_1.bed -b z", rem_spp, "_0.bed > j_rem_spp_intxn.bed")

        system(CALL_i_ixn)
        system(CALL_j_ixn)

        CALL_i <- paste0("cd ", TMP, "; cat z", i, "_0.bed z", i, "_2.bed | bedtools subtract -a i_rem_spp_intxn.bed -b - | cat - z", i, "_1.bed | bedtools sort > i1_term.bed")
        #print(CALL_i)
        CALL_j <- paste0("cd ", TMP, "; cat z", j, "_0.bed z", j, "_2.bed | bedtools subtract -a j_rem_spp_intxn.bed -b - | cat - z", j, "_1.bed | bedtools sort > j1_term.bed")
        #print(CALL_j)

        system(CALL_i)
        system(CALL_j)

        CALL_int <- paste0("cd ", TMP, "; bedtools intersect -a i1_term.bed -b j1_term.bed |  bedtools sort | bedtools merge -d 8 > het_result_", i, "_", j, ".bed")
        #print(CALL_int)
        system(CALL_int)

      } else if(length(rem_spp) > 1) {
        stop("Het ancestry not yet implemented for more than three species. NS = ", NS, ", i = ", i, ", j = ", j, ", rem_spp = ", rem_spp)
      }

      hets[[paste0(i, ",", j)]] <- read_tsv(
        file.path(TMP, paste0("het_result_", i, "_", j, ".bed")),
        col_names = FALSE,
        progress = FALSE,
        show_col_types = FALSE
      )
    }
  }

  # now, we just need to get all the segments back into a tibble
  # with columns: indiv chrom_f copy_num    start     stop
  # and then return that thing.
  bind_rows(
    bind_rows(two_doses, .id = "ancstr"),
    bind_rows(hets, .id = "ancstr")
  ) %>%
    extract(X1, into = c("indiv", "chrom"), regex = "(^.*)_-_(.*)$") %>%
    mutate(
      copy_num = case_when(
        ancstr == "1" ~ "A2B0C0",
        ancstr == "2" ~ "A0B2C0",
        ancstr == "3" ~ "A0B0C2",
        ancstr == "1,2" ~ "A1B1C0",
        ancstr == "2,3" ~ "A0B1C1",
        ancstr == "1,3" ~ "A1B0C1",
        TRUE ~ "SOMETHING_IS_WRONG"
      ),
      chrom_f = factor(chrom, levels = levels(X$chrom_f))
    ) %>%
    rename(
      start = X2,
      stop = X3
    ) %>%
    select(indiv, chrom_f, copy_num, start, stop) %>%
    arrange(indiv, chrom_f, start, stop) %>%
    mutate(stop = stop + 1L)  %>% # for consistency with _some_ of the results from extend_ancestral_segments_3
    filter(stop - start > 10) # toss out segments less than 10 bp long
}

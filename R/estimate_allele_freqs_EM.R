

#' Estimate allele frequencies in the different backgrounds by EM
#'
#' Takes the ancestral copy-number segs tibble and the variable
#' markers data and joins and filters them and then counts up
#' alleles in different ancestry copy number segments and then
#' simply uses R to implement an EM algorithm.
#' @param E the ancestral copy number segs tibble like that produced
#' by extend_ancestral_segments_3().
#' @param V the variable-markers genotypes. Must have indiv IDs and
#' chromosome names that are consistent with those in E.


estimate_allele_freqs_EM <- function(E, V) {

  # get the positions intersected into the segments
  EV <- V %>%
    mutate(chrom_f = factor(chrom, levels = levels(E$chrom_f))) %>%
    left_join(E, ., by = c("chrom_f", "indiv"), relationship = "many-to-many") %>%
    filter(pos >=start & pos <= stop)



  n_anc = nchar(EV$copy_num[1]) / 2  # just to have the number of ancestries

  # then count up the number of each of the genotypes (0, 1, 2) in each
  # background
  EV2 <- EV %>%
    rename(geno = n) %>%
    count(chrom_f, pos, copy_num, geno) %>%
    filter(!is.na(geno)) %>%
    mutate(
      isCertain = str_detect(copy_num, "2") | geno %in% c(0, 2),
      cn_vec = str_split(str_sub(copy_num, 2), pattern = "[A-F]"),
      cn_vec = map(cn_vec, as.integer)
    ) %>%
    mutate(
      counts0 = pmap(., function(geno, n, isCertain, cn_vec, ...) {
        if(isCertain) {ret <- n * ((geno == 0) * cn_vec + (geno == 1) * cn_vec / 2)}
        else ret <- rep(0L, n_anc)
        ret
      }),
      counts1 = pmap(., function(geno, n, isCertain, cn_vec, ...) {
        if(isCertain) {ret <- n * ((geno == 2) * cn_vec + (geno == 1) * cn_vec / 2)}
        else ret <- rep(0L, n_anc)
        ret
      })
    )

  # now, here we can count up the number of 0 and 1 alleles in the three n_anc backgrounds
  # that we know for certain (these are fixed).  And, after doing that, calculate the
  # allele freqs in each of the three backgrounds.
  EV3 <- EV2 %>%
    group_by(chrom_f, pos) %>%
    nest() %>%
    ungroup() %>%
    mutate(
      fc0 = map(data, function(x) colSums(matrix(unlist(x$counts0), byrow = TRUE, ncol = n_anc))), # fixed count 0
      fc1 = map(data, function(x) colSums(matrix(unlist(x$counts1), byrow = TRUE, ncol = n_anc))), # fixed count 1
      UC = map(data, function(x) x %>% filter(!isCertain)), # uncertain combos
    ) %>%
    mutate(
      p_init0 = map2(fc0, fc1, function(x, y) (x + 0.5) / (1 + x + y)),  # initial guess for p0
      p_init1 = map2(fc0, fc1, function(x, y) (y + 0.5) / (1 + x + y))   # initial guess for p1
    )


  # here is a helper function we need for computing some weights
  em_weights <- function(cn_vec, p_old0, p_old1, n, ...) {
    first <- min(which(cn_vec != 0))
    second <- max(which(cn_vec != 0))

    w01 <- p_old0[first] * p_old1[second]
    w10 <- p_old0[second] * p_old1[first]

    p01 <- w01 / (w01 + w10)
    p10 <- w10 / (w01 + w10)

    n0 <- rep(0, n_anc)
    n1 <- rep(0, n_anc)

    n0[c(first, second)] <- c(p01 * n, p10 * n)
    n1[c(first, second)] <- c(p10 * n, p01 * n)

    tibble(n0 = list(n0), n1 = list(n1))

  }

  # Here is a function we can apply to each row of EV3, and each time
  # it implements a simple EM algorithm to estimate the allele frequencies in each
  # background.
  EM <- function(fc0, fc1, UC, p_init0, p_init1, ...) {

    # first off, if there are no ambiguous cases, we just want to return
    # the estimates from the certain cases
    if(nrow(UC) == 0) {
      return(tibble(p_init0 = p_init0, p_init1 = p_init1, pnew0 = p_init0, pnew1 = p_init1, pdiff = 0, iter = 0, n0 = ((0.5 + fc0)), n1 = ((0.5 + fc1)), nu0 = 0, nu1 = 0))
    }
    MaxIter <- 100
    MinIter <- 5
    tol = 0.0001
    iter <- 0

    p_old0 <- p_init0
    p_old1 <- p_init1
    pdiff <- 0.2  # just initialize to a large value so the while loop gets entered

    # add some columns to UC to make the syntax easier
    UCplus <- UC %>% mutate(p_old0 = list(p_old0), p_old1 = list(p_old1))

    while(iter < MinIter && iter <= MaxIter && pdiff > tol) {
      fract_n <- pmap(UCplus, .f = em_weights) %>%
        bind_rows()

      nn0 <- colSums(matrix(unlist(fract_n$n0), byrow = TRUE, ncol = n_anc))
      nn1 <- colSums(matrix(unlist(fract_n$n1), byrow = TRUE, ncol = n_anc))

      denom <- (0.5 + fc0) + nn0 + (0.5 + fc1) + nn1
      p_new0 <- ((0.5 + fc0) + nn0) / denom
      p_new1 <- ((0.5 + fc1) + nn1) / denom

      pdiff <- sum( abs(p_new0 - p_old0) + abs(p_new1 - p_old1) )
      p_old0 <- p_new0
      p_old1 <- p_new1
      UCplus$p_old0 <- list(p_new0)
      UCplus$p_old1 <- list(p_new1)

      iter <- iter + 1
    }
    # return a tibble
    tibble(p_init0 = p_init0, p_init1 = p_init1, pnew0 = p_new0, pnew1 = p_new1, pdiff = pdiff, iter = iter, n0 = 0.5 + fc0, n1 = 0.5 + fc1, nu0 = nn0, nu1 =  nn1)
  }

  EV4 <- EV3 %>%
    mutate(newPs = pmap(., EM))

  # return the information about this stuff in a useful format after unnesting
  alle_freqs <- EV4 %>%
    select(chrom_f, pos, newPs) %>%
    mutate(Ancestry = list(c("A", "B", "C")), .before = newPs) %>%
    unnest(cols = c(Ancestry, newPs))

  geno_counts <- EV2 %>%
    select(chrom_f:n)

  list(
    alle_freqs = alle_freqs,
    geno_counts = geno_counts
  )

}


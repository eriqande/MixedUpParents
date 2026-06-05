MixedUpParents
================

# MixedUpParents

`MixedUpParents` is an R package for genetic parentage assignment in
admixed populations. It was developed for the paper **“MixedUpParents: A
New and Improved Method for Genetic Parentage Assignment in Admixed
Populations.”**

The method combines ancestry information from species-diagnostic markers
with genotype probabilities at markers that are variable across ancestry
backgrounds. It is designed for data sets in which individuals may have
mosaic ancestry and where parentage inference benefits from explicitly
accounting for those ancestry tracts.

## Installation

`MixedUpParents` is currently available only from GitHub:

``` r
remotes::install_github("eriqande/MixedUpParents")
```

The package uses `Rcpp` for computationally intensive
genotype-probability calculations.

## External Dependency

`MixedUpParents` requires `bedtools` as an external command-line
dependency. The `bedtools` executable must be installed and available on
your shell `PATH` before running workflows that extend ancestry
segments.

You can check whether R can find it with:

``` r
Sys.which("bedtools")
```

## Documentation

For a primer on how to use the package, see the online pkgdown
documentation:

<https://eriqande.github.io/MixedUpParents/>

In particular, the example vignette walks through a complete
simulated-data run:

<https://eriqande.github.io/MixedUpParents/articles/example-run-of-MixedUpParents.html>

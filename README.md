
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kifbetadiversity

<!-- badges: start -->

[![R-CMD-check](https://github.com/OdedHol/kifbetadiversity/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/OdedHol/kifbetadiversity/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

You sampled two habitats and want to know which one is more
heterogeneous. The sites in one habitat may look more alike than those
in the other, but is that difference real, or just sampling noise?

kifbetadiversity estimates Whittaker’s multiplicative beta diversity
(`beta_M = gamma / alpha - 1`) with a confidence interval, and tests
whether two habitats differ, following Kiflawi & Spencer (2004,
*Ecology* 85(10):2895-2900). Regional richness (gamma) is estimated with
the first-order jackknife, so the result accounts for species your
survey missed and for the uncertainty in gamma.

## Installation

In RStudio, run:

``` r
install.packages("remotes")   # first time only
remotes::install_github("OdedHol/kifbetadiversity", build_vignettes = TRUE)
library(kifbetadiversity)
```

On Windows, if R asks to install Rtools, you can answer no: the package
has no compiled code.

## From a field CSV to a test, in four steps

The package includes a small made-up survey (30 open and 30 tree sites)
to practise with. Replace the two file paths with your own CSVs.

``` r
library(kifbetadiversity)

# 1. One row per record, plus a list of every site surveyed, including
#    sites where nothing was found (so they still count).
records <- read.csv(system.file("extdata", "example_survey.csv",
                                package = "kifbetadiversity"))
sites   <- read.csv(system.file("extdata", "example_sites.csv",
                                package = "kifbetadiversity"))

# 2. Reshape into a sites x species presence/absence table
res <- as_site_by_species(records, site = "site", species = "species",
                          value = "count", group = "habitat",
                          all_sites = sites)

# 3. beta_M and its 95% confidence interval, per habitat
beta_diversity(res$x, group = res$group)
#>   group n_sites gamma_obs gamma_hat se_gamma_hat    alpha  se_alpha   beta_M
#> 1  open      30        25  29.83333     2.006932 7.800000 0.4723464 2.824786
#> 2  tree      30        38  45.73333     2.381409 6.433333 0.5026228 6.108808
#>   se_ln_beta_M beta_M_lwr beta_M_upr
#> 1    0.1225556   2.221601   3.591742
#> 2    0.1092602   4.931217   7.567613

# 4. Do the two habitats differ?
beta_diversity_test(res$x, group = res$group)
#>   group_1 group_2 beta_M_1 beta_M_2 Z_beta_M     p_beta_M
#> 1    open    tree 2.824786 6.108808 -4.69766 2.631584e-06
```

The vignette explains every column, the assumptions behind the test
(sites must be independent and comparably sampled) and common errors:

``` r
vignette("kifbetadiversity")
```

## Citation

If you use kifbetadiversity, please cite the method it implements:

> Kiflawi, M., & Spencer, M. (2004). Confidence intervals and hypothesis
> testing for beta diversity. *Ecology*, 85(10), 2895-2900.

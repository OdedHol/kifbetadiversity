
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

## Which shape is your data?

Field data usually comes in one of two shapes. Pick the matching
function to prepare it; after that, both paths are the same.

**Shape 1: one row per record** (each row is one species found at one
site)

    site     species  habitat  count
    open_01  sp04     open     3
    open_01  sp10     open     2
    tree_01  sp05     tree     3

``` r
res <- as_site_by_species(records, site = "site", species = "species",
                          value = "count", group = "habitat",
                          all_sites = sites)   # every site surveyed, incl. empty ones
```

**Shape 2: one row per site, one column per species** (often with extra
columns such as an ID, a habitat or a richness total)

    plot_num  habitat  site_id  richness  sp01  sp02  sp03  ...
    1         open     open_01  6         0     0     0     ...
    1         tree     tree_01  7         0     0     0     ...

``` r
res <- as_site_table(tab, site = "site_id", group = "habitat",
                     species = "sp01:sp45")   # only the species columns
```

Every column you give as `species` counts as a species, so leave out IDs
and totals like `plot_num` and `richness`. The site column must be
unique per row: `plot_num` repeats across habitats, so `site_id` is
used.

**Then, for both shapes:**

``` r
beta_diversity(res$x, group = res$group)        # beta and its 95% CI per habitat
beta_diversity_test(res$x, group = res$group)   # do the two habitats differ?
```

## A full example

The package includes a small made-up survey (30 open and 30 tree sites)
in both shapes. Replace the file paths with your own CSV.

``` r
library(kifbetadiversity)

tab <- read.csv(system.file("extdata", "example_site_table.csv",
                            package = "kifbetadiversity"))
res <- as_site_table(tab, site = "site_id", group = "habitat",
                     species = "sp01:sp45")

beta_diversity(res$x, group = res$group)
#>   group n_sites gamma_obs gamma_hat se_gamma_hat    alpha  se_alpha   beta_M
#> 1  open      30        25  29.83333     2.006932 7.800000 0.4723464 2.824786
#> 2  tree      30        38  45.73333     2.381409 6.433333 0.5026228 6.108808
#>   se_ln_beta_M beta_M_lwr beta_M_upr   beta_W beta_W_lwr beta_W_upr
#> 1    0.1225556   2.221601   3.591742 3.824786   3.221601   4.591742
#> 2    0.1092602   4.931217   7.567613 7.108808   5.931217   8.567613
beta_diversity_test(res$x, group = res$group)
#>   group_1 group_2 beta_M_1 beta_M_2 Z_beta_M     p_beta_M
#> 1    open    tree 2.824786 6.108808 -4.69766 2.631584e-06
```

`beta_M` is Whittaker’s beta as in Kiflawi & Spencer (2004),
`gamma_hat / alpha - 1`. `beta_W` is the conventional
`gamma_hat / alpha`, which is simply `beta_M + 1`. Both carry the same
information, and the test is the same for either.

To read your own file, use `tab <- read.csv(file.choose())` and change
the column names to match. The vignette explains every column, the
assumptions behind the test (sites must be independent and comparably
sampled) and common errors:

``` r
vignette("kifbetadiversity")
```

## Citation

If you use kifbetadiversity, please cite the method it implements:

> Kiflawi, M., & Spencer, M. (2004). Confidence intervals and hypothesis
> testing for beta diversity. *Ecology*, 85(10), 2895-2900.

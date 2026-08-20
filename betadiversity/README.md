
<!-- README.md is generated from README.Rmd. Please edit that file -->

# betadiversity

<!-- badges: start -->

<!-- badges: end -->

betadiversity computes Whittaker’s multiplicative (`beta_M`) and Lande’s
additive (`beta_A`) beta diversity from site x species presence/absence
data, with analytical confidence intervals and a two-sample hypothesis
test, following Kiflawi & Spencer (2004, *Ecology* 85(10):2895-2900).
True regional richness (gamma) is estimated via the first-order
jackknife estimator rather than assumed equal to observed richness,
which is what gives the confidence intervals and hypothesis test their
validity.

## Installation

betadiversity is not yet on CRAN or GitHub. Until then, install it from
a local clone with:

``` r
# install.packages("devtools")
devtools::install("path/to/betadiversity")
```

## Example

Raw field data usually comes as one row per observation (a site, a
species, and optionally a detection value and a habitat/category label).
`as_site_by_species()` reshapes that into the sites x species matrix
that `beta_diversity()` expects:

``` r
library(betadiversity)

records <- data.frame(
  site = c("s1", "s1", "s2", "s2", "s3", "s4", "s4"),
  species = c("sp1", "sp2", "sp1", "sp3", "sp2", "sp1", "sp3"),
  seen = c("Y", "Y", "yes", "N", "Y", "Y", "yes"),
  habitat = c("tree", "tree", "tree", "tree", "open", "open", "open")
)

res <- as_site_by_species(records, site = "site", species = "species",
                           value = "seen", group = "habitat",
                           all_sites = paste0("s", 1:4))
res$x
#>    sp1 sp2 sp3
#> s1   1   1   0
#> s2   1   0   0
#> s3   0   1   0
#> s4   1   0   1

beta_diversity(res$x, group = res$group)
#> Warning in .beta_diversity_one(x_g, conf): N = 2 sites is below 10, the
#> smallest sample size Kiflawi & Spencer (2004) validated CI coverage for;
#> results may be unreliable
#> Warning in .beta_diversity_one(x_g, conf): N = 2 sites is below 10, the
#> smallest sample size Kiflawi & Spencer (2004) validated CI coverage for;
#> results may be unreliable
#>   group n_sites gamma_obs gamma_hat se_gamma_hat alpha  se_alpha    beta_M
#> 1  open       2         3       4.5          0.5   1.5 0.3535534 2.0000000
#> 2  tree       2         2       2.5          0.5   1.5 0.3535534 0.6666667
#>   se_ln_beta_M beta_M_lwr beta_M_upr beta_A se_beta_A beta_A_lwr beta_A_upr
#> 1    0.3908680  0.9296566   4.302664      3 0.6123724  1.7997721   4.200228
#> 2    0.7728015  0.1465882   3.031926      1 0.6123724 -0.2002279   2.200228
```

Two independent site-sets can be compared directly with
`beta_diversity_test()`, which returns a two-sample Z-test on both
`beta_M` and `beta_A` (Eqs. 10a-10b of Kiflawi & Spencer 2004):

``` r
set.seed(1)
x <- matrix(rbinom(8 * 10, 1, 0.35), nrow = 8)
y <- matrix(rbinom(8 * 10, 1, 0.35), nrow = 8)
beta_diversity_test(x, y)
#> Warning in .beta_diversity_one(x, conf): N = 8 sites is below 10, the smallest
#> sample size Kiflawi & Spencer (2004) validated CI coverage for; results may be
#> unreliable
#> Warning in .beta_diversity_one(x, conf): N = 8 sites is below 10, the smallest
#> sample size Kiflawi & Spencer (2004) validated CI coverage for; results may be
#> unreliable
#>    Z_beta_M  p_beta_M   Z_beta_A  p_beta_A
#> 1 -1.498708 0.1339493 -0.7695154 0.4415874
```

See `vignette("betadiversity")` for a fuller walkthrough of the
jackknife estimator and what the two beta diversity measures mean.

## Citation

If you use betadiversity, please cite the method it implements:

> Kiflawi, M., & Spencer, M. (2004). Confidence intervals and hypothesis
> testing for beta diversity. *Ecology*, 85(10), 2895-2900.

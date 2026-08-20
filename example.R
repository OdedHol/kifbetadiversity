# Example usage of beta_diversity() and beta_diversity_test()
# (Kiflawi & Spencer 2004 beta diversity, Ecology 85:2895-2900)

devtools::load_all("betadiversity")

## ---- Toy site/species presence-absence pivot table ----
# 7 sites x 9 species. Species sp8 and sp9 are "uniques" (found at only one
# site each), so the jackknife correction to gamma is visibly non-zero.
site_by_species <- matrix(
  c(
    # sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8 sp9
      1,  1,  0,  1,  0,  0,  1,  0,  0,   # site1
      1,  1,  1,  0,  0,  1,  0,  0,  0,   # site2
      0,  1,  1,  1,  1,  0,  0,  0,  0,   # site3
      1,  0,  1,  1,  0,  0,  1,  1,  0,   # site4
      0,  1,  0,  0,  1,  1,  0,  0,  0,   # site5
      1,  1,  1,  1,  1,  0,  0,  0,  0,   # site6
      0,  0,  1,  0,  1,  1,  1,  0,  1    # site7
  ),
  nrow = 7, byrow = TRUE,
  dimnames = list(paste0("site", 1:7), paste0("sp", 1:9))
)

print(site_by_species)

## ---- Point estimates + confidence intervals ----
summary_tbl <- beta_diversity(site_by_species)
print(summary_tbl)

## ---- Manual sanity check against a hand calculation ----
N <- nrow(site_by_species)
gamma_obs_manual <- sum(colSums(site_by_species) > 0)
alpha_manual <- mean(rowSums(site_by_species))
stopifnot(summary_tbl$gamma_obs == gamma_obs_manual)
stopifnot(isTRUE(all.equal(summary_tbl$alpha, alpha_manual)))
cat("\nManual check: gamma_obs =", gamma_obs_manual,
    " alpha =", round(alpha_manual, 3), " -- matches beta_diversity() output\n")

## ---- Wiring check: gamma_hat / se_gamma_hat match Kiflawi & Spencer Eqs. 6a-6b ----
freq_manual <- colSums(site_by_species > 0)
u_manual <- sum(freq_manual == 1)
gamma_hat_manual <- gamma_obs_manual + u_manual * (N - 1) / N
k_i_manual <- rowSums(site_by_species[, freq_manual == 1, drop = FALSE] > 0)
se_gamma_hat_manual <- sqrt((N - 1) / N * (sum(k_i_manual^2) - u_manual^2 / N))
stopifnot(isTRUE(all.equal(summary_tbl$gamma_hat, gamma_hat_manual)))
stopifnot(isTRUE(all.equal(summary_tbl$se_gamma_hat, se_gamma_hat_manual)))
cat("Wiring check: gamma_hat and se_gamma_hat match a manual Eq. 6a/6b computation\n")

# Optional cross-check: gamma_hat (only) still matches vegan::specpool()'s jack1.
# se_gamma_hat is intentionally NOT expected to match vegan's jack1.se, which is
# based on a different derivation (Smith & van Belle 1984) than this package's
# Eq. 6b.
if (requireNamespace("vegan", quietly = TRUE)) {
  sp_check <- vegan::specpool(site_by_species)
  stopifnot(isTRUE(all.equal(summary_tbl$gamma_hat, sp_check$jack1)))
  cat("Optional vegan cross-check: gamma_hat matches vegan::specpool()'s jack1",
      "(se_gamma_hat intentionally differs from jack1.se)\n")
}

## ---- Sanity checks on CI ordering ----
stopifnot(summary_tbl$beta_M_lwr <= summary_tbl$beta_M,
          summary_tbl$beta_M <= summary_tbl$beta_M_upr)
stopifnot(summary_tbl$beta_A_lwr <= summary_tbl$beta_A,
          summary_tbl$beta_A <= summary_tbl$beta_A_upr)
stopifnot(summary_tbl$gamma_hat >= summary_tbl$gamma_obs)
cat("Sanity checks passed: gamma_hat >= gamma_obs, CIs bracket the point estimates\n")

## ---- Two-sample hypothesis test ----
# A second toy site-set (random, same mean incidence as the first) just to
# demonstrate usage. Note site_by_species was hand-built, not drawn from
# the same generative process as site_by_species_2, so a "significant"
# result here reflects that structural difference, not a bug.
set.seed(42)
site_by_species_2 <- matrix(
  rbinom(7 * 9, 1, mean(site_by_species)),
  nrow = 7,
  dimnames = list(paste0("siteB", 1:7), paste0("sp", 1:9))
)

test_tbl <- beta_diversity_test(site_by_species, site_by_species_2)
cat("\nTwo-sample test (usage demo):\n")
print(test_tbl)

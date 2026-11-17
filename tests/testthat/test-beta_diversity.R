test_that("point estimates match the defining algebraic relationships", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  result <- beta_diversity(x)

  expect_equal(result$beta_M, result$gamma_hat / result$alpha - 1)
})

test_that("confidence intervals bracket their point estimates", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  result <- beta_diversity(x)

  expect_true(result$beta_M_lwr <= result$beta_M && result$beta_M <= result$beta_M_upr)
})

test_that("jackknife richness estimate is never below the observed count", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  result <- beta_diversity(x)

  expect_gte(result$gamma_hat, result$gamma_obs)
})

test_that("output matches a hand-verified toy example", {
  # same 7-site x 9-species matrix used in example.R (N < 10, so a sample-size
  # warning is expected here alongside the point estimates)
  x <- matrix(
    c(1,1,0,1,0,0,1,0,0,
      1,1,1,0,0,1,0,0,0,
      0,1,1,1,1,0,0,0,0,
      1,0,1,1,0,0,1,1,0,
      0,1,0,0,1,1,0,0,0,
      1,1,1,1,1,0,0,0,0,
      0,0,1,0,1,1,1,0,1),
    nrow = 7, byrow = TRUE
  )
  expect_warning(result <- beta_diversity(x), "Only N = ")

  expect_equal(result$n_sites, 7)
  expect_equal(result$gamma_obs, 9)             # hand-counted: all 9 species appear somewhere
  expect_equal(result$alpha, 30 / 7)             # hand-counted: 30 total presences / 7 sites
  expect_equal(result$gamma_hat, 9 + 2 * 6 / 7)  # Eq. 6a: u = 2 unique species (sp8, sp9)
  expect_equal(result$se_gamma_hat, 1.106567, tolerance = 1e-6)  # Eq. 6b, hand-verified
})

test_that("gamma_hat matches vegan::specpool()'s jack1; se_gamma_hat intentionally diverges from jack1.se", {
  testthat::skip_if_not_installed("vegan")
  # toy matrix with u = 2 unique species (sp8, sp9), so se_gamma_hat is
  # actually nonzero and the two formulas have something to disagree on
  # (with u = 0 both formulas trivially agree that se_gamma_hat is 0).
  x <- matrix(
    c(1,1,0,1,0,0,1,0,0,
      1,1,1,0,0,1,0,0,0,
      0,1,1,1,1,0,0,0,0,
      1,0,1,1,0,0,1,1,0,
      0,1,0,0,1,1,0,0,0,
      1,1,1,1,1,0,0,0,0,
      0,0,1,0,1,1,1,0,1),
    nrow = 7, byrow = TRUE
  )
  sp <- vegan::specpool(x)
  expect_warning(result <- beta_diversity(x), "Only N = ")

  expect_equal(result$gamma_hat, sp$jack1)
  # se_gamma_hat now follows the paper's own Eq. 6b, not vegan's Smith & van
  # Belle (1984) jack1.se -- this divergence is intentional (regression guard
  # against silently reverting the fix).
  expect_false(isTRUE(all.equal(result$se_gamma_hat, sp$jack1.se)))
})

test_that("se_gamma_hat is exactly 0 when there are no unique species", {
  # every species present at exactly 2 sites: no species is a "unique"
  x <- matrix(0, nrow = 4, ncol = 3)
  x[1:2, 1] <- 1
  x[2:3, 2] <- 1
  x[3:4, 3] <- 1
  expect_warning(result <- beta_diversity(x), "Only N = ")

  expect_equal(result$se_gamma_hat, 0)
  expect_equal(result$gamma_hat, result$gamma_obs)
})

test_that("beta_diversity_test() on identical data gives a null result", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  result <- beta_diversity_test(x, x)

  expect_equal(result$Z_beta_M, 0)
  expect_equal(result$p_beta_M, 1)
})

## ---- input validation ----

test_that("non-numeric, non-logical input errors", {
  x <- matrix(c("a", "b", "c", "d"), nrow = 2)
  expect_error(beta_diversity(x), "only numbers or TRUE/FALSE")
})

test_that("negative values error", {
  x <- matrix(c(1, -1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1), nrow = 3)
  expect_error(beta_diversity(x), "negative")
})

test_that("infinite values error", {
  x <- matrix(c(1, Inf, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1), nrow = 3)
  expect_error(beta_diversity(x), "infinite")
})

test_that("NA values are treated as absence, with a warning", {
  set.seed(2)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  x_na <- x
  x_na[1, 1] <- NA
  x_zero <- x
  x_zero[1, 1] <- 0

  expect_warning(result_na <- beta_diversity(x_na), "blank/missing")
  result_zero <- beta_diversity(x_zero)
  expect_equal(result_na, result_zero)
})

test_that("fewer than 2 sites errors", {
  x <- matrix(c(1, 0, 1), nrow = 1)
  expect_error(beta_diversity(x), "at least 2 sites")
})

test_that("zero species (columns) errors", {
  x <- matrix(nrow = 3, ncol = 0)
  expect_error(beta_diversity(x), "at least 1 species")
})

test_that("2-9 sites warns about sample size but still returns a result", {
  set.seed(3)
  x <- matrix(rbinom(5 * 10, 1, 0.35), nrow = 5)
  expect_warning(result <- beta_diversity(x), "Only N = ")
  expect_equal(result$n_sites, 5)
})

test_that("10 or more sites does not warn about sample size", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_no_warning(beta_diversity(x))
})

test_that("abundance data is thresholded the same as presence/absence", {
  set.seed(4)
  counts <- matrix(rpois(12 * 10, 1.5), nrow = 12)
  incidence <- (counts > 0) * 1

  expect_equal(beta_diversity(counts), beta_diversity(incidence))
})

test_that("identical species lists at every site error, naming the real cause", {
  x <- matrix(1, nrow = 10, ncol = 4)  # gamma_hat == alpha exactly: no uniques, no turnover
  expect_error(beta_diversity(x), "same 4 species", class = "kifbeta_no_turnover")
  expect_error(beta_diversity(x), "no turnover")
})

test_that("a species present at every site is fine as long as there is turnover", {
  # sp1 is everywhere, but the other species vary between sites
  x <- cbind(1, diag(10))
  result <- beta_diversity(x)
  expect_equal(result$beta_M, (11 + 10 * 9 / 10) / 2 - 1)  # Eq. 6a: u = 10
})

test_that("no species recorded at all errors, naming that cause", {
  x <- matrix(0, nrow = 10, ncol = 4)
  expect_error(beta_diversity(x), "no species were recorded",
               class = "kifbeta_no_turnover")
})

## ---- grouping ----

test_that("ungrouped call is unchanged numerically and reports group == 'All'", {
  x <- matrix(
    c(1,1,0,1,0,0,1,0,0,
      1,1,1,0,0,1,0,0,0,
      0,1,1,1,1,0,0,0,0,
      1,0,1,1,0,0,1,1,0,
      0,1,0,0,1,1,0,0,0,
      1,1,1,1,1,0,0,0,0,
      0,0,1,0,1,1,1,0,1),
    nrow = 7, byrow = TRUE
  )
  expect_warning(result <- beta_diversity(x), "Only N = ")

  expect_equal(result$group, "All")
  expect_equal(result$gamma_obs, 9)
  expect_equal(result$alpha, 30 / 7)
})

test_that("group length mismatch errors", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, group = rep("a", 5)), "length")
})

test_that("more than 5 groups errors", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, group = letters[1:12]), "at most 5 groups")
})

test_that("NA and blank group labels form their own NA group, last, with a warning naming the sites", {
  set.seed(1)
  # 12 sites per group; 2 of each become unlabelled (one NA, one blank),
  # leaving 10 labelled sites per group and 4 unlabelled ones.
  x <- matrix(rbinom(24 * 10, 1, 0.35), nrow = 24,
              dimnames = list(paste0("site", 1:24), NULL))
  group <- rep(c("control", "npk"), each = 12)
  group_na <- group
  group_na[c(1, 13)] <- NA
  group_na[c(2, 14)] <- "  "
  unlabelled <- c(1, 2, 13, 14)

  warns <- character(0)
  result <- withCallingHandlers(
    beta_diversity(x, group = group_na),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("site1, site2, site13, site14", warns)))
  expect_true(any(grepl("shown as NA", warns)))

  expect_equal(result$group, c("control", "npk", NA))
  expect_false("" %in% result$group)

  cols <- setdiff(names(result), "group")
  direct_na <- suppressWarnings(
    beta_diversity(x[unlabelled, , drop = FALSE])
  )[, cols]
  expect_equal(as.numeric(result[3, cols]), as.numeric(direct_na))
  direct_control <- beta_diversity(x[3:12, , drop = FALSE])[, cols]
  expect_equal(as.numeric(result[1, cols]), as.numeric(direct_control))
})

test_that("a single unlabelled site is skipped like any undersized group", {
  set.seed(1)
  x <- matrix(rbinom(21 * 10, 1, 0.35), nrow = 21)
  group <- c(rep(c("control", "npk"), each = 10), NA)
  result <- suppressWarnings(beta_diversity(x, group = group))
  expect_equal(result$group, c("control", "npk"))
  expect_warning(
    expect_warning(beta_diversity(x, group = group), "fewer than 2 sites"),
    "blank or missing group label: row 21"
  )
})

test_that("a group with no turnover is skipped with a warning; other groups still return", {
  set.seed(1)
  x <- rbind(matrix(rbinom(10 * 6, 1, 0.4), nrow = 10),
             matrix(1, nrow = 10, ncol = 6))
  group <- rep(c("varied", "uniform"), each = 10)
  expect_warning(result <- beta_diversity(x, group = group),
                 "'uniform' was skipped")
  expect_equal(result$group, "varied")
})

test_that("an undersized group is skipped with a warning; other groups still return", {
  set.seed(1)
  x <- matrix(rbinom(21 * 10, 1, 0.35), nrow = 21)
  group <- c(rep("control", 10), rep("npk", 10), "fire")  # fire has 1 site

  expect_warning(result <- beta_diversity(x, group = group), "skipped")

  expect_setequal(result$group, c("control", "npk"))

  cols <- setdiff(names(result), "group")
  direct_control <- beta_diversity(x[group == "control", , drop = FALSE])[, cols]
  actual_control <- result[result$group == "control", cols]
  expect_equal(as.numeric(actual_control), as.numeric(direct_control))
})

test_that("all groups undersized errors", {
  set.seed(1)
  x <- matrix(rbinom(2 * 5, 1, 0.35), nrow = 2)
  group <- c("a", "b")  # 1 site each

  # each undersized group also emits its own "skipped" warning on the way to
  # the final error; that per-group behavior is already covered by the
  # "undersized group is skipped" test above, so it's suppressed here.
  suppressWarnings(
    expect_error(beta_diversity(x, group = group), "No group could be analysed")
  )
})

test_that("grouped results match calling beta_diversity() directly on each subset", {
  set.seed(1)
  x <- matrix(rbinom(24 * 10, 1, 0.35), nrow = 24)
  group <- rep(c("control", "npk"), each = 12)

  result <- beta_diversity(x, group = group)
  cols <- setdiff(names(result), "group")

  for (g in c("control", "npk")) {
    direct <- beta_diversity(x[group == g, , drop = FALSE])[, cols]
    actual <- result[result$group == g, cols]
    expect_equal(as.numeric(actual), as.numeric(direct))
  }
})

test_that("group as a factor behaves the same as group as a character vector", {
  set.seed(1)
  x <- matrix(rbinom(24 * 10, 1, 0.35), nrow = 24)
  group_chr <- rep(c("control", "npk"), each = 12)
  group_fct <- factor(group_chr)

  expect_equal(beta_diversity(x, group = group_chr), beta_diversity(x, group = group_fct))
})

## ---- conf validation ----

test_that("non-numeric conf errors", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, conf = "0.95"), "single number")
})

test_that("conf of length > 1 errors", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, conf = c(0.9, 0.95)), "single number")
})

test_that("conf = NA errors", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, conf = NA), "single number")
})

test_that("conf > 1 errors", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, conf = 1.5), "strictly between 0 and 1")
})

test_that("conf < 0 errors (regression test for the silently-inverted-CI bug)", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, conf = -0.2), "strictly between 0 and 1")
})

test_that("conf given as a percentage suggests the proportion form", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, conf = 95), "Did you mean conf = 0.95")
})

test_that("conf = 0 and conf = 1 (boundaries) error", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_error(beta_diversity(x, conf = 0), "strictly between 0 and 1")
  expect_error(beta_diversity(x, conf = 1), "strictly between 0 and 1")
})

test_that("valid non-default conf values work and change interval width as expected", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)

  r90 <- beta_diversity(x, conf = 0.90)
  r95 <- beta_diversity(x, conf = 0.95)
  r99 <- beta_diversity(x, conf = 0.99)

  width90 <- r90$beta_M_upr - r90$beta_M_lwr
  width95 <- r95$beta_M_upr - r95$beta_M_lwr
  width99 <- r99$beta_M_upr - r99$beta_M_lwr

  expect_true(width90 < width95)
  expect_true(width95 < width99)
})

## ---- output shape (beta_M only) ----

test_that("beta_diversity() returns beta_M columns only (no beta_A)", {
  set.seed(1)
  x <- matrix(rbinom(12 * 10, 1, 0.35), nrow = 12)
  expect_equal(names(beta_diversity(x)),
               c("group", "n_sites", "gamma_obs", "gamma_hat", "se_gamma_hat",
                 "alpha", "se_alpha", "beta_M", "se_ln_beta_M",
                 "beta_M_lwr", "beta_M_upr"))
  expect_false(any(grepl("beta_A", names(beta_diversity_test(x, x)))))
})

## ---- wide CSV read without row.names ----

test_that("a data frame with a text column names it and suggests row.names = 1", {
  csv <- "site,sp1,sp2,sp3\ns1,1,0,1\ns2,0,1,1\ns3,1,1,0\n"
  wide <- read.csv(text = csv)
  expect_error(beta_diversity(wide), "'site' contain text")
  expect_error(beta_diversity(wide), "row.names = 1")
  expect_warning(beta_diversity(read.csv(text = csv, row.names = 1)), "Only N = 3")
})

## ---- beta_diversity_test(group = ) shortcut ----

test_that("beta_diversity_test(group = ) equals the two-table form, with labels and betas", {
  set.seed(2)
  a <- matrix(rbinom(15 * 12, 1, 0.4), nrow = 15)
  b <- matrix(rbinom(15 * 12, 1, 0.2), nrow = 15)
  two <- beta_diversity_test(a, b)
  expect_equal(two$group_1, "x")
  expect_equal(two$group_2, "y")

  one <- beta_diversity_test(rbind(b, a), group = rep(c("tree", "open"), each = 15))
  expect_equal(one$group_1, "open")  # sorted order
  expect_equal(one$group_2, "tree")
  expect_equal(one$Z_beta_M, two$Z_beta_M)
  expect_equal(one$p_beta_M, two$p_beta_M)
  expect_equal(one$beta_M_1, beta_diversity(a)$beta_M)

  # Manual Eq. 10a recomputation
  s1 <- beta_diversity(a)
  s2 <- beta_diversity(b)
  z <- log(s1$beta_M / s2$beta_M) / sqrt(s1$se_ln_beta_M^2 + s2$se_ln_beta_M^2)
  expect_equal(two$Z_beta_M, z)
})

test_that("beta_diversity_test() needs exactly one of y or group, and exactly 2 groups", {
  set.seed(1)
  x <- matrix(rbinom(30 * 10, 1, 0.35), nrow = 30)
  expect_error(beta_diversity_test(x), "exactly one of `y` or `group`")
  expect_error(beta_diversity_test(x, x, group = rep(1:2, 15)),
               "exactly one of `y` or `group`")
  expect_error(beta_diversity_test(x, group = rep(c("a", "b", "c"), 10)),
               "exactly 2 groups, but `group` has 3: a, b, c")
})

test_that("beta_diversity_test(group = ) leaves unlabelled sites out, with a warning", {
  set.seed(1)
  x <- matrix(rbinom(24 * 10, 1, 0.35), nrow = 24)
  group <- rep(c("a", "b"), each = 12)
  group_na <- group
  group_na[c(1, 24)] <- c(NA, "")
  expect_warning(res <- beta_diversity_test(x, group = group_na),
                 "left out of the test: row 1, row 24")
  expected <- beta_diversity_test(x[2:12, ], x[13:23, ])
  expect_equal(res$Z_beta_M, expected$Z_beta_M)
})

#' Point estimates and CIs for a single site-by-species matrix (internal)
#'
#' Not exported. Contains the full Kiflawi & Spencer (2004) computation and
#' all of its input validation. \code{\link{beta_diversity}} is a thin
#' group-aware wrapper around this.
#'
#' @importFrom stats var qnorm pnorm
#' @noRd
.beta_diversity_one <- function(x, conf) {
  x <- as.matrix(x)

  if (!is.numeric(x) && !is.logical(x)) {
    stop("x must contain only numbers or TRUE/FALSE values (e.g. 0/1, ",
         "presence/absence, counts, or weights) -- not text. If your data ",
         "has text like \"yes\"/\"no\", convert it first with ",
         "as_site_by_species().")
  }
  if (anyNA(x)) {
    warning("x has some blank/missing values (NA); treating them as absence (0)")
    x[is.na(x)] <- 0
  }
  if (any(!is.finite(x))) {
    stop("x must not contain infinite values")
  }
  if (any(x < 0)) {
    stop("x must not contain negative values")
  }
  if (ncol(x) < 1) {
    stop("x must have at least 1 species (column)")
  }
  N <- nrow(x)
  if (N < 2) {
    stop("beta_diversity() requires at least 2 sites (rows); got ", N)
  }
  if (N < 10) {
    warning("N = ", N, " sites is below 10, the smallest sample size Kiflawi & ",
            "Spencer (2004) validated CI coverage for; results may be unreliable")
  }

  # Eq. 6a: first-order jackknife estimate of gamma richness. freq[j] = number
  # of sites species j occurs at; a "unique" species occurs at exactly one
  # site. u = total count of unique species.
  freq <- colSums(x > 0)
  gamma_obs <- sum(freq > 0)
  u <- sum(freq == 1)
  gamma_hat <- gamma_obs + u * (N - 1) / N

  # Eq. 6b: sampling variance of gamma_hat, from the per-site distribution of
  # unique-species counts. k_i = number of the u unique species found at site
  # i; by construction sum(k_i) == u.
  if (u == 0) {
    se_gamma_hat <- 0
  } else {
    k_i <- rowSums(x[, freq == 1, drop = FALSE] > 0)
    var_gamma_hat <- (N - 1) / N * (sum(k_i^2) - u^2 / N)
    se_gamma_hat <- sqrt(var_gamma_hat)
  }

  site_rich <- rowSums(x > 0)
  alpha <- mean(site_rich)
  # Sampling variance of alpha (mean local richness): algebraically equivalent
  # to Kiflawi & Spencer's Eq. 8 (species-pairwise covariance form), but
  # simplifies to the ordinary standard-error-of-the-mean of per-site richness.
  var_alpha <- var(site_rich) * (N - 1) / N^2
  se_alpha <- sqrt(var_alpha)

  if (gamma_hat <= alpha) {
    stop("beta_M is undefined: estimated gamma (", round(gamma_hat, 3),
         ") is not greater than alpha (", round(alpha, 3), "). This happens ",
         "when one or more species are present at every site. Consider using ",
         "beta_A, or removing ubiquitous species before calling beta_diversity().")
  }

  beta_M <- gamma_hat / alpha - 1
  beta_A <- gamma_hat - alpha

  # Delta-method variance of ln(beta_M), Eq. 7
  var_ln_beta_M <- (gamma_hat / (alpha * (gamma_hat - alpha)))^2 * var_alpha +
    (1 / (gamma_hat - alpha))^2 * se_gamma_hat^2
  se_ln_beta_M <- sqrt(var_ln_beta_M)

  se_beta_A <- sqrt(se_gamma_hat^2 + var_alpha)

  z <- qnorm(1 - (1 - conf) / 2)

  beta_M_ci <- beta_M * exp(c(-1, 1) * z * se_ln_beta_M)
  beta_A_ci <- beta_A + c(-1, 1) * z * se_beta_A

  data.frame(
    n_sites = N,
    gamma_obs = gamma_obs,
    gamma_hat = gamma_hat,
    se_gamma_hat = se_gamma_hat,
    alpha = alpha,
    se_alpha = se_alpha,
    beta_M = beta_M,
    se_ln_beta_M = se_ln_beta_M,
    beta_M_lwr = beta_M_ci[1],
    beta_M_upr = beta_M_ci[2],
    beta_A = beta_A,
    se_beta_A = se_beta_A,
    beta_A_lwr = beta_A_ci[1],
    beta_A_upr = beta_A_ci[2]
  )
}

#' Beta diversity point estimates and confidence intervals (Kiflawi & Spencer 2004)
#'
#' `r lifecycle::badge("experimental")`
#'
#' Computes Whittaker's multiplicative beta diversity (beta_M = gamma/alpha - 1)
#' and Lande's additive beta diversity (beta_A = gamma - alpha) from a
#' sites x species presence/absence matrix, along with analytical 95%
#' confidence intervals, following Kiflawi & Spencer (2004, Ecology 85:2895-2900).
#' If \code{group} is supplied, sites are split by group and every quantity
#' is computed separately per group.
#'
#' The true regional richness (gamma) is not observed directly; it is estimated
#' via the first-order jackknife richness estimator (Heltshe & Forrester 1983),
#' computed directly following Kiflawi & Spencer (2004) Eqs. 6a-6b. The point
#' estimate (\code{gamma_hat}) is numerically identical to
#' \code{vegan::specpool()}'s \code{jack1}; its standard error
#' (\code{se_gamma_hat}) follows the paper's own Eq. 6b derivation and differs
#' from \code{vegan::specpool()}'s \code{jack1.se}, which is based on a
#' different derivation (Smith & van Belle 1984).
#'
#' @param x A sites (rows) x species (columns) matrix or data.frame,
#'   numeric or logical. Values are treated as present if > 0, so either
#'   presence/absence (0/1, TRUE/FALSE) or raw abundance/count data may be
#'   supplied directly. \code{NA} values are treated as absence (with a
#'   warning). Errors if \code{x} is non-numeric, contains negative or
#'   non-finite values, has fewer than 2 sites or 1 species, or if the
#'   estimated gamma diversity is not greater than alpha (e.g. because a
#'   species is present at every site, making beta_M undefined). Warns if
#'   there are fewer than 10 sites, the smallest sample size Kiflawi &
#'   Spencer (2004) validated confidence-interval coverage for.
#' @param group Optional vector of length \code{nrow(x)} assigning each site
#'   to a group (e.g. treatment/plot type); results are then computed
#'   separately per group. Defaults to \code{NULL}, treating all sites as a
#'   single group. At most 5 distinct groups are supported (error above
#'   that). A group with fewer than 2 sites is skipped, with a warning,
#'   rather than failing the whole call. \code{NA} entries in \code{group}
#'   cause those sites to be dropped from the analysis, with a warning.
#' @param conf Confidence level for the intervals (default 0.95). Must be a
#'   single number strictly between 0 and 1; errors otherwise (e.g. if you
#'   pass a percentage like 95 instead of 0.95).
#'
#' @return A data.frame with one row per group (one row total if
#'   \code{group} is \code{NULL}), with columns: \code{group} (always
#'   \code{"All"} when ungrouped), \code{n_sites}, \code{gamma_obs},
#'   \code{gamma_hat}, \code{se_gamma_hat}, \code{alpha}, \code{se_alpha},
#'   \code{beta_M}, \code{se_ln_beta_M}, \code{beta_M_lwr}, \code{beta_M_upr},
#'   \code{beta_A}, \code{se_beta_A}, \code{beta_A_lwr}, \code{beta_A_upr}.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rbinom(8 * 10, 1, 0.35), nrow = 8,
#'              dimnames = list(paste0("site", 1:8), paste0("sp", 1:10)))
#' beta_diversity(x)
#'
#' set.seed(1)
#' x2 <- matrix(rbinom(20 * 10, 1, 0.35), nrow = 20)
#' treatment <- rep(c("control", "npk"), each = 10)
#' beta_diversity(x2, group = treatment)
#'
#' @family beta diversity functions
#' @export
beta_diversity <- function(x, group = NULL, conf = 0.95) {
  if (!is.numeric(conf) || length(conf) != 1 || is.na(conf)) {
    stop("conf must be a single number between 0 and 1 ",
         "(e.g. 0.95 for a 95% confidence interval); got ",
         if (length(conf) == 1) conf else paste0("a vector of length ", length(conf)))
  }
  if (conf <= 0 || conf >= 1) {
    if (conf >= 50 && conf <= 100) {
      stop("conf must be a proportion between 0 and 1, not a percentage. ",
           "Did you mean conf = ", conf / 100, "?")
    }
    stop("conf must be strictly between 0 and 1 ",
         "(e.g. 0.95 for a 95% confidence interval); got conf = ", conf)
  }

  x <- as.matrix(x)

  if (is.null(group)) {
    result <- .beta_diversity_one(x, conf)
    result$group <- "All"
    return(result[, c("group", setdiff(names(result), "group"))])
  }

  if (length(group) != nrow(x)) {
    stop("length(group) must equal nrow(x); got ", length(group),
         " and ", nrow(x))
  }

  if (anyNA(group)) {
    n_dropped <- sum(is.na(group))
    warning(n_dropped, " site(s) have a blank/missing group label and were ",
            "dropped from the analysis")
    x <- x[!is.na(group), , drop = FALSE]
    group <- group[!is.na(group)]
  }

  groups <- sort(unique(group))
  if (length(groups) > 5) {
    stop("beta_diversity() supports at most 5 groups; got ", length(groups))
  }

  results <- list()
  for (g in groups) {
    x_g <- x[group == g, , drop = FALSE]
    if (nrow(x_g) < 2) {
      warning("Group '", g, "' has fewer than 2 sites and was skipped")
      next
    }
    r <- .beta_diversity_one(x_g, conf)
    r$group <- as.character(g)
    results[[length(results) + 1]] <- r
  }

  if (length(results) == 0) {
    stop("No group had at least 2 sites; nothing to return")
  }

  out <- do.call(rbind, results)
  out[, c("group", setdiff(names(out), "group"))]
}

#' Two-sample hypothesis test for beta diversity (Kiflawi & Spencer 2004)
#'
#' `r lifecycle::badge("experimental")`
#'
#' Tests whether two independent site-sets differ in beta diversity, for both
#' the multiplicative (beta_M) and additive (beta_A) measures, following
#' Eqs. 10a-10b of Kiflawi & Spencer (2004).
#'
#' @param x,y Sites x species presence/absence matrices or data.frames for
#'   the two independent samples being compared.
#' @param conf Confidence level, passed through to the underlying
#'   \code{\link{beta_diversity}} calls (does not affect the test itself).
#'
#' @return A one-row data.frame with columns \code{Z_beta_M}, \code{p_beta_M},
#'   \code{Z_beta_A}, \code{p_beta_A} (two-tailed p-values).
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rbinom(8 * 10, 1, 0.35), nrow = 8)
#' y <- matrix(rbinom(8 * 10, 1, 0.35), nrow = 8)
#' beta_diversity_test(x, y)
#'
#' @family beta diversity functions
#' @export
beta_diversity_test <- function(x, y, conf = 0.95) {
  s1 <- beta_diversity(x, conf = conf)
  s2 <- beta_diversity(y, conf = conf)

  Z_beta_M <- log(s1$beta_M / s2$beta_M) /
    sqrt(s1$se_ln_beta_M^2 + s2$se_ln_beta_M^2)
  Z_beta_A <- (s1$beta_A - s2$beta_A) /
    sqrt(s1$se_beta_A^2 + s2$se_beta_A^2)

  data.frame(
    Z_beta_M = Z_beta_M,
    p_beta_M = 2 * pnorm(-abs(Z_beta_M)),
    Z_beta_A = Z_beta_A,
    p_beta_A = 2 * pnorm(-abs(Z_beta_A))
  )
}

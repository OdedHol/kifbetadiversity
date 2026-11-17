#' Signal that beta_M can't be estimated for a set of sites (internal)
#'
#' A classed error, so that a grouped run in \code{beta_diversity()} can skip
#' just the affected group instead of failing the whole call.
#' @noRd
.stop_no_turnover <- function(msg) {
  stop(structure(
    class = c("kifbeta_no_turnover", "error", "condition"),
    list(message = msg, call = NULL)
  ))
}

#' Turn blank / whitespace-only group labels into NA (internal)
#' @noRd
.clean_group <- function(group) {
  g <- trimws(as.character(group))
  g[!is.na(g) & g == ""] <- NA_character_
  g
}

#' Short, human-readable list of site names (internal)
#' @noRd
.list_sites <- function(sites, max_show = 10) {
  sites <- as.character(sites)
  if (length(sites) <= max_show) {
    return(paste(sites, collapse = ", "))
  }
  paste0(paste(sites[seq_len(max_show)], collapse = ", "),
         ", and ", length(sites) - max_show, " more")
}

#' Site names of the rows of x, falling back to row numbers (internal)
#' @noRd
.site_names <- function(x, rows) {
  rn <- rownames(x)
  if (is.null(rn)) paste0("row ", rows) else rn[rows]
}

#' Point estimates and CI for a single site-by-species matrix (internal)
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
    warning("Only N = ", N, " sites. Kiflawi & Spencer's (2004) simulations ",
            "found the beta_M confidence interval can be unreliable with few ",
            "sites (it was reasonably reliable from about 30 sites), so treat ",
            "this interval with caution.")
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
  # Eq. 8b: sampling variance of alpha (mean local richness). This is the
  # population-variance form, sum((s_i - alpha)^2) / N^2, i.e.
  # var(s_i) * (N - 1) / N^2 -- slightly smaller than the ordinary squared
  # standard error of the mean, var(s_i) / N.
  var_alpha <- var(site_rich) * (N - 1) / N^2
  se_alpha <- sqrt(var_alpha)

  # gamma_hat >= gamma_obs >= max(site richness) >= alpha always, so
  # gamma_hat == alpha only happens in the two cases below. beta_M is then 0
  # (or 0/0), and its log-scale interval (Eq. 7) is undefined.
  if (alpha == 0) {
    .stop_no_turnover(paste0(
      "beta_M can't be estimated: no species were recorded at any of the ",
      N, " sites."))
  }
  if (gamma_hat <= alpha) {
    .stop_no_turnover(paste0(
      "beta_M can't be estimated: every one of the ", N, " sites has exactly ",
      "the same ", gamma_obs, " species, so there is no turnover between ",
      "sites (beta_M = 0) and no confidence interval can be put around it."))
  }

  beta_M <- gamma_hat / alpha - 1

  # Delta-method variance of ln(beta_M), Eq. 7
  var_ln_beta_M <- (gamma_hat / (alpha * (gamma_hat - alpha)))^2 * var_alpha +
    (1 / (gamma_hat - alpha))^2 * se_gamma_hat^2
  se_ln_beta_M <- sqrt(var_ln_beta_M)

  z <- qnorm(1 - (1 - conf) / 2)
  beta_M_ci <- beta_M * exp(c(-1, 1) * z * se_ln_beta_M)

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
    # Conventional Whittaker beta, gamma / alpha = beta_M + 1. Same
    # information as beta_M; its interval is the beta_M interval shifted by 1.
    beta_W = beta_M + 1,
    beta_W_lwr = beta_M_ci[1] + 1,
    beta_W_upr = beta_M_ci[2] + 1
  )
}

#' Validate the conf argument (internal)
#' @noRd
.check_conf <- function(conf) {
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
  invisible(NULL)
}

#' Warn about columns that probably aren't species (internal)
#'
#' Flags columns whose name looks like an ID / summary / date column, whose
#' values equal the number (or total) of the other species at each site
#' (a richness or abundance column), or whose values are far out of scale
#' in an otherwise presence/absence table. Warns once, never errors.
#' @noRd
.check_species_columns <- function(m) {
  m <- as.matrix(m)
  p <- ncol(m)
  if (p < 2 || !(is.numeric(m) || is.logical(m))) {
    return(invisible(NULL))
  }
  storage.mode(m) <- "double"
  m[is.na(m)] <- 0
  labels <- colnames(m)
  if (is.null(labels)) labels <- paste0("column ", seq_len(p))
  reasons <- vector("list", p)

  if (!is.null(colnames(m))) {
    name_hit <- grepl(paste0("richness|abundance|total|^n_|_id$|^id$|num$|",
                             "^year|^date|^lat|^lon|^plot|^site"),
                      colnames(m), ignore.case = TRUE)
    for (j in which(name_hit)) {
      reasons[[j]] <- c(reasons[[j]], "its name looks like an ID, number or summary column")
    }
  }

  if (p >= 5) {
    present <- m > 0
    n_present <- rowSums(present)
    total <- rowSums(m)
    for (j in seq_len(p)) {
      col <- m[, j]
      if (max(col) <= 1) next
      if (all(col == n_present - present[, j])) {
        reasons[[j]] <- c(reasons[[j]], "it equals the number of other species at each site (a richness column?)")
      } else if (all(col == total - col)) {
        reasons[[j]] <- c(reasons[[j]], "it equals the total of the other columns at each site (an abundance column?)")
      }
    }

    is01 <- apply(m, 2, function(col) all(col %in% c(0, 1)))
    if (mean(is01) >= 0.9 && any(!is01)) {
      maxes <- apply(m, 2, max)
      for (j in which(!is01)) {
        if (maxes[j] >= 10 * max(1, stats::median(maxes[-j]))) {
          reasons[[j]] <- c(reasons[[j]], paste0(
            "its values go up to ", maxes[j], " while almost every other ",
            "column is 0/1"))
        }
      }
    }
  }

  flagged <- which(lengths(reasons) > 0)
  if (length(flagged) > 0) {
    detail <- vapply(flagged, function(j) {
      paste0("'", labels[j], "' (", paste(unique(reasons[[j]]), collapse = "; "), ")")
    }, character(1))
    warning("These columns are being treated as species, but may not be ",
            "species: ", paste(detail, collapse = ", "), ". Every column of ",
            "x counts as a species, so a column like this changes the ",
            "results. If they aren't species, leave them out, e.g. with ",
            "as_site_table(..., species = ) or by removing them from x.",
            call. = FALSE)
  }
  invisible(NULL)
}

#' Coerce x to a matrix, explaining text columns in a data frame (internal)
#' @noRd
.as_site_matrix <- function(x) {
  if (is.data.frame(x)) {
    is_text <- !vapply(x, function(col) is.numeric(col) || is.logical(col),
                       logical(1))
    if (any(is_text)) {
      text_cols <- names(x)[is_text]
      # Guess which text column holds site names (all different) and which
      # holds a habitat (few different values), to show a concrete fix.
      n_distinct <- vapply(x[text_cols], function(col) length(unique(col)),
                           integer(1))
      site_guess <- text_cols[n_distinct == nrow(x)][1]
      group_guess <- text_cols[n_distinct <= 5][1]
      hint <- if (!is.na(site_guess)) {
        paste0(" It looks like a table with one row per site, so prepare it ",
               "with as_site_table(), for example (a guess -- check the ",
               "column names): res <- as_site_table(x, site = \"", site_guess,
               "\"", if (!is.na(group_guess)) paste0(", group = \"", group_guess, "\""),
               ", species = <your species columns>), then ",
               "beta_diversity(res$x, group = res$group).")
      } else {
        ""
      }
      stop("x must contain only numbers or TRUE/FALSE values, but column(s) ",
           paste0("'", text_cols, "'", collapse = ", "), " contain text.",
           hint, " If your site names are the first column of a CSV, you can ",
           "also read it with read.csv(\"file.csv\", row.names = 1). If your ",
           "data is one row per record (site, species, ...), use ",
           "as_site_by_species() instead.")
    }
  }
  as.matrix(x)
}

#' Beta diversity point estimate and confidence interval (Kiflawi & Spencer 2004)
#'
#' `r lifecycle::badge("experimental")`
#'
#' Computes Whittaker's multiplicative beta diversity,
#' \eqn{\beta_M = \hat\gamma / \bar\alpha - 1}, from a sites x species
#' presence/absence matrix, with an analytical confidence interval, following
#' Kiflawi & Spencer (2004, Ecology 85:2895-2900). If \code{group} is
#' supplied, sites are split by group and every quantity is computed
#' separately per group.
#'
#' Observed regional richness undercounts the true species pool, because
#' rare species are missed by a finite survey. Gamma is therefore estimated
#' with the first-order jackknife (Heltshe & Forrester 1983), following
#' Kiflawi & Spencer (2004) Eqs. 6a-6b, and its uncertainty is carried into
#' the interval for beta_M (Eq. 7, delta method on \eqn{\ln\beta_M}). The
#' point estimate (\code{gamma_hat}) is numerically identical to
#' \code{vegan::specpool()}'s \code{jack1}; its standard error
#' (\code{se_gamma_hat}) follows the paper's own Eq. 6b and differs from
#' \code{vegan::specpool()}'s \code{jack1.se}, which uses a different
#' derivation (Smith & van Belle 1984).
#'
#' @param x A sites (rows) x species (columns) matrix or data.frame,
#'   numeric or logical. Values are treated as present if > 0, so either
#'   presence/absence (0/1, TRUE/FALSE) or raw abundance/count data may be
#'   supplied directly. \code{NA} values are treated as absence (with a
#'   warning). Errors if \code{x} contains text (naming the column), negative
#'   or non-finite values, or has fewer than 2 sites or 1 species. Warns if
#'   there are fewer than 10 sites, or if some columns look like they aren't
#'   species (e.g. a richness, ID or total column). Every column counts as a
#'   species, so prepare a table that has extra columns with
#'   \code{\link{as_site_table}}.
#' @param group Optional vector of length \code{nrow(x)} assigning each site
#'   to a group (e.g. habitat or treatment); results are then computed
#'   separately per group. Defaults to \code{NULL}, treating all sites as a
#'   single group. At most 5 labelled groups are supported. Sites with a
#'   blank or missing (\code{NA}) label are not dropped: they are analysed
#'   together as their own group, reported as \code{NA} in the last row,
#'   with a warning naming them, so you can fix the labels. A group with
#'   fewer than 2 sites, or one where beta_M can't be estimated (see
#'   Value), is skipped with a warning while the other groups are returned.
#' @param conf Confidence level for the interval (default 0.95). Must be a
#'   single number strictly between 0 and 1; errors otherwise (e.g. if you
#'   pass a percentage like 95 instead of 0.95).
#'
#' @return A data.frame with one row per group (one row if \code{group} is
#'   \code{NULL}) and columns:
#'   \describe{
#'     \item{group}{Group label (\code{"All"} when ungrouped; \code{NA} for
#'       unlabelled sites).}
#'     \item{n_sites}{Number of sites, N.}
#'     \item{gamma_obs}{Species observed across all sites.}
#'     \item{gamma_hat}{Jackknife estimate of the true regional richness
#'       (Eq. 6a): observed species plus an estimate of those missed.}
#'     \item{se_gamma_hat}{Standard error of \code{gamma_hat} (Eq. 6b).}
#'     \item{alpha}{Mean number of species per site.}
#'     \item{se_alpha}{Standard error of \code{alpha} (Eq. 8b).}
#'     \item{beta_M}{Whittaker's multiplicative beta, \code{gamma_hat / alpha - 1}.
#'       It can be read as the odds that a species of average commonness is
#'       absent from a randomly chosen site; 0 means every site has every
#'       species.}
#'     \item{se_ln_beta_M}{Standard error of \code{log(beta_M)} (Eq. 7).}
#'     \item{beta_M_lwr, beta_M_upr}{Confidence interval for \code{beta_M}.
#'       It is built on the log scale, so it is not symmetric around
#'       \code{beta_M}.}
#'     \item{beta_W, beta_W_lwr, beta_W_upr}{Conventional Whittaker beta,
#'       \code{gamma_hat / alpha} (= \code{beta_M + 1}), and its interval
#'       (the \code{beta_M} interval + 1). Same information as
#'       \code{beta_M}; reported for comparison with studies that use
#'       gamma / alpha.}
#'   }
#'   beta_M can't be estimated, and an error is raised (or the group is
#'   skipped, in a grouped call), when every site has exactly the same
#'   species or when no species were recorded at all.
#'
#' @references Kiflawi, M., & Spencer, M. (2004). Confidence intervals and
#'   hypothesis testing for beta diversity. \emph{Ecology}, 85(10), 2895-2900.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30,
#'             dimnames = list(paste0("site", 1:30), paste0("sp", 1:15)))
#' beta_diversity(x)
#'
#' habitat <- rep(c("open", "tree"), each = 15)
#' beta_diversity(x, group = habitat)
#'
#' @family beta diversity functions
#' @export
beta_diversity <- function(x, group = NULL, conf = 0.95) {
  .beta_diversity(x, group = group, conf = conf, check_columns = TRUE)
}

#' beta_diversity() body, with the column check optional (internal)
#' @noRd
.beta_diversity <- function(x, group = NULL, conf = 0.95, check_columns = TRUE) {
  .check_conf(conf)
  x <- .as_site_matrix(x)
  if (check_columns) .check_species_columns(x)

  if (is.null(group)) {
    result <- .beta_diversity_one(x, conf)
    result$group <- "All"
    return(result[, c("group", setdiff(names(result), "group"))])
  }

  if (length(group) != nrow(x)) {
    stop("group must give one label per site: length(group) is ",
         length(group), " but x has ", nrow(x), " sites (rows)")
  }

  group <- .clean_group(group)
  unlabelled <- is.na(group)
  if (any(unlabelled)) {
    warning(sum(unlabelled), " site(s) have a blank or missing group label: ",
            .list_sites(.site_names(x, which(unlabelled))), ". They are ",
            "analysed together as a separate group, shown as NA in the ",
            "results. Fix their labels in your data, or remove those sites, ",
            "before reporting results.")
  }

  labelled_groups <- sort(unique(group[!unlabelled]))
  if (length(labelled_groups) > 5) {
    stop("beta_diversity() supports at most 5 groups; got ",
         length(labelled_groups), ": ",
         paste(labelled_groups, collapse = ", "))
  }
  groups <- if (any(unlabelled)) c(labelled_groups, NA) else labelled_groups

  results <- list()
  for (g in groups) {
    in_g <- if (is.na(g)) unlabelled else (!unlabelled & group == g)
    label <- if (is.na(g)) "NA (unlabelled sites)" else paste0("'", g, "'")
    x_g <- x[in_g, , drop = FALSE]
    if (nrow(x_g) < 2) {
      warning("Group ", label, " has fewer than 2 sites and was skipped")
      next
    }
    r <- tryCatch(
      .beta_diversity_one(x_g, conf),
      kifbeta_no_turnover = function(e) {
        warning("Group ", label, " was skipped. ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    if (is.null(r)) next
    r$group <- if (is.na(g)) NA_character_ else as.character(g)
    results[[length(results) + 1]] <- r
  }

  if (length(results) == 0) {
    stop("No group could be analysed (each had fewer than 2 sites, or no ",
         "turnover between its sites); nothing to return. See the warnings ",
         "above for each group.")
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out[, c("group", setdiff(names(out), "group"))]
}

#' Two-sample test: do two habitats differ in beta diversity? (Kiflawi & Spencer 2004)
#'
#' `r lifecycle::badge("experimental")`
#'
#' Tests whether two independent sets of sites (e.g. two habitats) differ in
#' Whittaker's multiplicative beta diversity, using the Z-test of Kiflawi &
#' Spencer (2004), Eq. 10a:
#' \deqn{Z = \ln(\beta_{M,1} / \beta_{M,2}) /
#'   \sqrt{Var[\ln\beta_{M,1}] + Var[\ln\beta_{M,2}]}}
#'
#' Give the two samples either as two separate tables (\code{x} and
#' \code{y}), or as one table plus a \code{group} vector with exactly two
#' labels.
#'
#' @section Assumptions and limitations:
#' \itemize{
#'   \item \strong{Independent, randomly sampled sites.} The standard errors
#'     assume every site is an independent random sample from its habitat.
#'     Spatially clumped designs (e.g. several traps under the same tree),
#'     paired designs (e.g. a tree trap and an open trap at the same
#'     location) or block designs carry less independent information than
#'     their number of sites suggests. The standard errors are then too
#'     small and Z too large, so a "significant" result may not be. Consider
#'     using a coarser, more independent unit (e.g. pool the traps of one
#'     tree into one site).
#'   \item \strong{The two samples must be independent of each other.}
#'   \item \strong{Comparable sampling effort.} The jackknife reduces the
#'     effect of sampling effort on estimated gamma, but does not remove it,
#'     and the largest possible beta_M grows with the number of sites (with
#'     observed gamma, beta_M can't exceed N - 1). Compare habitats sampled with a similar number of sites and
#'     similar effort per site.
#'   \item \strong{Enough sites.} In the paper's simulations the beta_M
#'     interval was reasonably reliable from about 30 sites per sample, and
#'     could be poor with 10.
#'   \item \strong{Include empty sites.} Sites surveyed where nothing was
#'     found must be included (see the \code{all_sites} argument of
#'     \code{\link{as_site_by_species}}); leaving them out makes beta_M too
#'     low.
#' }
#'
#' @param x A sites x species matrix or data.frame (see
#'   \code{\link{beta_diversity}}). Either the first sample (with \code{y}),
#'   or all sites of both samples (with \code{group}).
#' @param y A sites x species matrix or data.frame for the second sample.
#'   Give either \code{y} or \code{group}, not both.
#' @param conf Confidence level, passed through to the underlying
#'   \code{\link{beta_diversity}} calls (does not affect the test itself).
#' @param group A vector of length \code{nrow(x)} with exactly two distinct
#'   labels, splitting the sites of \code{x} into the two samples. Sites with
#'   a blank or missing label are left out of the test, with a warning.
#'
#' @return A one-row data.frame with columns \code{group_1}, \code{group_2}
#'   (the two sample labels: \code{"x"} and \code{"y"}, or the two
#'   \code{group} labels in sorted order), \code{beta_M_1}, \code{beta_M_2},
#'   \code{Z_beta_M} (positive when \code{group_1} has the higher beta_M) and
#'   \code{p_beta_M} (two-tailed p-value).
#'
#' @references Kiflawi, M., & Spencer, M. (2004). Confidence intervals and
#'   hypothesis testing for beta diversity. \emph{Ecology}, 85(10), 2895-2900.
#'
#' @examples
#' set.seed(1)
#' open <- matrix(rbinom(30 * 15, 1, 0.45), nrow = 30)
#' tree <- matrix(rbinom(30 * 15, 1, 0.20), nrow = 30)
#' beta_diversity_test(open, tree)
#'
#' # Or: one table plus a group vector
#' both <- rbind(open, tree)
#' habitat <- rep(c("open", "tree"), each = 30)
#' beta_diversity_test(both, group = habitat)
#'
#' @family beta diversity functions
#' @export
beta_diversity_test <- function(x, y = NULL, conf = 0.95, group = NULL) {
  if (is.null(y) == is.null(group)) {
    stop("Give the two samples either as two tables, ",
         "beta_diversity_test(x, y), or as one table plus a group vector ",
         "with two labels, beta_diversity_test(x, group = habitat) -- ",
         "exactly one of `y` or `group`.")
  }

  if (!is.null(y)) {
    s1 <- beta_diversity(x, conf = conf)
    s2 <- beta_diversity(y, conf = conf)
    labels <- c("x", "y")
  } else {
    x <- .as_site_matrix(x)
    if (length(group) != nrow(x)) {
      stop("group must give one label per site: length(group) is ",
           length(group), " but x has ", nrow(x), " sites (rows)")
    }
    group <- .clean_group(group)
    unlabelled <- is.na(group)
    if (any(unlabelled)) {
      warning(sum(unlabelled), " site(s) have a blank or missing group ",
              "label and were left out of the test: ",
              .list_sites(.site_names(x, which(unlabelled))), ". Fix their ",
              "labels in your data to include them.")
      x <- x[!unlabelled, , drop = FALSE]
      group <- group[!unlabelled]
    }
    labels <- sort(unique(group))
    if (length(labels) != 2) {
      stop("beta_diversity_test() compares exactly 2 groups, but `group` ",
           "has ", length(labels), ": ", paste(labels, collapse = ", "),
           ". To compare two of them, subset your data first, e.g. ",
           "keep <- group %in% c(\"A\", \"B\").")
    }
    .check_species_columns(x)
    s1 <- .beta_diversity(x[group == labels[1], , drop = FALSE], conf = conf,
                          check_columns = FALSE)
    s2 <- .beta_diversity(x[group == labels[2], , drop = FALSE], conf = conf,
                          check_columns = FALSE)
  }

  # Eq. 10a
  Z_beta_M <- log(s1$beta_M / s2$beta_M) /
    sqrt(s1$se_ln_beta_M^2 + s2$se_ln_beta_M^2)

  data.frame(
    group_1 = labels[1],
    group_2 = labels[2],
    beta_M_1 = s1$beta_M,
    beta_M_2 = s2$beta_M,
    Z_beta_M = Z_beta_M,
    p_beta_M = 2 * pnorm(-abs(Z_beta_M))
  )
}

#' Check a column-name argument is a single non-NA string (internal)
#' @noRd
.check_col_arg <- function(val, arg_name) {
  if (is.null(val)) {
    return(invisible(NULL))
  }
  if (!is.character(val) || length(val) != 1 || is.na(val)) {
    stop("`", arg_name, "` must be the name of a column in your data, ",
         "given as a single piece of text (e.g. \"site\")")
  }
  invisible(NULL)
}

#' Warn about ids that differ only by case (internal)
#' @noRd
.warn_case_collisions <- function(ids, what) {
  uids <- unique(ids)
  if (length(uids) == length(unique(tolower(uids)))) {
    return(invisible(NULL))
  }
  collisions <- split(uids, tolower(uids))
  collisions <- collisions[vapply(collisions, length, integer(1)) > 1]
  detail <- vapply(collisions, paste, character(1), collapse = "/")
  warning(
    "Some ", what, " names differ only by case and are being kept as ",
    "distinct ", what, "s: ", paste(detail, collapse = "; "),
    ". If these should be the same ", what, ", standardize their case ",
    "before calling as_site_by_species()."
  )
  invisible(NULL)
}

#' Resolve a raw value column into a logical presence vector (internal)
#' @noRd
.resolve_presence <- function(data, value) {
  n <- nrow(data)
  if (is.null(value)) {
    return(rep(TRUE, n))
  }

  raw <- data[[value]]

  if (is.logical(raw)) {
    presence <- raw
  } else if (is.numeric(raw)) {
    presence <- raw > 0
  } else {
    chr <- trimws(tolower(as.character(raw)))
    true_tokens <- c("y", "yes", "true", "t", "present", "p")
    false_tokens <- c("n", "no", "false", "f", "absent", "a")
    presence <- rep(NA, length(chr))
    presence[chr %in% true_tokens] <- TRUE
    presence[chr %in% false_tokens] <- FALSE

    # Anything not matched by the tokens above may still be a plain number
    # (e.g. a weight/count stored as text); threshold those the same way
    # numeric columns are thresholded.
    unresolved <- is.na(presence) & !is.na(raw)
    if (any(unresolved)) {
      numeric_chr <- suppressWarnings(as.numeric(chr[unresolved]))
      presence[unresolved] <- numeric_chr > 0
    }

    unparsed <- unique(raw[!is.na(raw) & is.na(presence)])
    if (length(unparsed) > 0) {
      stop("Unrecognized value(s) in `", value, "`: ",
           paste(unparsed, collapse = ", "),
           ". Supported: TRUE/FALSE, positive counts/weights, or Y/N, ",
           "yes/no, present/absent (case-insensitive).")
    }
  }

  if (anyNA(presence)) {
    warning(sum(is.na(presence)), " row(s) have a blank/missing value (NA) ",
            "in '", value, "'; treating them as absence")
    presence[is.na(presence)] <- FALSE
  }

  presence
}

#' Add zero-detection sites / validate all_sites (internal)
#' @noRd
.apply_all_sites <- function(mat, all_sites) {
  if (is.null(all_sites)) {
    warning(
      "No `all_sites` list was given, so only sites where something was ",
      "actually found are included. If you surveyed a site and found ",
      "nothing there, it won't show up here -- which can make it look like ",
      "sites differ more than they really do. Pass all_sites = <every site ",
      "you surveyed> to include zero-detection sites too."
    )
    return(mat)
  }

  if (anyNA(all_sites)) {
    stop("`all_sites` can't contain blank or missing site names")
  }
  all_sites <- trimws(as.character(all_sites))
  if (anyDuplicated(all_sites)) {
    warning("`all_sites` lists the same site more than once (duplicate ",
            "site name); the repeats were removed")
    all_sites <- unique(all_sites)
  }

  observed <- rownames(mat)
  not_in_all_sites <- setdiff(observed, all_sites)
  if (length(not_in_all_sites) > 0) {
    stop("Site(s) appear in `data` but not in `all_sites`: ",
         paste(not_in_all_sites, collapse = ", "),
         ". `all_sites` must include every surveyed site.")
  }

  zero_sites <- setdiff(all_sites, observed)
  if (length(zero_sites) > 0) {
    zero_rows <- matrix(0, nrow = length(zero_sites), ncol = ncol(mat),
                         dimnames = list(zero_sites, colnames(mat)))
    mat <- rbind(mat, zero_rows)
  }

  mat[all_sites, , drop = FALSE]
}

#' Resolve the per-site group/category vector (internal)
#' @noRd
.resolve_group <- function(data, group, site_ids, mat_rownames) {
  if (is.null(group)) {
    return(NULL)
  }

  grp_raw <- trimws(as.character(data[[group]]))
  site_group <- unique(data.frame(site = site_ids, grp = grp_raw,
                                   stringsAsFactors = FALSE))

  dup_sites <- unique(site_group$site[duplicated(site_group$site)])
  if (length(dup_sites) > 0) {
    detail <- vapply(dup_sites, function(s) {
      vals <- unique(site_group$grp[site_group$site == s])
      paste0(s, " (", paste(vals, collapse = " vs "), ")")
    }, character(1))
    stop("Site(s) have inconsistent group/category labels across rows: ",
         paste(detail, collapse = "; "))
  }

  lookup <- stats::setNames(site_group$grp, site_group$site)
  group_vec <- unname(lookup[mat_rownames])

  if (anyNA(group_vec)) {
    warning(sum(is.na(group_vec)), " site(s) have no records in your data, ",
            "so they have no group/category label (shown as missing/NA)")
  }

  group_vec
}

#' Reshape raw long-format survey records into a site x species matrix
#'
#' `r lifecycle::badge("experimental")`
#'
#' Takes one-row-per-observation field data (a site ID, a species ID, and
#' optionally a detection value and/or a category label) and reshapes it into
#' the sites (rows) x species (columns) matrix expected by
#' \code{\link{beta_diversity}}, along with an aligned per-site group vector
#' when a category column is supplied.
#'
#' @param x A long-format data.frame, one row per observation.
#' @param site Column name (single string) in \code{x} identifying the
#'   site of each record. Required. Rows with \code{NA} site error.
#' @param species Column name (single string) in \code{x} identifying the
#'   species of each record. Required. Rows with \code{NA} species error.
#' @param value Optional column name (single string) giving the detection
#'   value for each record: logical, numeric (counts or weights, thresholded
#'   at > 0, matching \code{\link{beta_diversity}}'s own convention), or
#'   character/factor (parsed case-insensitively as
#'   TRUE/FALSE/Y/N/yes/no/present/absent, or as a plain number if not one of
#'   those tokens). \code{NA} values are treated as absence, with a warning.
#'   Unparseable non-NA values error, naming the offending value(s). Defaults
#'   to \code{NULL}, meaning every row in \code{x} is itself a detection
#'   (presence-only records).
#' @param group Optional column name (single string) giving a per-site
#'   category/habitat label (e.g. "under tree" vs "open field"), passed
#'   through as the \code{group} vector for \code{\link{beta_diversity}}. A
#'   site with inconsistent labels across its rows errors, naming the site
#'   and the conflicting values. Defaults to \code{NULL}.
#' @param all_sites Optional character vector of every surveyed site ID,
#'   including sites with zero detections. Since presence-only field data has
#'   no way to represent "surveyed but nothing seen," those sites are added
#'   as all-zero rows. If omitted, a warning explains that zero-detection
#'   sites (if any) will be silently missing, deflating N and biasing
#'   beta-diversity estimates. A site present in \code{x} but missing from
#'   \code{all_sites} errors, naming the site. Defaults to \code{NULL}.
#'
#' @details
#' Site and species IDs are automatically trimmed of leading/trailing
#' whitespace. Names differing only by case (e.g. "SiteA" vs "sitea") are
#' kept as distinct sites/species, with a warning. Duplicate (site, species)
#' records (e.g. repeat visits) are collapsed: a species is treated as
#' present at a site if detected on any of its records there, reported via
#' \code{message()} since this is expected, benign behavior for raw field
#' data rather than a problem with the input.
#'
#' @return A list with two elements: \code{x}, the sites x species matrix
#'   ready to pass as \code{beta_diversity(x = ...)}; and \code{group}, the
#'   per-site category vector aligned to \code{rownames(x)} (or \code{NULL}
#'   if \code{group} was not supplied).
#'
#' @examples
#' records <- data.frame(
#'   site = c("s1", "s1", "s2", "s2", "s3"),
#'   species = c("sp1", "sp2", "sp1", "sp3", "sp2"),
#'   seen = c("Y", "Y", "yes", "N", "Y"),
#'   habitat = c("tree", "tree", "open", "open", "tree")
#' )
#' res <- as_site_by_species(records, site = "site", species = "species",
#'                            value = "seen", group = "habitat")
#' res$x
#' res$group
#' beta_diversity(res$x, group = res$group)
#'
#' @family beta diversity functions
#' @export
as_site_by_species <- function(x, site, species, value = NULL,
                                group = NULL, all_sites = NULL) {
  if (!is.data.frame(x) || nrow(x) == 0) {
    stop("Your data table is empty (0 rows) -- as_site_by_species() needs ",
         "at least one record to reshape")
  }

  .check_col_arg(site, "site")
  .check_col_arg(species, "species")
  .check_col_arg(value, "value")
  .check_col_arg(group, "group")

  needed <- c(site, species, value, group)
  missing_cols <- setdiff(needed, names(x))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in `x`: ", paste(missing_cols, collapse = ", "),
         ". Available columns: ", paste(names(x), collapse = ", "))
  }

  site_ids <- x[[site]]
  species_ids <- x[[species]]

  if (anyNA(site_ids)) {
    stop(sum(is.na(site_ids)), " row(s) are missing a site name in the '",
         site, "' column; every record needs to say which site it's from")
  }
  if (anyNA(species_ids)) {
    stop(sum(is.na(species_ids)), " row(s) are missing a species name in ",
         "the '", species, "' column; every record needs to say which ",
         "species it's about")
  }

  site_ids <- trimws(as.character(site_ids))
  species_ids <- trimws(as.character(species_ids))

  .warn_case_collisions(site_ids, "site")
  .warn_case_collisions(species_ids, "species")

  presence <- .resolve_presence(x, value)

  n_dup <- sum(duplicated(data.frame(site_ids, species_ids)))
  if (n_dup > 0) {
    message(n_dup, " duplicate (site, species) record(s) found (e.g. repeat ",
            "visits); a species is treated as present at a site if detected ",
            "on any of its records there")
  }

  wide <- tidyr::pivot_wider(
    data.frame(site = site_ids, species = species_ids,
               presence = as.numeric(presence)),
    id_cols = site, names_from = species, values_from = presence,
    values_fn = max, values_fill = 0
  )
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$site
  mat <- mat[order(rownames(mat)), order(colnames(mat)), drop = FALSE]
  storage.mode(mat) <- "double"

  mat <- .apply_all_sites(mat, all_sites)

  group_vec <- .resolve_group(x, group, site_ids, rownames(mat))

  list(x = mat, group = group_vec)
}

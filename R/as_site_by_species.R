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
    # A blank cell (read.csv() reads it as "" in a text column) means the
    # same as a missing value.
    blank <- is.na(chr) | chr == ""
    true_tokens <- c("y", "yes", "true", "t", "present", "p")
    false_tokens <- c("n", "no", "false", "f", "absent", "a")
    presence <- rep(NA, length(chr))
    presence[chr %in% true_tokens] <- TRUE
    presence[chr %in% false_tokens] <- FALSE

    # Anything not matched by the tokens above may still be a plain number
    # (e.g. a weight/count stored as text); threshold those the same way
    # numeric columns are thresholded.
    unresolved <- is.na(presence) & !blank
    if (any(unresolved)) {
      numeric_chr <- suppressWarnings(as.numeric(chr[unresolved]))
      presence[unresolved] <- numeric_chr > 0
    }

    unparsed <- unique(as.character(raw)[!blank & is.na(presence)])
    if (length(unparsed) > 0) {
      stop("Unrecognized value(s) in `", value, "`: ",
           paste0("\"", unparsed, "\"", collapse = ", "),
           ". Supported: TRUE/FALSE, numbers (counts/weights; above 0 means ",
           "present), or Y/N, yes/no, T/F, present/absent, P/A ",
           "(case-insensitive). Blank cells count as absence.")
    }
  }

  if (anyNA(presence)) {
    warning(sum(is.na(presence)), " row(s) have a blank/missing value in '",
            value, "'; treating them as absence (not found)")
    presence[is.na(presence)] <- FALSE
  }

  presence
}

#' Split an all_sites argument into site names and (optional) group labels
#' (internal)
#' @noRd
.parse_all_sites <- function(all_sites, site, group) {
  if (is.null(all_sites)) {
    return(list(sites = NULL, groups = NULL))
  }
  if (!is.data.frame(all_sites)) {
    return(list(sites = all_sites, groups = NULL))
  }
  if (!site %in% names(all_sites)) {
    stop("`all_sites` is a table, so it needs a column named '", site,
         "' (the same name as your `site` column) listing every surveyed ",
         "site. Its columns are: ", paste(names(all_sites), collapse = ", "))
  }
  sites <- all_sites[[site]]
  if (is.null(group)) {
    return(list(sites = sites, groups = NULL))
  }
  if (!group %in% names(all_sites)) {
    stop("`all_sites` is a table and you used `group`, so it needs a column ",
         "named '", group, "' giving each site's group, so that sites ",
         "where nothing was found still get a group. Its columns are: ",
         paste(names(all_sites), collapse = ", "))
  }
  if (anyNA(sites)) {
    stop("`all_sites` can't contain blank or missing site names")
  }
  tab <- unique(data.frame(site = trimws(as.character(sites)),
                           grp = .clean_group(all_sites[[group]]),
                           stringsAsFactors = FALSE))
  conflict <- unique(tab$site[duplicated(tab$site)])
  if (length(conflict) > 0) {
    stop("`all_sites` gives more than one '", group, "' for site(s): ",
         .list_sites(conflict), ". Each site needs a single group.")
  }
  list(sites = sites, groups = stats::setNames(tab$grp, tab$site))
}

#' Add zero-detection sites / validate all_sites (internal)
#' @noRd
.apply_all_sites <- function(mat, all_sites) {
  if (is.null(all_sites)) {
    warning(
      "Calculating beta diversity without an explicit list of all surveyed ",
      "sites (`all_sites`). If any sites yielded zero species, they are ",
      "missing from the data, which will make beta diversity estimates too ",
      "low (an empty site lowers alpha, the average species per site). Pass ",
      "all_sites = <every site you surveyed> to include them."
    )
    return(mat)
  }

  if (anyNA(all_sites)) {
    stop("`all_sites` can't contain blank or missing site names")
  }
  all_sites <- trimws(as.character(all_sites))
  if (any(all_sites == "")) {
    stop("`all_sites` can't contain blank or missing site names")
  }
  if (anyDuplicated(all_sites)) {
    warning("`all_sites` lists the same site more than once (duplicate ",
            "site name); the repeats were removed")
    all_sites <- unique(all_sites)
  }

  observed <- rownames(mat)
  not_in_all_sites <- setdiff(observed, all_sites)
  if (length(not_in_all_sites) > 0) {
    stop("Site(s) appear in your records but not in `all_sites`: ",
         .list_sites(not_in_all_sites),
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
.resolve_group <- function(data, group, site_ids, mat_rownames,
                           table_groups = NULL) {
  if (is.null(group)) {
    return(NULL)
  }

  grp_raw <- .clean_group(data[[group]])
  site_group <- unique(data.frame(site = site_ids, grp = grp_raw,
                                   stringsAsFactors = FALSE))
  # A row with a blank label doesn't contradict another row of the same site
  # that has one.
  labelled <- site_group[!is.na(site_group$grp), , drop = FALSE]

  dup_sites <- unique(labelled$site[duplicated(labelled$site)])
  if (length(dup_sites) > 0) {
    detail <- vapply(dup_sites, function(s) {
      vals <- unique(labelled$grp[labelled$site == s])
      paste0(s, " (", paste(vals, collapse = " vs "), ")")
    }, character(1))
    stop("Site(s) have inconsistent group/category labels across rows: ",
         paste(detail, collapse = "; "))
  }

  lookup <- stats::setNames(labelled$grp, labelled$site)
  group_vec <- unname(lookup[mat_rownames])

  if (!is.null(table_groups)) {
    from_table <- unname(table_groups[mat_rownames])
    clash <- !is.na(group_vec) & !is.na(from_table) & group_vec != from_table
    if (any(clash)) {
      detail <- paste0(mat_rownames[clash], " (records: ", group_vec[clash],
                       " vs all_sites: ", from_table[clash], ")")
      stop("Your records and `all_sites` give different groups for ",
           "site(s): ", paste(detail, collapse = "; "))
    }
    group_vec[is.na(group_vec)] <- from_table[is.na(group_vec)]
  }

  if (anyNA(group_vec)) {
    missing_sites <- mat_rownames[is.na(group_vec)]
    warning(
      length(missing_sites), " site(s) have no group/category label: ",
      .list_sites(missing_sites), ". Their group is NA, and ",
      "beta_diversity() will report them as a separate NA group. ",
      "This happens when the label cell is blank, or when a site has no ",
      "records (nothing was found there). For sites with no records, give ",
      "`all_sites` as a table with a '", group, "' column so they get a group."
    )
  }

  group_vec
}

#' Reshape raw long-format survey records into a site x species matrix
#'
#' `r lifecycle::badge("experimental")`
#'
#' Takes one-row-per-record field data (a site ID, a species ID, and
#' optionally a detection value and/or a habitat label) and reshapes it into
#' the sites (rows) x species (columns) presence/absence matrix expected by
#' \code{\link{beta_diversity}}, along with an aligned per-site group vector
#' when a group column is supplied.
#'
#' @param x A long-format data.frame, one row per record (e.g. as read with
#'   \code{read.csv()}).
#' @param site Name of the column in \code{x} saying where each record is
#'   from (e.g. \code{"trap"}). Required. Rows with a blank or missing site
#'   are an error, naming the rows.
#' @param species Name of the column in \code{x} saying which species each
#'   record is (e.g. \code{"species"}, or a morphospecies code). Required.
#'   Rows with a blank or missing species are an error, naming the rows,
#'   since the package can't tell which species they were: fix or remove
#'   them first.
#' @param value Optional name of a column saying whether the species was
#'   found: logical, numeric (counts or weights; above 0 means present), or
#'   text (case-insensitive TRUE/FALSE, T/F, Y/N, yes/no, present/absent,
#'   P/A, or a number written as text). Blank or missing cells count as
#'   absence, with a warning. Other values are an error naming them.
#'   Defaults to \code{NULL}, meaning every row is a detection
#'   (presence-only records).
#' @param group Optional name of a column giving each site's habitat or
#'   treatment (e.g. \code{"habitat"}), returned as the \code{group} vector
#'   for \code{\link{beta_diversity}}. A site with two different labels
#'   across its rows is an error, naming the site and the labels. Blank
#'   labels count as missing. Defaults to \code{NULL}.
#' @param all_sites Every site you surveyed, including sites where nothing
#'   was found. Presence-only records can't show "surveyed, but nothing
#'   found", so without this those sites are missing, which makes beta
#'   diversity look too low. Either a vector of site names, or (when you use
#'   \code{group}) a data frame with a column named like your \code{site}
#'   column and one named like your \code{group} column, so that empty sites
#'   also get their habitat. Extra columns are ignored. A site in \code{x}
#'   but missing from \code{all_sites} is an error. If omitted, a warning
#'   explains the risk. Defaults to \code{NULL}.
#'
#' @details
#' Site and species IDs are trimmed of leading/trailing whitespace. Names
#' differing only by case (e.g. "SiteA" vs "sitea") are kept as distinct
#' sites/species, with a warning. Duplicate (site, species) records (e.g.
#' repeat visits, or one row per individual) are collapsed: a species is
#' present at a site if detected on any of its records there, reported via
#' \code{message()} since this is expected for raw field data.
#'
#' @return A list with two elements: \code{x}, the sites x species 0/1
#'   matrix ready to pass as \code{beta_diversity(x = ...)}; and
#'   \code{group}, the per-site group vector aligned to \code{rownames(x)}
#'   (\code{NULL} if \code{group} was not supplied; \code{NA} for a site
#'   with no label).
#'
#' @examples
#' records <- data.frame(
#'   site    = c("s1", "s1", "s2", "s2", "s3"),
#'   species = c("sp1", "sp2", "sp1", "sp3", "sp2"),
#'   seen    = c("Y", "Y", "yes", "N", "Y"),
#'   habitat = c("tree", "tree", "open", "open", "tree")
#' )
#' # s4 was surveyed but nothing was found there
#' sites <- data.frame(site = c("s1", "s2", "s3", "s4"),
#'                     habitat = c("tree", "open", "tree", "open"))
#' res <- as_site_by_species(records, site = "site", species = "species",
#'                           value = "seen", group = "habitat",
#'                           all_sites = sites)
#' res$x
#' res$group
#'
#' # A full example with enough sites per habitat:
#' path <- system.file("extdata", "example_survey.csv",
#'                     package = "kifbetadiversity")
#' head(read.csv(path))
#'
#' @family beta diversity functions
#' @export
as_site_by_species <- function(x, site, species, value = NULL,
                                group = NULL, all_sites = NULL) {
  if (!is.data.frame(x)) {
    stop("as_site_by_species() needs a data frame with one row per record ",
         "(e.g. from read.csv()), but got a ", class(x)[1], ". If your data ",
         "is already a sites x species table, skip this step and pass it ",
         "straight to beta_diversity().")
  }
  if (nrow(x) == 0) {
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

  site_ids <- trimws(as.character(x[[site]]))
  species_ids <- trimws(as.character(x[[species]]))

  no_site <- which(is.na(site_ids) | site_ids == "")
  if (length(no_site) > 0) {
    stop(length(no_site), " row(s) have a blank or missing site name in the '",
         site, "' column (row(s) ", .list_sites(no_site), "). Every record ",
         "needs to say which site it's from: fix or remove these rows.")
  }
  no_species <- which(is.na(species_ids) | species_ids == "")
  if (length(no_species) > 0) {
    stop(length(no_species), " row(s) have a blank or missing species name ",
         "in the '", species, "' column (row(s) ", .list_sites(no_species),
         "). The package can't tell which species these records are: fix ",
         "or remove these rows. (A species that was looked for but not ",
         "found doesn't need a row; use `value` or `all_sites` instead.)")
  }

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

  parsed <- .parse_all_sites(all_sites, site, group)
  mat <- .apply_all_sites(mat, parsed$sites)

  group_vec <- .resolve_group(x, group, site_ids, rownames(mat),
                              parsed$groups)

  list(x = mat, group = group_vec)
}

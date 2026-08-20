## ---- core reshape / value interpretation ----

test_that("logical value column reshapes into the correct wide matrix", {
  records <- data.frame(
    site = c("s1", "s1", "s2"),
    species = c("sp1", "sp2", "sp1"),
    seen = c(TRUE, FALSE, TRUE)
  )
  res <- as_site_by_species(records, site = "site", species = "species",
                             value = "seen", all_sites = c("s1", "s2"))

  expect_equal(dim(res$x), c(2, 2))
  expect_equal(res$x["s1", "sp1"], 1)
  expect_equal(res$x["s1", "sp2"], 0)
  expect_equal(res$x["s2", "sp1"], 1)
  expect_equal(res$x["s2", "sp2"], 0)
})

test_that("no value column supplied treats every row as a detection", {
  records <- data.frame(
    site = c("s1", "s1", "s2"),
    species = c("sp1", "sp2", "sp1")
  )
  res <- as_site_by_species(records, site = "site", species = "species",
                             all_sites = c("s1", "s2"))

  expect_equal(res$x["s1", "sp1"], 1)
  expect_equal(res$x["s1", "sp2"], 1)
  expect_equal(res$x["s2", "sp1"], 1)
  expect_equal(res$x["s2", "sp2"], 0)
})

test_that("numeric value column (counts/weights) thresholds at > 0", {
  records <- data.frame(
    site = c("s1", "s1", "s2"),
    species = c("sp1", "sp2", "sp1"),
    weight_g = c(12.4, 0, 3)
  )
  res <- as_site_by_species(records, site = "site", species = "species",
                             value = "weight_g", all_sites = c("s1", "s2"))

  expect_equal(res$x["s1", "sp1"], 1)
  expect_equal(res$x["s1", "sp2"], 0)
  expect_equal(res$x["s2", "sp1"], 1)
})

test_that("string value column accepts case-insensitive Y/N, yes/no, TRUE/FALSE", {
  records <- data.frame(
    site = rep("s1", 4),
    species = c("sp1", "sp2", "sp3", "sp4"),
    seen = c("Y", "no", "TRUE", "false")
  )
  res <- as_site_by_species(records, site = "site", species = "species",
                             value = "seen", all_sites = "s1")

  expect_equal(res$x["s1", "sp1"], 1)
  expect_equal(res$x["s1", "sp2"], 0)
  expect_equal(res$x["s1", "sp3"], 1)
  expect_equal(res$x["s1", "sp4"], 0)
})

test_that("unparseable string values in value error, listing the offending value", {
  records <- data.frame(
    site = c("s1", "s1"),
    species = c("sp1", "sp2"),
    seen = c("Y", "maybe")
  )
  expect_error(
    as_site_by_species(records, site = "site", species = "species", value = "seen"),
    "maybe"
  )
})

test_that("NA in value is treated as absence, with a warning", {
  records <- data.frame(
    site = c("s1", "s1"),
    species = c("sp1", "sp2"),
    seen = c(TRUE, NA)
  )
  expect_warning(
    res <- as_site_by_species(records, site = "site", species = "species",
                               value = "seen", all_sites = "s1"),
    "NA"
  )
  expect_equal(res$x["s1", "sp2"], 0)
})

## ---- duplicate (site, species) records ----

test_that("duplicate site/species rows collapse to present-if-detected-on-any-visit", {
  records <- data.frame(
    site = c("s1", "s1"),
    species = c("sp1", "sp1"),
    seen = c(FALSE, TRUE)
  )
  expect_message(
    res <- as_site_by_species(records, site = "site", species = "species",
                               value = "seen", all_sites = "s1"),
    "duplicate"
  )
  expect_equal(res$x["s1", "sp1"], 1)
})

test_that("duplicate collapsing emits a message, not a warning", {
  records <- data.frame(
    site = c("s1", "s1"),
    species = c("sp1", "sp1"),
    seen = c(TRUE, TRUE)
  )
  expect_no_warning(
    expect_message(
      as_site_by_species(records, site = "site", species = "species",
                          value = "seen", all_sites = "s1"),
      "duplicate"
    )
  )
})

test_that("no message is emitted when there are no duplicate records", {
  records <- data.frame(
    site = c("s1", "s2"),
    species = c("sp1", "sp1"),
    seen = c(TRUE, TRUE)
  )
  expect_no_message(
    as_site_by_species(records, site = "site", species = "species",
                        value = "seen", all_sites = c("s1", "s2"))
  )
})

## ---- implicit absence / all_sites ----

test_that("a zero-detection site is silently dropped when all_sites is omitted, with a warning", {
  records <- data.frame(site = "s1", species = "sp1")
  expect_warning(
    res <- as_site_by_species(records, site = "site", species = "species"),
    "all_sites"
  )
  expect_false("s2" %in% rownames(res$x))
})

test_that("all_sites adds correctly-shaped all-zero rows for undetected sites", {
  records <- data.frame(site = "s1", species = "sp1")
  res <- as_site_by_species(records, site = "site", species = "species",
                             all_sites = c("s1", "s2", "s3"))

  expect_equal(nrow(res$x), 3)
  expect_equal(unname(res$x["s2", ]), 0)
  expect_equal(unname(res$x["s3", ]), 0)
})

test_that("a site present in data but missing from all_sites errors, naming the site", {
  records <- data.frame(site = c("s1", "s2"), species = c("sp1", "sp1"))
  expect_error(
    as_site_by_species(records, site = "site", species = "species", all_sites = "s1"),
    "s2"
  )
})

test_that("duplicate IDs within all_sites are deduplicated, with a warning", {
  records <- data.frame(site = "s1", species = "sp1")
  expect_warning(
    res <- as_site_by_species(records, site = "site", species = "species",
                               all_sites = c("s1", "s1", "s2")),
    "duplicate"
  )
  expect_equal(nrow(res$x), 2)
})

## ---- group / category ----

test_that("a group column with consistent per-site labels returns an aligned group vector", {
  records <- data.frame(
    site = c("s1", "s1", "s2"),
    species = c("sp1", "sp2", "sp1"),
    habitat = c("tree", "tree", "open")
  )
  res <- as_site_by_species(records, site = "site", species = "species",
                             group = "habitat", all_sites = c("s1", "s2"))

  expect_equal(res$group, c("tree", "open")[match(rownames(res$x), c("s1", "s2"))])
})

test_that("a site with contradictory group labels errors, naming the site and the values", {
  records <- data.frame(
    site = c("s1", "s1"),
    species = c("sp1", "sp2"),
    habitat = c("tree", "open")
  )
  expect_error(
    as_site_by_species(records, site = "site", species = "species",
                        group = "habitat", all_sites = "s1"),
    "s1"
  )
})

test_that("group is NULL in the return value when not supplied", {
  records <- data.frame(site = "s1", species = "sp1")
  res <- as_site_by_species(records, site = "site", species = "species",
                             all_sites = "s1")
  expect_null(res$group)
})

test_that("group is NA, with a warning, for an all_sites-only zero-detection site", {
  records <- data.frame(site = "s1", species = "sp1", habitat = "tree")
  expect_warning(
    res <- as_site_by_species(records, site = "site", species = "species",
                               group = "habitat", all_sites = c("s1", "s2")),
    "group"
  )
  names(res$group) <- rownames(res$x)
  expect_true(is.na(res$group[["s2"]]))
})

test_that("the returned list feeds directly into beta_diversity(group = )", {
  set.seed(1)
  sites <- paste0("s", 1:20)
  records <- do.call(rbind, lapply(sites, function(s) {
    sp <- paste0("sp", sample(1:15, sample(3:8, 1)))
    data.frame(site = s, species = sp,
               habitat = if (as.integer(sub("s", "", s)) <= 10) "tree" else "open")
  }))
  res <- as_site_by_species(records, site = "site", species = "species",
                             group = "habitat", all_sites = sites)

  summary_tbl <- beta_diversity(res$x, group = res$group)
  expect_setequal(summary_tbl$group, c("tree", "open"))
})

## ---- validation / hygiene ----

test_that("a missing required column errors, naming the missing column and available columns", {
  records <- data.frame(location = "s1", taxon = "sp1")
  expect_error(
    as_site_by_species(records, site = "site", species = "species"),
    "site"
  )
})

test_that("an empty data.frame errors", {
  records <- data.frame(site = character(0), species = character(0))
  expect_error(
    as_site_by_species(records, site = "site", species = "species"),
    "empty"
  )
})

test_that("NA in the site column errors", {
  records <- data.frame(site = c("s1", NA), species = c("sp1", "sp2"))
  expect_error(
    as_site_by_species(records, site = "site", species = "species"),
    "site"
  )
})

test_that("NA in the species column errors", {
  records <- data.frame(site = c("s1", "s2"), species = c("sp1", NA))
  expect_error(
    as_site_by_species(records, site = "site", species = "species"),
    "species"
  )
})

test_that("leading/trailing whitespace in site or species IDs is trimmed automatically", {
  records <- data.frame(
    site = c(" s1", "s1 "),
    species = c("sp1", " sp1")
  )
  res <- as_site_by_species(records, site = "site", species = "species",
                             all_sites = "s1")

  expect_equal(nrow(res$x), 1)
  expect_equal(ncol(res$x), 1)
  expect_equal(rownames(res$x), "s1")
  expect_equal(colnames(res$x), "sp1")
})

test_that("site names differing only by case are kept distinct, with a warning", {
  records <- data.frame(
    site = c("SiteA", "sitea"),
    species = c("sp1", "sp2")
  )
  expect_warning(
    res <- as_site_by_species(records, site = "site", species = "species",
                               all_sites = c("SiteA", "sitea")),
    "case"
  )
  expect_equal(nrow(res$x), 2)
})

test_that("reshaping a long-format equivalent of the canonical 7x9 toy matrix matches it exactly", {
  wide <- matrix(
    c(1,1,0,1,0,0,1,0,0,
      1,1,1,0,0,1,0,0,0,
      0,1,1,1,1,0,0,0,0,
      1,0,1,1,0,0,1,1,0,
      0,1,0,0,1,1,0,0,0,
      1,1,1,1,1,0,0,0,0,
      0,0,1,0,1,1,1,0,1),
    nrow = 7, byrow = TRUE,
    dimnames = list(paste0("site", 1:7), paste0("sp", 1:9))
  )

  long <- do.call(rbind, lapply(seq_len(nrow(wide)), function(i) {
    present <- wide[i, ] > 0
    if (!any(present)) return(NULL)
    data.frame(site = rownames(wide)[i], species = colnames(wide)[present])
  }))

  res <- as_site_by_species(long, site = "site", species = "species",
                             all_sites = rownames(wide))

  expect_equal(res$x[rownames(wide), colnames(wide)], wide,
               ignore_attr = TRUE)

  expect_warning(
    summary_tbl <- beta_diversity(res$x),
    "below 10"
  )
  expect_equal(summary_tbl$gamma_obs, 9)
})

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
    "blank/missing"
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
    "Only N = "
  )
  expect_equal(summary_tbl$gamma_obs, 9)
})

## ---- CSV-shaped input (read.csv() quirks real users hit) ----

test_that("blank value cells in a CSV count as absence, with a warning (not a crash)", {
  records <- read.csv(text = paste(
    "site,species,seen",
    "s1,sp1,Y",
    "s1,sp2,",
    "s2,sp1,  ",
    "s2,sp3,yes",
    sep = "\n"))
  expect_type(records$seen, "character")  # read.csv reads blanks as ""
  expect_warning(
    res <- as_site_by_species(records, site = "site", species = "species",
                              value = "seen", all_sites = c("s1", "s2")),
    "2 row\\(s\\) have a blank/missing value"
  )
  expect_equal(res$x["s1", "sp2"], 0)
  expect_equal(res$x["s2", "sp1"], 0)
  expect_equal(res$x["s2", "sp3"], 1)
})

test_that("blank species or site names in a CSV error, giving the row numbers", {
  records <- read.csv(text = "site,species\ns1,sp1\ns1,\ns2,sp2\n")
  expect_error(
    as_site_by_species(records, site = "site", species = "species"),
    "row\\(s\\) 2\\)"
  )
  records <- read.csv(text = "site,species\ns1,sp1\n,sp1\n")
  expect_error(
    as_site_by_species(records, site = "site", species = "species"),
    "blank or missing site name"
  )
})

test_that("a matrix passed to as_site_by_species() gets a clear message, not '0 rows'", {
  m <- matrix(1, 2, 2)
  expect_error(as_site_by_species(m, site = "site", species = "species"),
               "needs a data frame")
})

test_that("all_sites as a table gives empty sites their group (hand-verified beta_M)", {
  # open: o1-o5 have 3 species each, o6 was surveyed but nothing was found.
  # tree: t1-t3, just to make a second group.
  records <- read.csv(text = paste(
    "site,species,habitat",
    "o1,a,open", "o1,b,open", "o1,c,open",
    "o2,a,open", "o2,b,open", "o2,d,open",
    "o3,a,open", "o3,c,open", "o3,e,open",
    "o4,b,open", "o4,c,open", "o4,f,open",
    "o5,a,open", "o5,b,open", "o5,c,open",
    "t1,a,tree", "t2,g,tree", "t3,h,tree",
    sep = "\n"))
  sites <- read.csv(text = paste(
    c("site,habitat,notes",
      paste0("o", 1:6, ",open,"), paste0("t", 1:3, ",tree,")),
    collapse = "\n"))

  res <- as_site_by_species(records, site = "site", species = "species",
                            group = "habitat", all_sites = sites)
  expect_equal(unname(res$group[rownames(res$x) == "o6"]), "open")
  expect_equal(sum(res$x["o6", ]), 0)

  out <- suppressWarnings(beta_diversity(res$x, group = res$group))
  open <- out[out$group == "open", ]
  # Hand calculation: N = 6, alpha = 15 / 6 = 2.5, gamma_obs = 6 (a-f),
  # u = 3 (d, e, f), gamma_hat = 6 + 3 * 5/6 = 8.5, beta_M = 8.5/2.5 - 1 = 2.4.
  # Without o6 it would be N = 5, alpha = 3, gamma_hat = 8.4, beta_M = 1.8:
  # leaving the empty site out makes beta_M too low.
  expect_equal(open$n_sites, 6)
  expect_equal(open$alpha, 2.5)
  expect_equal(open$gamma_hat, 8.5)
  expect_equal(open$beta_M, 2.4)
})

test_that("all_sites table errors: missing columns, conflicting groups", {
  records <- data.frame(site = c("s1", "s2"), species = "sp1",
                        habitat = c("tree", "open"))
  expect_error(
    as_site_by_species(records, site = "site", species = "species",
                       group = "habitat",
                       all_sites = data.frame(trap = c("s1", "s2"))),
    "needs a column named 'site'"
  )
  expect_error(
    as_site_by_species(records, site = "site", species = "species",
                       group = "habitat",
                       all_sites = data.frame(site = c("s1", "s2"))),
    "named 'habitat'"
  )
  expect_error(
    as_site_by_species(records, site = "site", species = "species",
                       group = "habitat",
                       all_sites = data.frame(site = c("s1", "s2"),
                                              habitat = c("open", "open"))),
    "s1 \\(records: tree vs all_sites: open\\)"
  )
  expect_error(
    as_site_by_species(records, site = "site", species = "species",
                       group = "habitat",
                       all_sites = data.frame(site = c("s1", "s1", "s2"),
                                              habitat = c("tree", "open", "open"))),
    "more than one 'habitat'"
  )
})

test_that("all_sites table without group just supplies the site list", {
  records <- data.frame(site = "s1", species = "sp1")
  res <- as_site_by_species(records, site = "site", species = "species",
                            all_sites = data.frame(site = c("s1", "s2")))
  expect_equal(rownames(res$x), c("s1", "s2"))
})

test_that("a blank habitat cell in a CSV becomes NA (not a group named \"\"), with a warning naming the site", {
  records <- read.csv(text = "site,species,habitat\ns1,sp1,tree\ns2,sp2,\ns3,sp1,open\n")
  expect_warning(
    res <- as_site_by_species(records, site = "site", species = "species",
                              group = "habitat", all_sites = c("s1", "s2", "s3")),
    "label: s2"
  )
  expect_equal(res$group, c("tree", NA, "open"))
})

test_that("a blank label on one row doesn't conflict with the site's label on another row", {
  records <- read.csv(text = "site,species,habitat\ns1,sp1,tree\ns1,sp2,\n")
  res <- as_site_by_species(records, site = "site", species = "species",
                            group = "habitat", all_sites = "s1")
  expect_equal(res$group, "tree")
})

test_that("the no-all_sites warning says missing empty sites make beta too low", {
  records <- data.frame(site = "s1", species = "sp1")
  expect_warning(as_site_by_species(records, site = "site", species = "species"),
                 "too low")
})

## ---- as_site_table() ----

# Shaped like a real field file: info columns first (a tree number that
# repeats across habitats, a text column, a richness total), then counts.
site_table_csv <- paste(
  "treenum,habitat,Specie,site_id,richness,aran_1,gna_2,liny_1,zod_1",
  "2,open,Ceratonia,2_open,2,0,1,0,3",
  "2,tree,Ceratonia,2_tree,3,1,2,1,0",
  "4,open,Ziziphus,4_open,0,0,0,0,0",
  "4,tree,Ziziphus,4_tree,1,0,0,0,1",
  sep = "\n")

test_that("species by name, position and \"first:last\" range give the same matrix", {
  tab <- read.csv(text = site_table_csv)
  by_range <- as_site_table(tab, site = "site_id", group = "habitat",
                            species = "aran_1:zod_1")
  by_pos <- as_site_table(tab, site = "site_id", group = "habitat",
                          species = 6:9)
  by_name <- as_site_table(tab, site = "site_id", group = "habitat",
                           species = c("aran_1", "gna_2", "liny_1", "zod_1"))
  expect_equal(by_pos, by_range)
  expect_equal(by_name, by_range)
  expect_equal(rownames(by_range$x), tab$site_id)
  expect_equal(colnames(by_range$x), c("aran_1", "gna_2", "liny_1", "zod_1"))
  expect_equal(by_range$group, tab$habitat)
  expect_equal(unname(by_range$x["2_open", ]), c(0, 1, 0, 3))
  # a site where nothing was found is kept, as a row of zeros
  expect_equal(unname(rowSums(by_range$x)["4_open"]), 0)
})

test_that("species = NULL uses every other column, and errors on text columns", {
  tab <- read.csv(text = site_table_csv)
  expect_error(as_site_table(tab, site = "site_id", group = "habitat"),
               "'Specie' contain text")
  num_only <- tab[c("site_id", "habitat", "aran_1", "gna_2", "liny_1", "zod_1")]
  res <- as_site_table(num_only, site = "site_id", group = "habitat")
  expect_equal(colnames(res$x), c("aran_1", "gna_2", "liny_1", "zod_1"))
})

test_that("info columns left in the species selection get a warning naming them", {
  tab <- read.csv(text = site_table_csv)
  expect_warning(
    as_site_table(tab, site = "site_id", group = "habitat",
                  species = c("treenum", "richness", "aran_1", "gna_2", "liny_1", "zod_1")),
    "'treenum'.*'richness'"
  )
})

test_that("a repeated site column (e.g. a tree number) errors, naming the value and rows", {
  tab <- read.csv(text = site_table_csv)
  expect_error(as_site_table(tab, site = "treenum", species = "aran_1:zod_1"),
               "'2' \\(rows 1, 2\\)")
  expect_error(as_site_table(tab, site = "treenum", species = "aran_1:zod_1"),
               "unique per row")
})

test_that("blank site names error with row numbers", {
  tab <- read.csv(text = site_table_csv)
  tab$site_id[3] <- ""
  expect_error(as_site_table(tab, site = "site_id", species = "aran_1:zod_1"),
               "row\\(s\\) 3\\)")
})

test_that("bad species selections give specific errors", {
  tab <- read.csv(text = site_table_csv)
  expect_error(as_site_table(tab, site = "site_id", species = "aran_1:nope"),
               "not found.*nope")
  expect_error(as_site_table(tab, site = "site_id", species = "zod_1:aran_1"),
               "comes after")
  expect_error(as_site_table(tab, site = "site_id", species = c("aran_1", "site_id")),
               "includes the site column")
  expect_error(as_site_table(tab, site = "site_id", species = 99),
               "between 1 and 9")
})

test_that("a blank group label becomes NA, with a warning naming the site", {
  tab <- read.csv(text = site_table_csv)
  tab$habitat[2] <- " "
  expect_warning(
    res <- as_site_table(tab, site = "site_id", group = "habitat",
                         species = "aran_1:zod_1"),
    "label: 2_tree"
  )
  expect_equal(res$group, c("open", NA, "open", "tree"))
})

test_that("non-data-frame input gets a clear message", {
  expect_error(as_site_table(matrix(1, 2, 2), site = "s"), "needs a data frame")
})

test_that("the hand-verified open example gives the same beta in both shapes", {
  # Same o1-o6 open sites as the all_sites test above, now one row per site.
  tab <- read.csv(text = paste(
    "site,habitat,a,b,c,d,e,f",
    "o1,open,1,1,1,0,0,0", "o2,open,1,1,0,1,0,0", "o3,open,1,0,1,0,1,0",
    "o4,open,0,1,1,0,0,1", "o5,open,1,1,1,0,0,0", "o6,open,0,0,0,0,0,0",
    sep = "\n"))
  res <- as_site_table(tab, site = "site", group = "habitat", species = "a:f")
  out <- suppressWarnings(beta_diversity(res$x, group = res$group))
  expect_equal(out$beta_M, 2.4)  # hand calculation in the all_sites test
  expect_equal(out$beta_W, 3.4)
})

test_that("the two example CSVs (record shape and site shape) give identical results", {
  records <- read.csv(system.file("extdata", "example_survey.csv",
                                  package = "kifbetadiversity"))
  sites <- read.csv(system.file("extdata", "example_sites.csv",
                                package = "kifbetadiversity"))
  tab <- read.csv(system.file("extdata", "example_site_table.csv",
                              package = "kifbetadiversity"))
  long <- as_site_by_species(records, site = "site", species = "species",
                             value = "count", group = "habitat",
                             all_sites = sites)
  wide <- as_site_table(tab, site = "site_id", group = "habitat",
                        species = "sp01:sp45")
  expect_equal(beta_diversity(wide$x, group = wide$group),
               beta_diversity(long$x, group = long$group))
})

test_that("as_site_by_species() on a one-row-per-site table points to as_site_table()", {
  tab <- read.csv(text = site_table_csv)
  expect_error(as_site_by_species(tab, site = "site_id", species = "species"),
               "use as_site_table")
})

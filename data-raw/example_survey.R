# Generates the synthetic example data shipped in inst/extdata/.
# Made-up data, for demonstrating the package only -- not a real survey.
#
# Two habitats, 30 sites each. "open" sites draw from a small pool of
# common species (sites look alike); "tree" sites draw from a larger pool
# of patchier species (sites differ more). One site per habitat was
# surveyed but nothing was found, so it appears only in example_sites.csv.
#
# Run from the package root: Rscript data-raw/example_survey.R

set.seed(2004)
n_sites <- 30
species <- sprintf("sp%02d", 1:45)

occupancy <- list(
  open = c(runif(18, 0.25, 0.6), rep(0.01, 27)),
  tree = c(runif(18, 0.15, 0.45), runif(27, 0.02, 0.10))
)

records <- list()
sites <- list()
for (habitat in names(occupancy)) {
  site_ids <- sprintf("%s_%02d", habitat, 1:n_sites)
  sites[[habitat]] <- data.frame(site = site_ids, habitat = habitat)
  for (s in site_ids[-n_sites]) {           # last site: nothing found
    found <- species[runif(length(species)) < occupancy[[habitat]]]
    if (length(found) == 0) next
    records[[s]] <- data.frame(site = s, species = found,
                               habitat = habitat,
                               count = rpois(length(found), 2) + 1)
  }
}

records <- do.call(rbind, records)
sites <- do.call(rbind, sites)
rownames(records) <- NULL
rownames(sites) <- NULL

write.csv(records, "inst/extdata/example_survey.csv", row.names = FALSE)
write.csv(sites, "inst/extdata/example_sites.csv", row.names = FALSE)

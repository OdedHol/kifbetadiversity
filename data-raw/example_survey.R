# Generates the synthetic example data shipped in inst/extdata/.
# Made-up data, for demonstrating the package only -- not a real survey.
#
# Two habitats, 30 sites each. "open" sites draw from a small pool of
# common species (sites look alike); "tree" sites draw from a larger pool
# of patchier species (sites differ more). One site per habitat was
# surveyed but nothing was found, so it appears only in example_sites.csv.
#
# The same survey is also written in the other common shape,
# example_site_table.csv: one row per site, one column per species (counts),
# plus info columns. plot_num repeats across habitats on purpose (like a
# tree number shared by a tree site and an open site), so it can't be used
# as the site name; richness is a summary column that must not be treated
# as a species.
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

# One row per site, species as columns
site_table <- sites
site_table$plot_num <- as.integer(sub(".*_", "", site_table$site))
counts <- matrix(0L, nrow = nrow(sites), ncol = length(species),
                 dimnames = list(sites$site, species))
counts[cbind(records$site, records$species)] <- as.integer(records$count)
site_table <- data.frame(plot_num = site_table$plot_num,
                         habitat = site_table$habitat,
                         site_id = site_table$site,
                         richness = rowSums(counts > 0),
                         counts, check.names = FALSE)
rownames(site_table) <- NULL
write.csv(site_table, "inst/extdata/example_site_table.csv", row.names = FALSE)

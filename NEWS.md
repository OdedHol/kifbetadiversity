# kifbetadiversity 0.2.0

* New `as_site_table()` prepares data with one row per site and one column
  per species, the other common field-data shape. Choose the species columns
  by name, position or a `"first:last"` range. It returns the same
  `list(x, group)` as `as_site_by_species()`.
* Results include `beta_W`, `beta_W_lwr` and `beta_W_upr`: conventional
  Whittaker beta, `gamma_hat / alpha` (= `beta_M + 1`), and its interval.
* A warning when columns treated as species look like IDs, dates, richness
  or abundance totals, or have out-of-scale values. These would otherwise
  silently change the results.
* Errors recognise the wrong data shape and say which function to use.
* README and vignette start with "Which shape is your data?", and a second
  example CSV, `example_site_table.csv`, holds the same survey in the
  one-row-per-site shape.

# kifbetadiversity 0.1.0

First shared version.

* Package renamed from betadiversity to kifbetadiversity, and installable
  with `remotes::install_github("OdedHol/kifbetadiversity")`.
* Only Whittaker's multiplicative beta (`beta_M`) is reported. The additive
  measure (`beta_A`) and its test were removed.
* `beta_diversity_test()` accepts one table plus `group = ` with two labels,
  and its output includes both group labels and both `beta_M` values.
* `as_site_by_species()`: `all_sites` can be a table of site + group, so
  sites where nothing was found keep their habitat.
* Blank value cells in a CSV count as absence instead of causing an error.
* Blank or missing group labels are reported as their own `NA` group, with a
  warning naming the sites, instead of being dropped (or, for blanks,
  becoming a group named `""`).
* When every site has the same species, the error says so, and a grouped
  run skips just that group.
* Clearer errors for a wide CSV read without `row.names = 1`, for a
  non-data-frame passed to `as_site_by_species()`, and for records with no
  species or site name (giving the row numbers).
* The warning without `all_sites` now says that missing empty sites make
  beta diversity too low.
* Example data in `inst/extdata/`, a rewritten vignette with an
  "Assumptions and limitations" section, and the same section in
  `?beta_diversity_test`.

# betadiversity 0.0.0.9000

* Initial development version.

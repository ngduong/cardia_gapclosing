# Checks the installation by recomputing the point estimates for imputation 1
# (every intervention, outcome definition, scale and outcome-model curve) and
# comparing them with tests/expected_imputation1.csv.  Takes a few minutes on
# one core.  Exits with an error if any estimate differs by more than 1e-6.

for (f in list.files("R", full.names = TRUE)) source(f)
plan(sequential)
cat(sprintf("%s | mgcv %s | eventglm %s | survival %s\n", R.version.string,
            packageVersion("mgcv"), packageVersion("eventglm"), packageVersion("survival")))

expected <- read_csv("tests/expected_imputation1.csv", show_col_types = FALSE)
dat1 <- attach_imputation(read_base_data(), read_le8_imputations(m_max = 1), 1)
worst <- 0
for (i in seq_len(nrow(distinct(expected, outcome, scale, le8_curve, intervention)))) {
  key <- distinct(expected, outcome, scale, le8_curve, intervention)[i, ]
  res <- estimate_intervention(dat1, INTERVENTIONS[[key$intervention]], key$outcome, key$scale,
                               le8_by_race = key$le8_curve == "by_race")
  got <- res$disparity %>% select(weighting, gap_closed)
  ref <- expected %>% filter(outcome == key$outcome, scale == key$scale, le8_curve == key$le8_curve,
                             intervention == key$intervention) %>% select(weighting, expected = gap_closed)
  d <- inner_join(got, ref, by = "weighting") %>% mutate(diff = abs(gap_closed - expected))
  worst <- max(worst, d$diff)
  cat(sprintf("%-9s %-4s %-7s %-40s max |diff| = %.1e\n", key$outcome, key$scale, key$le8_curve, key$intervention, max(d$diff)))
}
cat(sprintf("\nlargest difference from the reference values: %.2e -> %s\n", worst,
            if (worst < 1e-6) "installation reproduces the reference values" else "MISMATCH"))
if (worst >= 1e-6) quit(status = 2)

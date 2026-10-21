# Pools the saved analyses into result tables.
#   Rscript summarize_results.R            # every folder under results/
#   Rscript summarize_results.R <tag> ...  # selected folders
#
# Estimates are averaged over imputations, and their variances combine the
# within-imputation bootstrap variance with the between-imputation variance
# (Rubin's rules).  Intervals are Wald 95% intervals.  Written to
#   results/summary_<tag>.csv          disparities, gap closed, proportion closed
#   results/summary_<tag>_by_race.csv  race-specific outcomes, observed and under each intervention

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
source("R/specifications.R")

args <- commandArgs(trailingOnly = TRUE)
tags <- if (length(args)) args else list.dirs("results", recursive = FALSE, full.names = FALSE)

to_wide <- function(d, ids) {
  d %>%
    pivot_wider(id_cols = all_of(ids), names_from = quantity, values_from = c(estimate, se, lower, upper),
                names_glue = "{quantity}_{.value}") %>%
    arrange(match(intervention, names(INTERVENTIONS)))
}

for (tag in tags) {
  files <- list.files(file.path("results", tag), pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) next
  fits <- lapply(files, readRDS)
  disparity <- to_wide(bind_rows(lapply(fits, `[[`, "pooled_disparity")),
                       c("intervention", "outcome", "scale", "weighting"))
  by_race <- to_wide(bind_rows(lapply(fits, `[[`, "pooled_means")),
                     c("intervention", "outcome", "scale", "weighting", "race_label"))
  write_csv(disparity, file.path("results", paste0("summary_", tag, ".csv")))
  write_csv(by_race, file.path("results", paste0("summary_", tag, "_by_race.csv")))

  # risk differences in percentage points; restricted mean times in years
  k <- if (disparity$scale[1] == "risk") 100 else 1
  cat(sprintf("\n%s (weighted; gap closed in %s, proportion in %%)\n", tag, if (k == 100) "pp" else "years"))
  print(as.data.frame(disparity %>% filter(weighting == "weighted") %>%
    transmute(intervention,
              gap_closed = sprintf("%.2f (%.2f, %.2f)", k * gap_closed_estimate, k * gap_closed_lower, k * gap_closed_upper),
              proportion_closed = sprintf("%.1f (%.1f, %.1f)", 100 * proportion_closed_estimate,
                                          100 * proportion_closed_lower, 100 * proportion_closed_upper))),
    row.names = FALSE)
}

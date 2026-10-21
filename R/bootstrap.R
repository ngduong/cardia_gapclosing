# Bootstrap within imputation and Rubin's rules across imputations.
#
# Participants are resampled within race.  Every replicate repeats the whole
# estimation (pseudo-observations, weights, LE8 models, intervention).
# Replicates use independent L'Ecuyer-CMRG streams seeded from the imputation
# number, so results do not depend on the number of workers.

bootstrap_intervention <- function(dat, spec, outcome, scale, le8_by_race = FALSE,
                                   B = 500, seed = 2026, race = RACE) {
  point <- estimate_intervention(dat, spec, outcome, scale, le8_by_race, with_distribution = TRUE)

  one_replicate <- function(b) {
    resample <- dat %>% group_by(across(all_of(race))) %>%
      slice_sample(prop = 1, replace = TRUE) %>% ungroup()
    out <- tryCatch(estimate_intervention(resample, spec, outcome, scale, le8_by_race),
                    error = function(e) NULL)
    rm(resample); gc(verbose = FALSE)
    if (is.null(out)) NULL else list(means = mutate(out$means, b = b), disparity = mutate(out$disparity, b = b))
  }

  workers <- max(1L, nbrOfWorkers())
  seeds <- future:::make_rng_seeds(B, seed = seed)
  run <- function(idx) future_lapply(idx, one_replicate, future.seed = seeds[idx],
                                     future.chunk.size = ceiling(length(idx) / workers),
                                     future.packages = c("dplyr", "tidyr", "tibble", "survival",
                                                         "eventglm", "mgcv", "splines"))
  replicates <- NULL
  for (attempt in 1:3) {                       # a crashed worker aborts the call; rerun with the same seeds
    replicates <- tryCatch(run(seq_len(B)), error = function(e) NULL)
    if (!is.null(replicates)) break
  }
  if (is.null(replicates)) {                   # isolate the failure to one batch
    replicates <- vector("list", B)
    for (batch in split(seq_len(B), ceiling(seq_len(B) / workers))) {
      r <- tryCatch(run(batch), error = function(e) NULL)
      if (!is.null(r)) replicates[batch] <- r
    }
  }
  replicates <- Filter(Negate(is.null), replicates)
  if (length(replicates) < 0.9 * B)
    warning(sprintf("only %d of %d bootstrap replicates succeeded", length(replicates), B))

  rep_means <- bind_rows(lapply(replicates, `[[`, "means"))
  rep_disp  <- bind_rows(lapply(replicates, `[[`, "disparity"))
  se_means <- rep_means %>% group_by(intervention, outcome, scale, weighting, across(all_of(race))) %>%
    summarize(se_estimate = sd(estimate), se_observed = sd(observed), .groups = "drop")
  se_disp <- rep_disp %>% group_by(intervention, outcome, scale, weighting) %>%
    summarize(se_factual_disparity = sd(factual_disparity),
              se_counterfactual_disparity = sd(counterfactual_disparity),
              se_gap_closed = sd(gap_closed), se_proportion_closed = sd(proportion_closed), .groups = "drop")
  list(means = left_join(point$means, se_means, by = c("intervention", "outcome", "scale", "weighting", race)),
       disparity = left_join(point$disparity, se_disp, by = c("intervention", "outcome", "scale", "weighting")),
       intervention = point$intervention, distribution = point$distribution,
       replicates = list(means = rep_means, disparity = rep_disp, requested = B, succeeded = length(replicates)))
}

# Run every imputation and pool with Rubin's rules.
run_imputations <- function(base, le8_long, spec, outcome, scale, le8_by_race = FALSE,
                            M = 50, B = 500, replicate_dir = NULL, log_file = NULL) {
  per_imp <- vector("list", M)
  for (m in seq_len(M)) {
    log_line(log_file, sprintf("  %s: imputation %d of %d", spec$name, m, M))
    res <- bootstrap_intervention(attach_imputation(base, le8_long, m), spec, outcome, scale,
                                  le8_by_race, B = B, seed = 2026 + 100 * m)
    if (!is.null(replicate_dir))
      saveRDS(res$replicates, file.path(replicate_dir, sprintf("%s_imp%02d.rds", spec$name, m)))
    per_imp[[m]] <- list(
      means = mutate(res$means, imp = m),
      disparity = mutate(res$disparity, imp = m),
      params = as_tibble(c(list(imp = m, label = res$intervention$label), res$intervention$params)),
      distribution = mutate(res$distribution$summary, imp = m),
      bins = if (m == 1) res$distribution$bins else NULL,
      replicates_succeeded = res$replicates$succeeded)
    gc()
  }
  means <- bind_rows(lapply(per_imp, `[[`, "means"))
  disparity <- bind_rows(lapply(per_imp, `[[`, "disparity"))
  list(means = means, disparity = disparity,
       params = bind_rows(lapply(per_imp, `[[`, "params")),
       distribution = bind_rows(lapply(per_imp, `[[`, "distribution")),
       bins = per_imp[[1]]$bins,
       replicates_succeeded = sapply(per_imp, `[[`, "replicates_succeeded"),
       pooled_means = pool_rubin(means, c("intervention", "outcome", "scale", "weighting", RACE),
                                 c(estimate = "se_estimate", observed = "se_observed")),
       pooled_disparity = pool_rubin(disparity, c("intervention", "outcome", "scale", "weighting"),
                                     c(factual_disparity = "se_factual_disparity",
                                       counterfactual_disparity = "se_counterfactual_disparity",
                                       gap_closed = "se_gap_closed", proportion_closed = "se_proportion_closed")))
}

pool_rubin <- function(d, groups, estimates) {
  bind_rows(lapply(names(estimates), function(est) {
    se <- estimates[[est]]
    d %>% group_by(across(all_of(groups))) %>%
      summarize(quantity = est, M = n_distinct(imp), estimate = mean(.data[[est]]),
                within = mean(.data[[se]]^2), between = coalesce(var(.data[[est]]), 0), .groups = "drop") %>%
      mutate(se = sqrt(within + (1 + 1 / M) * between),
             lower = estimate - qnorm(0.975) * se, upper = estimate + qnorm(0.975) * se)
  }))
}

# Doubly robust estimation of race-specific outcomes under an intervention.
#
# For each race, the estimate is
#     plug-in + sum(w r (Y - Q)) / sum(w r)
# where Q is the outcome regression evaluated at the intervention score, r the
# density ratio of the post-intervention to the observed LE8 distribution, and
# w the allowability weight.  The ratio is taken from a censored-Gaussian
# (Tobit) model for LE8, which places a point mass at the 100-point ceiling, so
# the capped intervention has an exact ratio:
#     a < 100:  f(a - s) 1{a - s < tau} / f(a) + 1{a >= tau}
#     a = 100:  P(d(LE8) = 100) / P(LE8 = 100)

# Tobit model for LE8: LE8 = min(LE8*, 100), LE8* ~ N(mu, sigma^2).
fit_le8_tobit <- function(dat, race = RACE, ceiling = SCORE_BOUNDS[2]) {
  d <- dat
  d$.y <- d[[TREATMENT]]; d$.uncensored <- as.numeric(d$.y < ceiling - 1e-9)
  f <- as.formula(paste0("survival::Surv(.y, .uncensored) ~ ", race, " * (sex + a03ed + c03incom) + ",
                         race, ":splines::ns(examage, 3) + ", race, ":splines::ns(deprv80n00, 3)"))
  fit <- survreg(f, data = d, dist = "gaussian")
  list(mu = as.numeric(predict(fit, newdata = d, type = "lp")), sigma = fit$scale)
}

# Density ratio r = g / f at the observed scores.
density_ratio <- function(dat, intervention, tobit, ceiling = SCORE_BOUNDS[2]) {
  a <- dat[[TREATMENT]]; s <- intervention$shift; tau <- intervention$tau
  mu <- tobit$mu; sg <- tobit$sigma
  f <- dnorm(a, mu, sg)
  g <- dnorm(a - s, mu, sg) * ((a - s < tau) & (a - s < ceiling)) + f * (a >= tau)
  at_ceiling <- a >= ceiling - 1e-9
  p_ceiling <- pnorm((ceiling - mu) / sg, lower.tail = FALSE)
  p_ceiling_post <- if (is.infinite(tau)) {
    ifelse(s >= 0, pnorm((ceiling - s - mu) / sg, lower.tail = FALSE), 0)
  } else {
    p_ceiling + ifelse(s > 0, pmax(0, pnorm((tau - mu) / sg) - pnorm((ceiling - s - mu) / sg)), 0)
  }
  f[at_ceiling] <- p_ceiling[at_ceiling]; g[at_ceiling] <- p_ceiling_post[at_ceiling]
  g / f
}

truncate_within_race <- function(r, race, probs = c(0.001, 0.999)) {
  for (g in unique(race)) {
    k <- race == g; q <- quantile(r[k], probs)
    r[k] <- pmin(pmax(r[k], q[1]), q[2])
  }
  r
}

outcome_formula <- function(le8_by_race = FALSE, race = RACE) {
  as.formula(paste0("pseudo_outcome ~ s(", TREATMENT, if (le8_by_race) ", by = .race_f" else "", ") + ",
                    race, " + s(examage) + sex + a03ed + c03incom + s(deprv80n00)"))
}

# Race-specific outcomes under the intervention for one weighting scheme.
estimate_outcomes <- function(dat, intervention, tobit, weighted = TRUE,
                              le8_by_race = FALSE, truncation = c(0.001, 0.999), race = RACE) {
  w <- if (weighted) dat$allow_weight else rep(1, nrow(dat))
  w <- w / mean(w)
  dat$.w <- w; dat$.race_f <- factor(dat[[race]])
  fit <- gam(outcome_formula(le8_by_race), data = dat, weights = .w)
  y <- dat$pseudo_outcome
  resid <- y - as.numeric(predict(fit, newdata = dat))
  post <- dat; post[[TREATMENT]] <- intervention$score
  q_post <- as.numeric(predict(fit, newdata = post))
  r <- truncate_within_race(density_ratio(dat, intervention, tobit), dat[[race]], truncation)

  per_race <- function(g) {
    k <- dat[[race]] == g
    plug <- weighted.mean(q_post[k], w[k])
    tibble(!!race := g,
           observed = weighted.mean(y[k], w[k]),
           plug_in = plug,
           ipw = weighted.mean(y[k], w[k] * r[k]),
           estimate = plug + weighted.mean(resid[k], w[k] * r[k]))
  }
  bind_rows(per_race("Black"), per_race("White")) %>%
    mutate(weighting = if (weighted) "weighted" else "unweighted")
}

# Full estimate for one dataset: pseudo-outcomes, weights, LE8 models, the
# intervention and both weighting schemes.
estimate_intervention <- function(dat, spec, outcome, scale, le8_by_race = FALSE,
                                  with_distribution = FALSE, race = RACE) {
  dat <- add_pseudo_outcome(dat, outcome, scale) %>%
    drop_na(all_of(c(race, "pseudo_outcome", TREATMENT, ALLOWABLE, NON_ALLOWABLE))) %>%
    allowability_weights()
  le8_fit <- fit_le8_model(dat)
  intervention <- build_intervention(dat, le8_fit, spec)
  tobit <- fit_le8_tobit(dat)
  means <- bind_rows(estimate_outcomes(dat, intervention, tobit, TRUE,  le8_by_race),
                     estimate_outcomes(dat, intervention, tobit, FALSE, le8_by_race)) %>%
    mutate(intervention = spec$name, outcome = outcome, scale = scale, .before = 1)
  disparity <- means %>%
    group_by(intervention, outcome, scale, weighting) %>%
    summarize(across(c(observed, plug_in, ipw, estimate),
                     ~ .x[.data[[race]] == "Black"] - .x[.data[[race]] == "White"]), .groups = "drop") %>%
    rename(factual_disparity = observed, counterfactual_disparity = estimate) %>%
    mutate(gap_closed = factual_disparity - counterfactual_disparity,
           proportion_closed = gap_closed / factual_disparity,
           gap_closed_plug_in = factual_disparity - plug_in,
           gap_closed_ipw = factual_disparity - ipw) %>%
    select(-plug_in, -ipw)
  list(means = means, disparity = disparity, intervention = intervention,
       distribution = if (with_distribution) summarize_le8_distribution(dat, intervention) else NULL)
}

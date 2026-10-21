# Hypothetical interventions on the LE8 score.
#
# Every intervention is a capped additive shift,
#     d(a) = min(max(a + s, 0), 100)   applied when a < tau, otherwise d(a) = a,
# where the shift s may depend on race and covariates and tau = Inf for
# universal rules.  Each constructor returns the per-participant shift, the
# threshold and the parameters that define the rule.

# Conditional mean of LE8 given race and covariates, used to define the
# equalization rules.  Race interacts with every covariate.
fit_le8_model <- function(dat, race = RACE, non_allowable = NON_ALLOWABLE) {
  f <- as.formula(paste0(
    TREATMENT, " ~ ", race, " + s(examage, by = factor(", race, ")) + ", race, " * sex + ",
    paste(sapply(non_allowable, function(v)
      if (v == "deprv80n00") paste0("s(", v, ", by = factor(", race, "))") else paste0(race, " * ", v)),
      collapse = " + ")))
  gam(f, data = dat, family = gaussian())
}

le8_model_mean <- function(fit, dat) as.numeric(predict(fit, newdata = dat))

# Expected LE8 under the White covariate process, for every participant.
white_predicted_le8 <- function(dat, race = RACE) {
  f <- as.formula(paste(TREATMENT, "~", paste(c(ALLOWABLE, NON_ALLOWABLE), collapse = " + ")))
  fit <- lm(f, data = dat[dat[[race]] == "White", ])
  as.numeric(predict(fit, newdata = dat))
}

counterfactual_score <- function(a, shift, tau = Inf, bounds = SCORE_BOUNDS) {
  ifelse(a < tau, pmin(pmax(a + shift, bounds[1]), bounds[2]), a)
}

# Scalar added to the Black shift so that the capped Black mean equals a target.
solve_parity <- function(a, shift, tau, black, target) {
  optimize(function(k) {
    s <- shift; s[black] <- s[black] + k
    (mean(counterfactual_score(a, s, tau)[black]) - target)^2
  }, c(-50, 50))$minimum
}

# 1. Universal increase: every participant gains delta.
intervention_universal <- function(dat, le8_fit, delta = 10, ...) {
  list(type = "universal", shift = rep(delta, nrow(dat)), tau = Inf,
       label = sprintf("universal_%g", delta), params = list(delta = delta))
}

# 2. CVH-sensitive increase: participants below tau gain delta.
intervention_threshold <- function(dat, le8_fit, delta = 10, tau = 68, ...) {
  list(type = "threshold", shift = rep(delta, nrow(dat)), tau = tau,
       label = sprintf("threshold_%g_%g", delta, tau), params = list(delta = delta, tau = tau))
}

# 3. Marginal equalization: Black participants gain the Black-White difference
#    in mean LE8.
intervention_marginal_equalization <- function(dat, le8_fit, race = RACE, ...) {
  a <- dat[[TREATMENT]]; black <- dat[[race]] == "Black"
  delta <- mean(a[!black]) - mean(a[black])
  list(type = "marginal_equalization", shift = ifelse(black, delta, 0), tau = Inf,
       label = "marginal_equalization", params = list(delta = delta))
}

# 4. Conditional equalization: each Black participant gains the difference
#    between the White-predicted and own-predicted LE8 at their covariates,
#    with negative differences set to zero.
intervention_conditional_equalization <- function(dat, le8_fit, race = RACE, ...) {
  black <- dat[[race]] == "Black"
  gap <- white_predicted_le8(dat) - le8_model_mean(le8_fit, dat)
  list(type = "conditional_equalization", shift = ifelse(black, pmax(gap, 0), 0), tau = Inf,
       label = "conditional_equalization",
       params = list(mean_shift = mean(pmax(gap, 0)[black]), share_floored = mean(gap[black] < 0)))
}

# 5. Conditional equalization plus a constant supplement chosen so that the
#    capped Black mean equals the White mean.
intervention_conditional_supplement <- function(dat, le8_fit, race = RACE, ...) {
  a <- dat[[TREATMENT]]; black <- dat[[race]] == "Black"
  gap <- pmax(white_predicted_le8(dat) - le8_model_mean(le8_fit, dat), 0)
  base <- ifelse(black, gap, 0)
  supplement <- solve_parity(a, base, Inf, black, mean(a[!black]))
  list(type = "conditional_supplement", shift = ifelse(black, gap + supplement, 0), tau = Inf,
       label = "conditional_equalization_supplement",
       params = list(supplement = supplement, share_floored = mean((white_predicted_le8(dat) - le8_model_mean(le8_fit, dat))[black] < 0)))
}

# 6-9. Combined rules: a race-blind increase of delta for everyone eligible,
#      plus for Black participants either a constant (marginal parity) or the
#      conditional gap (not floored) and a constant, the constant solved so
#      that the capped Black mean equals the capped White mean.
intervention_combined <- function(dat, le8_fit, delta = 10, tau = Inf,
                                  conditional = FALSE, race = RACE, ...) {
  a <- dat[[TREATMENT]]; black <- dat[[race]] == "Black"
  shift <- rep(delta, nrow(dat))
  if (conditional)
    shift[black] <- (white_predicted_le8(dat) - le8_model_mean(le8_fit, dat))[black] + delta
  target <- mean(counterfactual_score(a, shift, tau)[!black])
  kappa <- solve_parity(a, shift, tau, black, target)
  shift[black] <- shift[black] + kappa
  type <- if (is.infinite(tau)) "combined_universal" else "combined_threshold"
  list(type = type, shift = shift, tau = tau,
       label = sprintf("%s_%s_%g%s", type, if (conditional) "conditional" else "marginal", delta,
                       if (is.finite(tau)) sprintf("_%g", tau) else ""),
       params = list(delta = delta, tau = tau, conditional = conditional, kappa = kappa,
                     white_mean_post = target,
                     black_mean_post = mean(counterfactual_score(a, shift, tau)[black])))
}

build_intervention <- function(dat, le8_fit, spec) {
  constructor <- switch(spec$type,
    universal                = intervention_universal,
    threshold                = intervention_threshold,
    marginal_equalization    = intervention_marginal_equalization,
    conditional_equalization = intervention_conditional_equalization,
    conditional_supplement   = intervention_conditional_supplement,
    combined                 = intervention_combined,
    stop("unknown intervention type: ", spec$type))
  obj <- do.call(constructor, c(list(dat = dat, le8_fit = le8_fit), spec$params))
  obj$score <- counterfactual_score(dat[[TREATMENT]], obj$shift, obj$tau)
  obj
}

# Observed and post-intervention LE8 distributions by race (allowability weighted).
summarize_le8_distribution <- function(dat, intervention, race = RACE, breaks = seq(0, 100, 10)) {
  w <- dat$allow_weight
  one <- function(x, timing) {
    tibble(!!race := dat[[race]], a = x, w = w) %>%
      group_by(across(all_of(race))) %>%
      summarize(timing = timing,
                mean = weighted.mean(a, w), sd = sqrt(sum(w * (a - weighted.mean(a, w))^2) / sum(w)),
                p10 = Hmisc::wtd.quantile(a, w, 0.10), p25 = Hmisc::wtd.quantile(a, w, 0.25),
                p50 = Hmisc::wtd.quantile(a, w, 0.50), p75 = Hmisc::wtd.quantile(a, w, 0.75),
                p90 = Hmisc::wtd.quantile(a, w, 0.90),
                share_high = weighted.mean(a >= 80, w), n = n(), .groups = "drop")
  }
  bins <- function(x, timing) {
    tibble(!!race := dat[[race]], bin = cut(x, breaks, include.lowest = TRUE, right = FALSE), w = w) %>%
      group_by(across(all_of(c(race, "bin")))) %>% summarize(timing = timing, weight = sum(w), .groups = "drop") %>%
      group_by(across(all_of(race))) %>% mutate(share = weight / sum(weight)) %>% ungroup()
  }
  list(summary = bind_rows(one(dat[[TREATMENT]], "observed"), one(intervention$score, "intervention")),
       bins = bind_rows(bins(dat[[TREATMENT]], "observed"), bins(intervention$score, "intervention")))
}

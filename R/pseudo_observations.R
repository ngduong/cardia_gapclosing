# Pseudo-observations for the 35-year outcome.
#
# Within strata of race and baseline CVH (LE8 >= 80), jackknife pseudo-values
# are computed with inverse-probability-of-censoring weights from a Cox model
# for dropout (eventglm, Hajek weights).  Two outcome definitions:
#   "censored"   non-CVD death treated as censoring
#   "competing"  non-CVD death treated as a competing event (no MACE)
# and two scales:
#   "risk"  cumulative incidence of MACE by 35 years
#   "rmst"  restricted mean survival time (censored) or restricted mean time
#           lost to MACE (competing) by 35 years

OUTCOME_SPEC <- list(
  censored  = list(time = "t_cvd_or_cvd_death", status = "cvd_or_cvd_death_event", competing = FALSE),
  competing = list(time = "mace_time",          status = "mace_status",            competing = TRUE))

add_pseudo_outcome <- function(dat, outcome = c("censored", "competing"),
                               scale = c("risk", "rmst"), horizon = HORIZON,
                               censoring_formula = CENSORING_FORMULA,
                               strata = c(RACE, "le8_high")) {
  outcome <- match.arg(outcome); scale <- match.arg(scale)
  spec <- OUTCOME_SPEC[[outcome]]
  if (!"le8_high" %in% names(dat))
    dat$le8_high <- as.integer(dat[[TREATMENT]] >= 80)
  dat$.row <- seq_len(nrow(dat))
  needed <- unique(c(".row", strata, spec$time, spec$status, all.vars(censoring_formula)))
  cc <- drop_na(dat, all_of(needed))
  if (spec$competing) cc[[spec$status]] <- as.factor(cc[[spec$status]])
  else                cc[[spec$status]] <- as.numeric(cc[[spec$status]])
  f <- as.formula(sprintf("survival::Surv(%s, %s) ~ 1", spec$time, spec$status))

  one_stratum <- function(d, g) {
    args <- list(formula = f, data = d, time = horizon, model.censoring = "coxph",
                 formula.censoring = censoring_formula, ipcw.method = "hajek")
    if (scale == "risk") args$survival <- FALSE
    if (spec$competing) args$cause <- "MACE"
    fit <- tryCatch(do.call(if (scale == "risk") eventglm::cumincglm else eventglm::rmeanglm, args),
                    error = function(e) NULL)
    y <- if (is.null(fit) || length(fit$y) != nrow(d)) NA_real_ else as.numeric(fit$y)
    tibble(.row = d$.row, pseudo = y)
  }
  po <- cc %>% group_by(across(all_of(strata))) %>% group_modify(one_stratum) %>% ungroup()
  dat$pseudo_outcome <- NA_real_
  dat$pseudo_outcome[po$.row] <- po$pseudo
  dat$.row <- NULL
  dat
}

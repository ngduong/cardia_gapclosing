# Hypothetical interventions on cardiovascular health and the Black-White disparity in cardiovascular events

Code for the analysis of the CARDIA study described in the accompanying manuscript. Nine hypothetical interventions on early-adulthood Life's Essential 8 (LE8) scores are evaluated as modified treatment policies, and the share of the Black-White disparity in 35-year risk of major adverse cardiovascular events that each would close is estimated with a doubly robust gap-closing estimator.

## Layout

| File | Contents |
|---|---|
| `R/setup.R` | Packages, constants, data access |
| `R/pseudo_observations.R` | Jackknife pseudo-observations for 35-year risk and restricted mean time, with non-CVD death either censored or treated as a competing event |
| `R/allowability_weights.R` | Odds weights standardizing White participants to the Black age-sex distribution |
| `R/interventions.R` | The nine interventions as capped additive shifts of the LE8 score |
| `R/estimator.R` | Censored-Gaussian density ratio and the doubly robust estimator |
| `R/bootstrap.R` | Within-race bootstrap and Rubin's rules across imputations |
| `R/specifications.R` | Primary and sensitivity intervention parameters |
| `run_analysis.R` | One analysis: outcome definition x intervention set x scale |
| `run_all.R` | Every analysis in sequence |
| `summarize_results.R` | Pooled tables with Wald intervals |
| `check_installation.R` | Recomputes imputation 1 and compares with `tests/expected_imputation1.csv` |

## Data

The analysis uses CARDIA data, which are not publicly available. For ease of conceptualization, we describe the analytical dataset layout and variables in the analytical dataset:

- `cdm_data_clean_ipcw_imp_v2.csv`: one row per participant with baseline covariates (`examage`, `sex`, `race_label`, `a03ed`, `c03incom`, `deprv80n00`, `a10smoke`, `a07mlalc`), follow-up for MACE with non-CVD death censored (`t_cvd_or_cvd_death`, `cvd_or_cvd_death_event`) and with competing events (`mace_time`, `mace_status`).
- `le8_score_all.csv`: multiply imputed LE8 scores, long by imputation (`id`, `exam`, `imp`, `le8_imputed`).

## Method summary

- Outcome: 35-year cumulative incidence of MACE (or restricted mean time), as pseudo-observations within strata of race and baseline CVH estimated using a Cox model for dropout.
- Outcome regression: Gaussian GAM of the pseudo-outcome on the LE8 score, race, age, sex, education, income and neighborhood deprivation.
- LE8 model: GAM with race interacting with every covariate; a censored-Gaussian (Tobit) model is used to estimate the density ratio, which is truncated at the 100-point ceiling. IPW weights are truncated at the 0.1st and 99.9th percentiles within race.
- Inference: nonparametric bootstrap within race (500 replicates) in each of 50 imputations, pooled with Rubin's rules. We construct 95% Wald-type confidence intervals.

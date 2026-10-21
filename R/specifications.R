# Intervention specifications.
#
# Primary analysis: nine interventions with a 10-point increase and a
# CVH-sensitive threshold of 68 (the lowest quartile).  Sensitivity analysis:
# the same rules with a 5-point increase, and with a threshold of 80 (the
# high-CVH boundary).

make_spec <- function(name, number, type, ...) list(name = name, number = number, type = type, params = list(...))

INTERVENTIONS <- list(
  make_spec("universal_10",                           1, "universal", delta = 10),
  make_spec("cvh_sensitive_10_68",                    2, "threshold", delta = 10, tau = 68),
  make_spec("marginal_equalization",                  3, "marginal_equalization"),
  make_spec("conditional_equalization",               4, "conditional_equalization"),
  make_spec("conditional_equalization_supplement",    5, "conditional_supplement"),
  make_spec("universal_marginal_equalization",        6, "combined", delta = 10, tau = Inf, conditional = FALSE),
  make_spec("universal_conditional_supplement",       7, "combined", delta = 10, tau = Inf, conditional = TRUE),
  make_spec("cvh_sensitive_marginal_equalization",    8, "combined", delta = 10, tau = 68,  conditional = FALSE),
  make_spec("cvh_sensitive_conditional_supplement",   9, "combined", delta = 10, tau = 68,  conditional = TRUE),
  # sensitivity: 5-point increase
  make_spec("universal_5",                            1, "universal", delta = 5),
  make_spec("cvh_sensitive_5_68",                     2, "threshold", delta = 5, tau = 68),
  make_spec("universal_marginal_equalization_5",      6, "combined", delta = 5, tau = Inf, conditional = FALSE),
  make_spec("universal_conditional_supplement_5",     7, "combined", delta = 5, tau = Inf, conditional = TRUE),
  make_spec("cvh_sensitive_marginal_equalization_5",  8, "combined", delta = 5, tau = 68,  conditional = FALSE),
  make_spec("cvh_sensitive_conditional_supplement_5", 9, "combined", delta = 5, tau = 68,  conditional = TRUE),
  # sensitivity: threshold of 80
  make_spec("cvh_sensitive_10_80",                    2, "threshold", delta = 10, tau = 80),
  make_spec("cvh_sensitive_marginal_equalization_80", 8, "combined", delta = 10, tau = 80, conditional = FALSE),
  make_spec("cvh_sensitive_conditional_supplement_80", 9, "combined", delta = 10, tau = 80, conditional = TRUE))
names(INTERVENTIONS) <- sapply(INTERVENTIONS, `[[`, "name")

PRIMARY     <- names(INTERVENTIONS)[1:9]
SENSITIVITY <- names(INTERVENTIONS)[10:18]

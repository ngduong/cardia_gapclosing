# Allowability weights: odds weights that standardize White participants to
# the Black distribution of the allowable covariates (age and sex).  Black
# participants receive weight 1; weights are scaled to mean 1 within race.

allowability_weights <- function(dat, race = RACE, allowable = ALLOWABLE) {
  black <- as.integer(dat[[race]] == "Black")
  fit <- glm(as.formula(paste("black ~", paste(allowable, collapse = " + "))),
             data = cbind(dat, black = black), family = binomial())
  ps <- predict(fit, newdata = dat, type = "response")
  p_black <- mean(black)
  w <- ifelse(black == 1L, 1, ps / (1 - ps + 1e-12) * (1 - p_black) / p_black)
  dat$allow_weight <- w
  dat %>%
    group_by(across(all_of(race))) %>%
    mutate(allow_weight = allow_weight / mean(allow_weight)) %>%
    ungroup()
}

# Packages, analysis constants, data access and logging.

# Hmisc is called as Hmisc::wtd.quantile() rather than attached, since
# attaching it masks dplyr's summarize().
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(survival); library(eventglm); library(mgcv); library(splines)
  library(janitor); library(future); library(future.apply)
})

HORIZON        <- 35                                   # years of follow-up
TREATMENT      <- "le8_imputed"                        # Life's Essential 8 score, 0-100
RACE           <- "race_label"                         # "Black" / "White"
ALLOWABLE      <- c("examage", "sex")                  # standardized between groups
NON_ALLOWABLE  <- c("a03ed", "c03incom", "deprv80n00") # education, income, neighborhood deprivation
SCORE_BOUNDS   <- c(0, 100)
CENSORING_FORMULA <- ~ examage + sex + a03ed + deprv80n00 + c03incom +
                       a10smoke + a07mlalc + I(a07mlalc^2)

DATA_DIR <- "data"
BASE_DATA_FILE <- file.path(DATA_DIR, "cdm_data_clean_ipcw_imp_v2.csv")
LE8_DATA_FILE  <- file.path(DATA_DIR, "le8_score_all.csv")

# Baseline covariates and follow-up, one row per participant.  Columns derived
# from LE8 are dropped here and re-attached per imputation.
read_base_data <- function(file = BASE_DATA_FILE) {
  dat <- read_csv(file, show_col_types = FALSE)
  dat <- dat[, !grepl("^le8_", names(dat))]
  dat$mace_status <- factor(as.character(dat$mace_status),
                            levels = c("censoring", "MACE", "NonCVD death"))
  dat
}

# Multiply imputed LE8 scores at the baseline exam, long by imputation.
read_le8_imputations <- function(file = LE8_DATA_FILE, m_max = Inf) {
  read_csv(file, show_col_types = FALSE) %>%
    clean_names() %>%
    filter(exam == 1, imp <= m_max, !is.na(.data[[TREATMENT]])) %>%
    select(imp, id, all_of(TREATMENT))
}

attach_imputation <- function(base, le8_long, m) {
  base %>%
    select(-any_of(TREATMENT)) %>%
    left_join(filter(le8_long, imp == m) %>% select(-imp), by = "id")
}

log_line <- function(file, ...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(msg, "\n"); if (!is.null(file)) cat(msg, "\n", file = file, append = TRUE)
  flush.console()
}

# One multisession plan for the whole run; R caps open connections, so the
# worker count is reduced when fewer connections are free.
start_workers <- function(requested) {
  free  <- tryCatch(parallelly::freeConnections(), error = function(e) NA_integer_)
  cores <- tryCatch(parallelly::availableCores(), error = function(e) requested)
  w <- max(1L, min(requested, if (is.na(free)) requested else free - 12L, cores))
  plan(multisession, workers = w)
  w
}

# Runs every analysis in sequence.  From the repository folder in RStudio:
#   rstudioapi::jobRunScript("run_all.R", workingDir = getwd())
# or from a shell:  Rscript run_all.R
#
# Each analysis runs as a separate R process (run_analysis.R) and is restarted
# once if any intervention is missing.  Rerunning this script skips what has
# already been saved.  Progress: run_all.log, and results/<tag>/log.txt.

cores <- 118
B     <- 500
M     <- 50
analyses <- list(
  list(outcome = "censored",  interventions = "primary",     scale = "risk"),
  list(outcome = "competing", interventions = "primary",     scale = "risk"),
  list(outcome = "censored",  interventions = "sensitivity", scale = "risk"),
  list(outcome = "competing", interventions = "sensitivity", scale = "risk"),
  list(outcome = "censored",  interventions = "primary",     scale = "rmst"),
  list(outcome = "competing", interventions = "primary",     scale = "rmst"))
# Sensitivity to race-specific LE8 curves in the outcome model:
# analyses <- c(analyses, list(list(outcome = "censored",  interventions = "primary", scale = "risk", le8_curve = "by_race"),
#                              list(outcome = "competing", interventions = "primary", scale = "risk", le8_curve = "by_race")))

rscript <- file.path(R.home("bin"), "Rscript")
stamp <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  message(line); cat(line, "\n", file = "run_all.log", append = TRUE, sep = "")
}
stopifnot(file.exists("run_analysis.R"), dir.exists("R"))
stamp("start: ", cores, " workers, B = ", B, ", M = ", M)

check <- system(sprintf("%s check_installation.R > check_installation.log 2>&1", shQuote(rscript)))
stamp(tail(readLines("check_installation.log"), 1))
if (check != 0) stop("installation check failed; see check_installation.log")

incomplete <- character(0)
for (a in analyses) {
  curve <- if (is.null(a$le8_curve)) "common" else a$le8_curve
  tag <- sprintf("%s_%s_%s%s_B%d", a$outcome, a$interventions, a$scale, if (curve == "by_race") "_byrace" else "", B)
  cmd <- sprintf("%s run_analysis.R --outcome=%s --interventions=%s --scale=%s --le8_curve=%s --B=%d --M=%d --cores=%d >> %s 2>&1",
                 shQuote(rscript), a$outcome, a$interventions, a$scale, curve, B, M, cores,
                 shQuote(file.path("results", paste0(tag, "_console.log"))))
  dir.create("results", showWarnings = FALSE)
  n <- 0
  for (attempt in 1:2) {
    stamp("start ", tag, " (attempt ", attempt, ")")
    system(cmd)
    n <- length(list.files(file.path("results", tag), pattern = "\\.rds$"))
    stamp("end ", tag, ": ", n, " of 9 interventions saved")
    if (n >= 9) break
  }
  if (n < 9) incomplete <- c(incomplete, tag)
}
stamp(if (length(incomplete)) paste("finished with missing results:", paste(incomplete, collapse = ", "))
      else "all analyses complete")

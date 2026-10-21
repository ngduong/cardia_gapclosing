# Runs one analysis: a set of interventions for one outcome definition and
# one scale, with the bootstrap, across the imputations.  Results go to
# results/<tag>/, one .rds per intervention plus the bootstrap replicates.
# Interventions already saved are skipped, so an interrupted run can be resumed.
#
#   Rscript run_analysis.R --outcome=censored --interventions=primary --scale=risk
#
# Options (defaults in brackets):
#   --outcome        censored | competing      [censored]
#   --interventions  primary | sensitivity     [primary]
#   --scale          risk | rmst               [risk]
#   --le8_curve      common | by_race          [common]   LE8 term of the outcome model
#   --B              bootstrap replicates per imputation [500]
#   --M              imputations               [50]
#   --cores          parallel workers          [8]

for (f in list.files("R", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(".*=", "", hit[1]) else default
}
OUTCOME <- match.arg(opt("outcome", "censored"), c("censored", "competing"))
SET     <- match.arg(opt("interventions", "primary"), c("primary", "sensitivity"))
SCALE   <- match.arg(opt("scale", "risk"), c("risk", "rmst"))
CURVE   <- match.arg(opt("le8_curve", "common"), c("common", "by_race"))
B       <- as.integer(opt("B", 500)); M <- as.integer(opt("M", 50)); CORES <- as.integer(opt("cores", 8))

tag <- sprintf("%s_%s_%s%s_B%d", OUTCOME, SET, SCALE, if (CURVE == "by_race") "_byrace" else "", B)
out_dir <- file.path("results", tag); rep_dir <- file.path(out_dir, "replicates")
dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(out_dir, "log.txt")

base <- read_base_data()
le8  <- read_le8_imputations(m_max = M)
workers <- start_workers(CORES)
log_line(log_file, sprintf("%s: B = %d, M = %d, %d workers", tag, B, M, workers))

for (name in if (SET == "primary") PRIMARY else SENSITIVITY) {
  target <- file.path(out_dir, paste0(name, ".rds"))
  if (file.exists(target)) { log_line(log_file, "skip ", name); next }
  log_line(log_file, "start ", name); t0 <- Sys.time()
  res <- tryCatch(
    run_imputations(base, le8, INTERVENTIONS[[name]], OUTCOME, SCALE, CURVE == "by_race",
                    M = M, B = B, replicate_dir = rep_dir, log_file = log_file),
    error = function(e) { log_line(log_file, sprintf("failed %s: %s", name, conditionMessage(e))); NULL })
  if (!is.null(res)) {
    res$settings <- list(outcome = OUTCOME, scale = SCALE, le8_curve = CURVE, B = B, M = M,
                         truncation = c(0.001, 0.999), finished = format(Sys.time()))
    saveRDS(res, target)
  }
  log_line(log_file, sprintf("%s %s (%.1f min)", if (is.null(res)) "failed" else "done", name,
                             as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
log_line(log_file, "finished")

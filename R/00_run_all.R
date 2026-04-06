# 00_run_all.R
# Master script — toggle run_XX flags and set pipeline_config, then source.
# Figure scripts (07-11) are hardcoded to MODE = "core"; run individually
# or via the toggles below.


# Configuration

MODE <- "core"  # "core" or "LOO"

pipeline_config <- list(
  mode        = MODE,
  antigens    = c("K", "O", "KO"),    # Any subset of c("K", "O", "KO")
  count_types = c("adj", "raw"),       # Any subset of c("adj", "raw")
  subsets     = NULL                    # NULL = all; or c("Full", "Carba", ...)
)

# Script toggles
run_01_model_run       <- FALSE
run_02_model_process   <- FALSE


message("\n", rep("=", 70))
message("KNNS WEIGHTING ANALYSIS PIPELINE")
message(rep("=", 70))
message("Mode:        ", MODE)
message("Antigens:    ", paste(pipeline_config$antigens, collapse = ", "))
message("Count types: ", paste(pipeline_config$count_types, collapse = ", "))
message("Subsets:     ", if (is.null(pipeline_config$subsets)) "ALL"
                         else paste(pipeline_config$subsets, collapse = ", "))
message("Started:     ", Sys.time())
message(rep("=", 70), "\n")

if (run_01_model_run) {
  message("\n>>> Running 01_model_run.R <<<\n")
  source("R/01_model_run.R")
  diagnostics_01 <- run_model_fitting(pipeline_config)
  message("\n>>> Completed 01_model_run.R <<<\n")

  if (any(!diagnostics_01$diagnostics_pass)) {
    message("WARNING: Some models failed diagnostics. Review diagnostics_01.")
  }
}

if (run_02_model_process) {
  message("\n>>> Running 02_model_process.R <<<\n")
  source("R/02_model_process.R")
  run_model_processing(pipeline_config)
  message("\n>>> Completed 02_model_process.R <<<\n")
}


message("\n", rep("=", 70))
message("PIPELINE COMPLETE")
message(rep("=", 70))
message("Mode: ", MODE)
message("Finished: ", Sys.time())
message(rep("=", 70), "\n")

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
run_03_weighing        <- FALSE
run_04_model_plots     <- FALSE
run_05_diagnostics     <- FALSE
run_06_weight_analysis <- FALSE

run_07_ranked_density    <- FALSE
run_08_KO_density        <- FALSE
run_09_prevalence_panel  <- FALSE
run_10_comparison_panel  <- FALSE
run_11_stanton           <- FALSE


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

if (run_03_weighing) {
  message("\n>>> Running 03_weighing.R <<<\n")
  source("R/03_weighing.R")
  message("\n>>> Completed 03_weighing.R <<<\n")
}

if (run_04_model_plots) {
  message("\n>>> Running 04_model_plots.R <<<\n")
  source("R/04_model_plots.R")
  message("\n>>> Completed 04_model_plots.R <<<\n")
}

if (run_05_diagnostics) {
  message("\n>>> Running 05_diagnostics_table.R <<<\n")
  source("R/05_diagnostics_table.R")
  message("\n>>> Completed 05_diagnostics_table.R <<<\n")
}

if (run_06_weight_analysis) {
  message("\n>>> Running 06_weight_analysis.R <<<\n")
  source("R/06_weight_analysis.R")
  message("\n>>> Completed 06_weight_analysis.R <<<\n")
}

if (run_07_ranked_density) {
  message("\n>>> Running 07_ranked_density_plots.R <<<\n")
  source("R/07_ranked_density_plots.R")
  message("\n>>> Completed 07_ranked_density_plots.R <<<\n")
}

if (run_08_KO_density) {
  message("\n>>> Running 08_KO_density_plots.R <<<\n")
  source("R/08_KO_density_plots.R")
  message("\n>>> Completed 08_KO_density_plots.R <<<\n")
}

if (run_09_prevalence_panel) {
  message("\n>>> Running 09_prevalence_panel_plot.R <<<\n")
  source("R/09_prevalence_panel_plot.R")
  message("\n>>> Completed 09_prevalence_panel_plot.R <<<\n")
}

if (run_10_comparison_panel) {
  message("\n>>> Running 10_comparison_panel_plot.R <<<\n")
  source("R/10_comparison_panel_plot.R")
  message("\n>>> Completed 10_comparison_panel_plot.R <<<\n")
}

if (run_11_stanton) {
  message("\n>>> Running 11_stanton_comparison.R <<<\n")
  source("R/11_stanton_comparison.R")
  message("\n>>> Completed 11_stanton_comparison.R <<<\n")
}


message("\n", rep("=", 70))
message("PIPELINE COMPLETE")
message(rep("=", 70))
message("Mode: ", MODE)
message("Finished: ", Sys.time())
message(rep("=", 70), "\n")

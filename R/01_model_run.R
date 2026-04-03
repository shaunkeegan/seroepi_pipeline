# 01_model_run.R
# Fit hierarchical Bayesian binomial models for K/O locus prevalence.
# Output: brms model .rds files + JSON sidecars in models/{core|LOO}/

library(tidyverse)
library(brms)
library(jsonlite)

source("R/utils.R")


#' Fit Bayesian binomial models for all configured antigen/count_type combos.
#'
#' @param config Named list: mode ("core" or "LOO"), antigens (character
#'   vector), count_types (character vector), subsets (NULL or character
#'   vector).
#' @return Data frame of diagnostics for all fitted models (invisibly).
run_model_fitting <- function(config) {

  dirs <- resolve_directories(config$mode)

  message("=== MODEL FITTING: ", toupper(config$mode), " MODE ===")
  message("Input:  ", dirs$data_dir)
  message("Output: ", dirs$model_dir, "\n")

  if (!dir.exists(dirs$model_dir)) {
    dir.create(dirs$model_dir, recursive = TRUE)
  }

  file_counter <- get_max_counter(dirs$model_dir)
  message("File counter initialised at: ", file_counter)

  all_files <- list_site_count_files(dirs$data_dir)

  message("Config: antigens = [", paste(config$antigens, collapse = ", "), "]",
          " | count_types = [", paste(config$count_types, collapse = ", "), "]",
          " | subsets = ",
          if (is.null(config$subsets)) "ALL" else paste0("[", paste(config$subsets, collapse = ", "), "]"))

  diagnostics_log <- data.frame(
    file_id               = character(),
    antigen               = character(),
    subset                = character(),
    count_type            = character(),
    max_rhat              = numeric(),
    min_neff_ratio        = numeric(),
    divergent_transitions = integer(),
    diagnostics_pass      = logical(),
    stringsAsFactors      = FALSE
  )

  for (i in seq_along(all_files)) {

    filename <- all_files[i]
    meta <- extract_file_metadata(filename, config$mode)

    if (!meta$antigen %in% config$antigens) {
      message("Skipping ", filename, " (antigen ", meta$antigen, " not in config)")
      next
    }
    if (!is.null(config$subsets) && !meta$subset %in% config$subsets) {
      message("Skipping ", filename, " (subset ", meta$subset, " not in config)")
      next
    }

    message("\n", rep("=", 60))
    message("Processing file ", i, " of ", length(all_files), ": ", filename)
    message(rep("=", 60))
    message("  Subset: ", meta$subset, " | Antigen: ", meta$antigen,
            " | Days: ", meta$days, " | Data: ", meta$what_data)
    if (!is.na(meta$study)) message("  LOO Study: ", meta$study)

    site_prevalence <- read.delim(
      file.path(dirs$data_dir, filename), sep = "\t"
    )
    message("  Loaded ", nrow(site_prevalence), " rows, ",
            n_distinct(site_prevalence$locus), " loci, ",
            n_distinct(site_prevalence$Site), " sites")

    for (type in config$count_types) {

      message("\n  --- Fitting ", toupper(type), " model ---")

      data_csv <- prepare_model_data(site_prevalence, type)
      message("  Data prepared: ", nrow(data_csv), " rows")

      fit_with_interaction <- brm(
        formula = bf(
          event | trials(n) ~ 0 + locus + subgroup +
            (1 | Site) +
            (1 | locus:subgroup)
        ),
        data = data_csv,
        family = binomial(link = "logit"),
        prior = c(
          prior(student_t(3, 0, 5), class = "sd", group = "Site"),
          prior(student_t(3, 0, 5), class = "sd", group = "locus:subgroup")
        ),
        iter = 6000,
        chains = 4,
        cores = 4,
        control = list(adapt_delta = 0.999999, max_treedepth = 55)
      )

      diagnostics <- extract_model_diagnostics(fit_with_interaction)

      if (!diagnostics$diagnostics_pass) {
        warning("Model diagnostics FAILED for ", meta$antigen, " ",
                meta$subset, " ", type,
                ": rhat=", round(diagnostics$max_rhat, 4),
                ", neff=", round(diagnostics$min_neff_ratio, 4),
                ", divergent=", diagnostics$divergent_transitions)
      } else {
        message("  Diagnostics PASSED (rhat=",
                round(diagnostics$max_rhat, 4),
                ", neff=", round(diagnostics$min_neff_ratio, 4),
                ", divergent=", diagnostics$divergent_transitions, ")")
      }

      file_counter   <- file_counter + 1L
      file_id        <- make_file_id(meta$date, file_counter)
      model_filepath <- file.path(dirs$model_dir, paste0(file_id, ".rds"))

      saveRDS(fit_with_interaction, file = model_filepath)

      write_json_sidecar(model_filepath, c(
        list(
          file_type      = "model",
          antigen        = meta$antigen,
          subset         = meta$subset,
          what_data      = meta$what_data,
          count_type     = type,
          days           = as.integer(meta$days),
          data_date      = meta$date,
          held_out_study = if (is.na(meta$study)) NULL else meta$study,
          source_file    = filename
        ),
        diagnostics
      ))

      message("  Saved: ", model_filepath, " [", file_id, "]")

      diagnostics_log <- rbind(diagnostics_log, data.frame(
        file_id               = file_id,
        antigen               = meta$antigen,
        subset                = meta$subset,
        count_type            = type,
        max_rhat              = diagnostics$max_rhat,
        min_neff_ratio        = diagnostics$min_neff_ratio,
        divergent_transitions = diagnostics$divergent_transitions,
        diagnostics_pass      = diagnostics$diagnostics_pass,
        stringsAsFactors      = FALSE
      ))
    }
  }

  message("\n", rep("=", 60))
  message("MODEL FITTING COMPLETE")
  message(rep("=", 60))
  message("Models fitted: ", nrow(diagnostics_log))
  message("Output: ", dirs$model_dir)

  failed <- diagnostics_log[!diagnostics_log$diagnostics_pass, ]
  if (nrow(failed) > 0) {
    message("\nWARNING: ", nrow(failed), " of ", nrow(diagnostics_log),
            " models failed diagnostics:")
    print(failed)
  } else if (nrow(diagnostics_log) > 0) {
    message("\nAll ", nrow(diagnostics_log), " models passed diagnostics.")
    message("  - All R-hat values <= 1.01")
    message("  - All effective sample size ratios >= 0.1")
    message("  - No divergent transitions")
  }

  invisible(diagnostics_log)
}


if (!exists("pipeline_config")) {
  if (!exists("MODE")) MODE <- "core"
  run_model_fitting(list(
    mode        = MODE,
    antigens    = c("K", "O", "KO"),
    count_types = c("adj", "raw"),
    subsets     = NULL
  ))
}

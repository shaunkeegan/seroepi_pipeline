# 02_model_process.R
# Extract posterior distributions from fitted brms models.
# Input:  models/{core|LOO}/
# Output: outputs/posteriors/{core|LOO}/, outputs/estimates/{core|LOO}/

library(dplyr)
library(stringr)
library(tidyr)
library(brms)
library(jsonlite)

source("R/utils.R")


#' Extract posteriors from all configured models.
#'
#' @param config Named list: mode ("core" or "LOO"), antigens (character
#'   vector), count_types (character vector), subsets (NULL or character
#'   vector).
#' @return NULL (invisibly). Outputs saved to disk.
run_model_processing <- function(config) {

  dirs <- resolve_directories(config$mode, outputs = TRUE)

  message("=== POSTERIOR EXTRACTION: ", toupper(config$mode), " MODE ===")
  message("Input:  ", dirs$data_dir)
  message("Models: ", dirs$model_dir)
  message("Output: ", dirs$estimates_dir, " & ", dirs$posteriors_dir, "\n")

  if (!dir.exists(dirs$estimates_dir))  dir.create(dirs$estimates_dir, recursive = TRUE)
  if (!dir.exists(dirs$posteriors_dir)) dir.create(dirs$posteriors_dir, recursive = TRUE)

  post_counter <- get_max_counter(dirs$posteriors_dir)
  est_counter  <- get_max_counter(dirs$estimates_dir)
  message("Posterior counter initialised at: ", post_counter)
  message("Estimates counter initialised at: ", est_counter)

  all_files <- list_site_count_files(dirs$data_dir)

  message("Config: antigens = [", paste(config$antigens, collapse = ", "), "]",
          " | count_types = [", paste(config$count_types, collapse = ", "), "]",
          " | subsets = ",
          if (is.null(config$subsets)) "ALL" else paste0("[", paste(config$subsets, collapse = ", "), "]"))

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

    for (type in config$count_types) {

      message("\n  --- Processing ", toupper(type), " ---")

      model_id <- find_file_by_metadata(
        dirs$model_dir,
        file_type      = "model",
        antigen        = meta$antigen,
        subset         = meta$subset,
        what_data      = meta$what_data,
        count_type     = type,
        days           = as.integer(meta$days),
        data_date      = meta$date,
        held_out_study = if (is.na(meta$study)) NULL else meta$study
      )

      model_filepath <- file.path(dirs$model_dir, paste0(model_id, ".rds"))
      fit_with_interaction <- readRDS(model_filepath)
      message("  Loaded model: ", model_id)

      site_prevalence <- read.delim(file.path(dirs$data_dir, filename), sep = "\t")
      data_csv <- prepare_model_data(site_prevalence, type)
      message("  Data prepared: ", nrow(data_csv), " rows")

      message("  Extracting posteriors...")

      post_epred <- posterior_epred(
        fit_with_interaction,
        newdata = data_csv,
        re_formula = NULL,
        allow_new_levels = TRUE,
        summary = FALSE
      )

      message("  Computing global posteriors...")
      posterior_locus_df <- compute_posterior_prevalence(
        post_epred, data_csv, grouping_vars = "locus"
      )

      global_estimates <- posterior_locus_df %>%
        group_by(locus) %>%
        summarise(
          mean   = mean(prob),
          median = median(prob),
          lower  = quantile(prob, 0.025, na.rm = TRUE),
          upper  = quantile(prob, 0.975, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(subgroup = "Global")

      message("  Computing regional posteriors...")
      posterior_locus_subgroup_df <- compute_posterior_prevalence(
        post_epred, data_csv, grouping_vars = c("locus", "subgroup")
      )

      subgroup_estimates <- posterior_locus_subgroup_df %>%
        group_by(locus, subgroup) %>%
        summarise(
          mean   = mean(prob),
          median = median(prob),
          lower  = quantile(prob, 0.025, na.rm = TRUE),
          upper  = quantile(prob, 0.975, na.rm = TRUE),
          .groups = "drop"
        )

      shared_meta <- list(
        antigen        = meta$antigen,
        subset         = meta$subset,
        what_data      = meta$what_data,
        count_type     = type,
        days           = as.integer(meta$days),
        data_date      = meta$date,
        held_out_study = if (is.na(meta$study)) NULL else meta$study,
        source_file    = filename,
        model_file     = model_id
      )

      post_counter     <- post_counter + 1L
      post_global_id   <- make_file_id(meta$date, post_counter)
      post_global_path <- file.path(dirs$posteriors_dir, paste0(post_global_id, ".csv.gz"))
      write.csv(posterior_locus_df, gzfile(post_global_path), row.names = FALSE)
      write_json_sidecar(post_global_path,
        c(list(file_type = "posterior", level = "global"), shared_meta))

      est_counter     <- est_counter + 1L
      est_global_id   <- make_file_id(meta$date, est_counter)
      est_global_path <- file.path(dirs$estimates_dir, paste0(est_global_id, ".csv"))
      write.csv(global_estimates, est_global_path, row.names = FALSE)
      write_json_sidecar(est_global_path,
        c(list(file_type = "estimates", level = "global"), shared_meta))

      post_counter       <- post_counter + 1L
      post_subgroup_id   <- make_file_id(meta$date, post_counter)
      post_subgroup_path <- file.path(dirs$posteriors_dir, paste0(post_subgroup_id, ".csv.gz"))
      write.csv(posterior_locus_subgroup_df, gzfile(post_subgroup_path), row.names = FALSE)
      write_json_sidecar(post_subgroup_path,
        c(list(file_type = "posterior", level = "subgroup"), shared_meta))

      est_counter        <- est_counter + 1L
      est_subgroup_id    <- make_file_id(meta$date, est_counter)
      est_subgroup_path  <- file.path(dirs$estimates_dir, paste0(est_subgroup_id, ".csv"))
      write.csv(subgroup_estimates, est_subgroup_path, row.names = FALSE)
      write_json_sidecar(est_subgroup_path,
        c(list(file_type = "estimates", level = "subgroup"), shared_meta))

      message("  Saved: posterior global [", post_global_id, "]",
              " | estimates global [", est_global_id, "]")
      message("         posterior subgroup [", post_subgroup_id, "]",
              " | estimates subgroup [", est_subgroup_id, "]")
    }
  }

  message("\n", rep("=", 60))
  message("POSTERIOR EXTRACTION COMPLETE")
  message(rep("=", 60))
  message("Posteriors: ", dirs$posteriors_dir)
  message("Estimates:  ", dirs$estimates_dir)

  invisible(NULL)
}


if (!exists("pipeline_config")) {
  if (!exists("MODE")) MODE <- "core"
  run_model_processing(list(
    mode        = MODE,
    antigens    = c("K", "O", "KO"),
    count_types = c("adj", "raw"),
    subsets     = NULL
  ))
}

# utils.R
# Utility functions for opaque ID management and JSON sidecar read/write.

library(jsonlite)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)


# ID management

#' Scan a directory for existing opaque-ID files and return the highest counter.
#'
#' Files are expected to follow the pattern YYYYMMDD-NNNNNN.<ext>.
#' Returns 0L if no matching files exist (clean directory).
#'
#' @param dir  Directory to scan.
#' @return     Integer: the highest counter found, or 0L if none.
get_max_counter <- function(dir) {
  existing <- list.files(dir, pattern = "^\\d{8}-\\d{6}\\.")
  if (length(existing) == 0) return(0L)
  nums <- suppressWarnings(
    as.integer(sub("^\\d{8}-(\\d{6})\\..*", "\\1", existing))
  )
  max(nums[!is.na(nums)], 0L)
}

#' Construct an opaque file ID string.
#'
#' @param date_prefix  8-digit date string, e.g. "20251205".
#' @param counter      Integer counter value.
#' @return             Character string, e.g. "20251205-000001".
make_file_id <- function(date_prefix, counter) {
  sprintf("%s-%06d", date_prefix, as.integer(counter))
}


# JSON sidecar writing

#' Write a JSON metadata sidecar alongside a data file.
#'
#' The sidecar is written to the same directory as `filepath`, with the same
#' base name and a .json extension. For .csv.gz files, both extensions are
#' stripped before appending .json.
#'
#' @param filepath  Full path to the data file.
#' @param metadata  Named list of metadata fields.
#' @return          Path to the written JSON file (invisibly).
write_json_sidecar <- function(filepath, metadata) {
  json_path <- if (grepl("\\.csv\\.gz$", filepath)) {
    sub("\\.csv\\.gz$", ".json", filepath)
  } else {
    paste0(tools::file_path_sans_ext(filepath), ".json")
  }

  writeLines(
    jsonlite::toJSON(metadata, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    json_path
  )

  invisible(json_path)
}


# JSON sidecar querying

#' Find exactly one file in a directory matching all supplied metadata filters.
#'
#' Queries all .json sidecar files in `dir` and returns the base file ID
#' (YYYYMMDD-NNNNNN, without extension) of the unique matching file.
#' Stops with an informative error if zero or multiple matches are found.
#'
#' @param dir   Directory to search.
#' @param ...   Named filter arguments, e.g. antigen = "K", file_type = "model".
#' @return      Character: base file ID of the matching file.
find_file_by_metadata <- function(dir, ...) {
  results <- find_all_files_by_metadata(dir, ...)
  filters <- list(...)

  filter_str <- paste(names(filters), unlist(lapply(filters, function(x)
    if (is.null(x)) "NULL" else as.character(x)
  )), sep = "=", collapse = ", ")

  if (length(results) == 0) {
    stop("No file found in '", dir, "' matching: ", filter_str)
  }
  if (length(results) > 1) {
    stop("Multiple files (", length(results), ") found in '", dir,
         "' matching: ", filter_str,
         "\nMatches: ", paste(results, collapse = ", "))
  }

  results[[1]]
}

#' Find all files in a directory matching all supplied metadata filters.
#'
#' Queries all .json sidecar files in `dir` and returns base file IDs
#' (YYYYMMDD-NNNNNN, without extension) for all matching files.
#'
#' @param dir   Directory to search.
#' @param ...   Named filter arguments, e.g. file_type = "posterior".
#' @return      Character vector of matching base file IDs (may be empty).
find_all_files_by_metadata <- function(dir, ...) {
  filters <- list(...)
  json_files <- list.files(dir, pattern = "^\\d{8}-\\d{6}\\.json$",
                           full.names = TRUE)

  if (length(json_files) == 0) return(character(0))

  matches <- character(0)
  for (jf in json_files) {
    meta <- tryCatch(jsonlite::read_json(jf), error = function(e) NULL)
    if (is.null(meta)) next
    if (.metadata_matches(meta, filters)) {
      matches <- c(matches, sub("\\.json$", "", basename(jf)))
    }
  }

  matches
}

#' Internal: check whether a metadata list satisfies all filter criteria.
#'
#' NULL filter values match JSON null (R NULL). Comparisons are coerced to
#' character to handle integer/string differences (e.g. days: 28 vs "28").
#'
#' @param meta     Named list from jsonlite::read_json().
#' @param filters  Named list of filter values.
#' @return         Logical scalar.
.metadata_matches <- function(meta, filters) {
  for (key in names(filters)) {
    filter_val <- filters[[key]]
    meta_val   <- meta[[key]]

    if (is.null(filter_val) && is.null(meta_val)) next
    if (is.null(filter_val) || is.null(meta_val)) return(FALSE)
    # coerce to character - JSON integers arrive as numeric
    if (!identical(as.character(meta_val), as.character(filter_val))) return(FALSE)
  }
  TRUE
}


# Manifest regeneration

#' Regenerate manifest.csv in a directory from all JSON sidecar files.
#'
#' The manifest is a convenience view only - the JSON sidecars are the source
#' of truth. This function can be called at any time to rebuild or refresh it.
#' Skips directories containing no JSON sidecar files.
#'
#' @param dir  Directory to process.
#' @return     The manifest data frame (invisibly), or NULL if skipped.
rebuild_manifest <- function(dir) {
  if (!dir.exists(dir)) {
    message("Directory does not exist: ", dir, " - skipping")
    return(invisible(NULL))
  }

  json_files <- list.files(dir, pattern = "^\\d{8}-\\d{6}\\.json$",
                           full.names = TRUE)

  if (length(json_files) == 0) {
    message("No JSON sidecars found in: ", dir, " - skipping")
    return(invisible(NULL))
  }

  manifest <- purrr::map_dfr(json_files, function(jf) {
    meta <- jsonlite::read_json(jf)
    # NULL fields break data.frame coercion
    meta[vapply(meta, is.null, logical(1))] <- NA_character_
    as.data.frame(
      c(list(filename = sub("\\.json$", "", basename(jf))), meta),
      stringsAsFactors = FALSE
    )
  })

  out_path <- file.path(dir, "manifest.csv")
  write.csv(manifest, out_path, row.names = FALSE)
  message("Manifest written: ", out_path, " (", nrow(manifest), " entries)")
  invisible(manifest)
}


# Directory resolution

#' Build standard directory paths for a given analysis mode.
#'
#' @param mode    Character: "core" or "LOO".
#' @param outputs Logical: if TRUE, also return estimates_dir and posteriors_dir.
#' @return Named list: data_dir, model_dir [, estimates_dir, posteriors_dir].
resolve_directories <- function(mode, outputs = FALSE) {
  if (!mode %in% c("core", "LOO")) {
    stop("Invalid mode '", mode, "'. Must be 'core' or 'LOO'.")
  }

  dirs <- list(
    data_dir  = if (mode == "core") "data/core/site_counts/" else "data/LOO/",
    model_dir = if (mode == "core") "models/core/"           else "models/LOO/"
  )

  if (outputs) {
    dirs$estimates_dir  <- paste0("outputs/estimates/",  mode, "/")
    dirs$posteriors_dir <- paste0("outputs/posteriors/", mode, "/")
  }

  dirs
}


# File listing and validation

#' List and validate site count files in a directory.
#'
#' Lists all files, excludes consolidated summaries and site_info files,
#' validates against the expected filename regex, and warns on mismatches.
#'
#' @param data_dir Character: path to site count directory.
#' @return Character vector of validated filenames (relative to data_dir).
list_site_count_files <- function(data_dir) {
  all_files <- list.files(path = data_dir, recursive = TRUE)

  all_files <- all_files[!grepl("filterN10_all\\.tsv$", all_files)]
  all_files <- all_files[!grepl("^site_info", all_files)]

  message("Found ", length(all_files), " files in ", data_dir)

  pattern <- paste0(
    "^((Carba|ESBL|Fatal)_)?",
    "(min10_)?",
    "Neonatal_shareable_\\d{8}_\\d+_\\d+_filterN10_",
    "(K|O|KO)_site_counts_all\\.tsv$"
  )

  bad_files <- all_files[!grepl(pattern, all_files)]

  if (length(bad_files) > 0) {
    warning("Files not matching required format:\n",
            paste("  - ", bad_files, collapse = "\n"))
  } else {
    message("All files match the required naming format.")
  }

  all_files
}


# Filename metadata extraction

#' Extract metadata fields from a site count filename.
#'
#' Parses subset, date, days, what_data, antigen, and study from filename.
#'
#' @param filename Character: the filename to parse.
#' @param mode     Character: "core" or "LOO" - controls LOO study extraction.
#' @return Named list: subset, date, days, what_data, antigen, study.
extract_file_metadata <- function(filename, mode) {
  subset    <- str_extract(filename, "^(Carba|ESBL|Fatal)")
  subset    <- ifelse(is.na(subset), "Full", subset)
  date      <- str_extract(filename, "\\d{8}")
  days      <- str_extract(filename, "(?<=\\d{8}_)\\d+")
  what_data <- ifelse(str_detect(filename, "min10"), "min10", "ALL")
  antigen   <- str_extract(filename, "_(K|O|KO)_") %>% str_replace_all("_", "")

  study <- NA_character_
  if (mode == "LOO" && str_detect(filename, "LOO")) {
    study <- str_extract(filename, "(?<=LOO_)[^.]+")
  }

  list(
    subset    = subset,
    date      = date,
    days      = days,
    what_data = what_data,
    antigen   = antigen,
    study     = study
  )
}


# Data preparation

#' Prepare site prevalence data for brms model fitting or posterior extraction.
#'
#' Selects columns for the given count type, renames to model-standard names,
#' expands via expand.grid to include zero-event locus-site combinations,
#' joins n-values and events, replaces NA events with 0, and drops rows
#' with NA n.
#'
#' @param site_prevalence Data frame: raw site count data as loaded from TSV.
#' @param count_type      Character: "adj" or "raw".
#' @return Tibble ready for brm() or posterior_epred().
prepare_model_data <- function(site_prevalence, count_type) {
  site_prevalence_type <- site_prevalence %>%
    select(locus, Site, Region, contains(count_type)) %>%
    rename(
      subgroup = Region,
      event    = paste0(count_type, "_count"),
      n        = paste0(count_type, "_sum")
    )

  site_subgroup_n <- site_prevalence_type %>%
    distinct(subgroup, Site, n)

  all_combinations <- expand.grid(
    locus    = unique(site_prevalence_type$locus),
    subgroup = unique(site_prevalence_type$subgroup),
    Site     = unique(site_prevalence_type$Site),
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    as_tibble()

  all_combinations %>%
    left_join(site_subgroup_n, by = c("subgroup", "Site")) %>%
    left_join(
      site_prevalence_type %>% select(subgroup, Site, locus, event),
      by = c("subgroup", "Site", "locus")
    ) %>%
    mutate(event = if_else(is.na(event), 0, event)) %>%
    tidyr::drop_na(n)
}


# Model diagnostics

#' Extract MCMC diagnostics from a fitted brms model.
#'
#' @param fit A brmsfit object.
#' @return Named list: max_rhat, min_neff_ratio, divergent_transitions,
#'         diagnostics_pass (logical).
extract_model_diagnostics <- function(fit) {
  max_rhat <- max(rhat(fit), na.rm = TRUE)
  min_neff <- min(neff_ratio(fit), na.rm = TRUE)

  # nuts_params() preferred; fall back to internal slot if unavailable
  div_trans <- tryCatch({
    np <- nuts_params(fit)
    sum(np$Value[np$Parameter == "divergent__"])
  }, error = function(e) {
    sum(fit$fit@sim$sampler_diagnostics$divergent__)
  })

  pass <- (max_rhat <= 1.01) && (min_neff >= 0.1) && (div_trans == 0)

  list(
    max_rhat              = max_rhat,
    min_neff_ratio        = min_neff,
    divergent_transitions = as.integer(div_trans),
    diagnostics_pass      = pass
  )
}


# Posterior prevalence computation

#' Compute posterior prevalence from posterior expected predictions.
#'
#' Vectorised computation using rowSums instead of nested for-loops.
#'
#' @param post_epred     Matrix: draws x observations from posterior_epred().
#' @param data_csv       Data frame: the prepared data (with locus, subgroup, n).
#' @param grouping_vars  Character vector: columns to group by.
#'        "locus" for global; c("locus", "subgroup") for regional.
#' @return Data frame: locus, [subgroup], draw_id, prob.
compute_posterior_prevalence <- function(post_epred, data_csv, grouping_vars) {
  n_draws <- nrow(post_epred)
  groups  <- data_csv %>% distinct(across(all_of(grouping_vars)))

  result_list <- vector("list", nrow(groups))

  for (k in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(data_csv))
    for (g in grouping_vars) {
      idx <- idx & (data_csv[[g]] == groups[[g]][k])
    }

    total_n          <- sum(data_csv$n[idx], na.rm = TRUE)
    prevalence_draws <- rowSums(post_epred[, idx, drop = FALSE]) / total_n

    row <- data.frame(
      draw_id = seq_len(n_draws),
      prob    = prevalence_draws,
      stringsAsFactors = FALSE
    )
    for (g in grouping_vars) {
      row[[g]] <- groups[[g]][k]
    }

    result_list[[k]] <- row
  }

  result <- do.call(rbind, result_list)
  rownames(result) <- NULL

  if (!"subgroup" %in% grouping_vars) {
    result$subgroup <- "Global"
  }

  col_order <- c("locus", "subgroup", "draw_id", "prob")
  col_order <- col_order[col_order %in% names(result)]
  result[, col_order]
}

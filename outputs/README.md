# outputs/

Pipeline outputs from `02_model_process.R`. All output files (`.csv`, `.rds`, `.json`) are excluded from version control; this README preserves the directory structure in git.

## Subdirectories

| Directory | Contents |
|-----------|----------|
| `estimates/core/` | Summary posterior statistics for the full dataset |
| `estimates/LOO/` | Summary posterior statistics for each leave-one-out model |
| `posteriors/core/` | Full posterior draw files for the full dataset |
| `posteriors/LOO/` | Full posterior draw files for each leave-one-out model |

## File identity system

Output files use opaque timestamped IDs (`YYYYMMDD-NNNNNN`) rather than descriptive names. Each output file has a paired JSON sidecar that stores all metadata (antigen, subset, count type, mode, etc.).

To locate a file by metadata, use `find_file_by_metadata()` from `R/utils.R`.

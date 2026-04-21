# seroepi_pipeline

Bayesian hierarchical pipeline for estimating *Klebsiella pneumoniae* serotype (K/O locus) prevalence from multi-site genomic surveillance data.

This pipeline is an update of the analysis code from [Stanton, Keegan et al. 2026](https://github.com/klebgenomics/KlebSero_paper), with updates including improved metadata handling. Please cite the original paper if you use this pipeline:

> Stanton TD, Keegan SP, Abdulahi JA, et al. Distribution of capsule and O types in *Klebsiella pneumoniae* causing neonatal sepsis in Africa and South Asia: A meta-analysis of genome-predicted serotype prevalence to inform potential vaccine coverage. *PLoS Med.* 2026;23(1):e1004879. doi:10.1371/journal.pmed.1004879

---

## Overview

The pipeline fits binomial mixed-effects models using [brms](https://paul-buerkner.github.io/brms/) to estimate the prevalence of each K/O locus globally and by region. It supports two analysis modes:

- **core** -- full multi-site dataset
- **LOO** -- leave-one-out cross-validation (one study held out per run)

---

## Requirements

R packages required:

- `brms`
- `tidyverse`
- `jsonlite`

No formal dependency management file is included. Install packages manually before running.

---

## Input data

Place site count files in the appropriate directory:

| Mode | Directory |
|------|-----------|
| core | `data/core/site_counts/` |
| LOO  | `data/LOO/` |

**core files** are tab-separated with the naming convention:

```
[Subset_][min10_]Neonatal_shareable_YYYYMMDD_DD_N_filterN10_[K|O|KO]_site_counts_all.tsv
```

**LOO files** follow the same convention with a `_LOO_[StudyID]` suffix and use CSV format.

Expected columns: `locus`, `Site`, `Country`, `Region`, `raw_count`, `adj_count`, `raw_sum`, `adj_sum`, `raw_prop`, `adj_prop`.

---

## Running the pipeline

Open [R/00_run_all.R](R/00_run_all.R) and set the configuration at the top of the file:

```r
MODE        <- "core"          # "core" or "LOO"
antigens    <- c("K", "O")    # antigens to model
count_types <- "adj"           # "adj" (KNNS-adjusted) or "raw"
subsets     <- NULL            # NULL for all, or e.g. c("Carba", "ESBL", "Fatal")

run_01_model_run     <- TRUE   # fit models
run_02_model_process <- TRUE   # extract posteriors and summaries
```

Then source the file:

```r
source("R/00_run_all.R")
```

The two stages can be run independently by toggling the boolean flags.

---

## Outputs

All outputs are written under `outputs/` with an opaque identifier (`YYYYMMDD-NNNNNN`). Each output file has a JSON sidecar that records its metadata (antigen, subset, count type, analysis mode, convergence diagnostics, etc.).

| Path | Contents |
|------|----------|
| `outputs/posteriors/{core\|LOO}/` | Full posterior draws (`.csv.gz`), one row per draw |
| `outputs/estimates/{core\|LOO}/`  | Summary statistics per locus: mean, median, 95% credible interval |

Summaries are produced at two levels: global (all sites pooled) and by subgroup (region).

---

## Model specification

The model is fit with `brms::brm()`:

```
event | trials(n) ~ 0 + locus + subgroup + (1 | Site) + (1 | locus:subgroup)
```

- Family: binomial (logit link)
- Priors: Student-t(3, 0, 5) on random-effect standard deviations
- MCMC: 6000 iterations, 4 chains
- Convergence criteria: R-hat <= 1.01, ESS ratio >= 0.1, zero divergent transitions

---

## Repository structure

```
R/
  00_run_all.R        entry point and configuration
  01_model_run.R      model fitting
  02_model_process.R  posterior extraction and summarisation
  utils.R             shared utility functions
data/
  core/site_counts/   full dataset input files
  LOO/                leave-one-out input files
outputs/
  posteriors/         posterior draw files
  estimates/          summary statistic files
```

---

## License

MIT License. See [LICENSE](LICENSE).

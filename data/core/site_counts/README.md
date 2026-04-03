# data/core/site_counts/

Per-locus site count TSVs for the full (`core`) analysis. One file per antigen × clinical subset × count variant combination. These files are excluded from version control.

## File naming convention

```
[Subset_][min10_]Neonatal_shareable_YYYYMMDD_DD_N_filterN10_[K|O|KO]_site_counts_all.tsv
```

| Component | Values | Meaning |
|-----------|--------|---------|
| `Subset_` | `Carba_`, `ESBL_`, `Fatal_`, *(absent)* | Clinical subset; absent = Full dataset |
| `min10_` | present or absent | `min10` variant filters loci with fewer than 10 observations at a site; absent = all data |
| `YYYYMMDD` | e.g. `20251205` | Data freeze date |
| `DD_N` | e.g. `28_10` | Number of sites (`DD`) and number of countries (`N`) |
| `filterN10` | fixed | Indicates site-level minimum count filter applied |
| `K\|O\|KO` | `K`, `O`, or `KO` | Antigen type: K locus, O locus, or both |

### Examples

| Filename | Subset | min10 | Antigen |
|----------|--------|-------|---------|
| `Neonatal_shareable_20251205_28_10_filterN10_K_site_counts_all.tsv` | Full | No | K |
| `Carba_min10_Neonatal_shareable_20251205_28_10_filterN10_O_site_counts_all.tsv` | Carba | Yes | O |
| `ESBL_Neonatal_shareable_20251205_28_10_filterN10_KO_site_counts_all.tsv` | ESBL | No | KO |

## Column descriptions

| Column | Type | Description |
|--------|------|-------------|
| `locus` | character | Serotype locus identifier (e.g. `KL1`, `O1`) |
| `Site` | character | Study site identifier |
| `Country` | character | Country of study site |
| `Region` | character | WHO region — one of: Eastern Africa, Southern Africa, Southern Asia, Western Africa |
| `raw_count` | integer | Observed count for this locus at this site (unweighted) |
| `adj_count` | numeric | Adjusted count for this locus at this site |
| `raw_sum` | integer | Total isolates at this site (raw denominator) |
| `adj_sum` | numeric | Total isolates at this site (adjusted denominator) |
| `raw_prop` | numeric | `raw_count / raw_sum` |
| `adj_prop` | numeric | `adj_count / adj_sum` |

> **Important:** The denominator for pooled regional prevalence must be computed from `distinct(Site, Region, adj_sum)` across all sites — not by summing within a locus-filtered subset. See CLAUDE.md for details.

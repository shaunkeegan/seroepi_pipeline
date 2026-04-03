# data/LOO/

Per-locus site count CSVs for the leave-one-out (`LOO`) analysis. Each file is identical in structure to the `core` site count files but with one study removed. These files are excluded from version control.

## File naming convention

```
[Subset_][min10_]Neonatal_shareable_YYYYMMDD_DD_filterN10_[K|O|KO]_site_counts_all_LOO_[StudyID].csv
```

| Component | Values | Meaning |
|-----------|--------|---------|
| `Subset_` | `Carba_`, `ESBL_`, `Fatal_`, `Full_`, *(absent)* | Clinical subset |
| `min10_` | present or absent | `min10` variant; absent = all data |
| `YYYYMMDD` | e.g. `20241002` | Data freeze date |
| `DD` | e.g. `28` | Number of sites in the full dataset before holdout |
| `filterN10` | fixed | Site-level minimum count filter applied |
| `K\|O\|KO` or `OlocusType` | `K`, `O`, `KO`, `OlocusType` | Antigen type |
| `LOO_[StudyID]` | e.g. `LOO_AIIMS`, `LOO_Kilifi` | Identifies which study was held out |

### Example study IDs held out

`AIIMS`, `AKU`, `BARNARDS`, `BabyGERMS`, `Botswana`, `CHRF`, `GBS`, `Kilifi`, `MBIRA`, `MLW`, `NeoBAC`, `NeoOBS`, `SPINZ`

## Column descriptions

Identical to `data/core/site_counts/` — see that README for the full column list. The held-out study's rows are absent from the file.

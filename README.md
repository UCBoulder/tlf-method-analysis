# TLF method paper: analysis code and data

Code and input data behind "A Method for Tryptophan-Like Fluorescence as an Alternative Fecal Indicator in
Natural Waters" (Thomas, Ross, Vlah, Knopp; ACS ES&T Water, manuscript ew-2026-014215, revision of
November 2026). Run every script from the repository root; outputs go to `output/`.

| Script | Reproduces |
|---|---|
| `scripts/recompute_ia_table.py` | Table 2 grab-and-read rows (IA and R² on log10 pairs, EPA Alternative Methods Calculator convention) from the published panel data |
| `scripts/generate_three_way_comparison.R` | The four grab-and-read panels: per-unit Colilert and MF calibrations on the paired record |
| `scripts/ia_common_set.R` | Table 2 common-set rows: all three comparisons on one sample set |
| `scripts/reference_uncertainty_propagation.R` | Section 3.1: Quanti-Tray interval propagation, expected and reference-tolerant balanced accuracy at the 82 and 126 CFU/100 mL cuts |
| `scripts/nowcast_insitu_paper.py` | In-situ arms (Section 3.2) from the public decision log, `validation.thelume.ai/api/nowcast/decisions` |
| `scripts/boulder_insitu_cv.py` | Boulder fleet regression in sample and under day-blocked cross-validation (n=155) |
| `scripts/ia_drinking_water.py` | Drinking-water rows of Table 2 against the compartment bag test (ref. 10) and the replicate-bag comparison |
| `scripts/make_fig_altcalc.py`, `make_fig_insitu.py`, `make_fig_geomean.py` | Figures 1 to 3 |

Requirements: Python 3 with numpy, pandas and matplotlib; R with tidyverse, patchwork, ggthemes, yardstick and jsonlite.

## Data (`data/`)
- `paper_clean_bucket.rds`, `paper_clean_stream.rds`: the paired grab-and-read record (chamber and field), sensor channels with Colilert and MF results.
- `quantitray_error_intervals.csv`: Quanti-Tray/2000 MPN table with 95% confidence limits.
- `bucket_combo_extended.csv`: extended January 2026 chamber record.
- `validation_panels.json`: the panel data served on `validation.thelume.ai/TLF-Method`.
- `nowcast_insitu_pairs.csv`, `nowcast_insitu_stats.json`: the in-situ pairs and statistics the paper reports, from the decision log.
- `nowcast_gm_windows.csv`: 30-day window comparison from the hourly record as committed on 9 September 2026. The hourly series endpoint has changed since, so `nowcast_insitu_paper.py` run today recomputes the pairs exactly but not these windows.
- `boulder_field_refit_20260923.json`: Boulder Creek grabs from the field-refit snapshot of 23 September 2026.
- `cbt_final_samples.csv`, `cbt_replicate_bags.csv`: drinking-water predictions and replicate bag tests from ref. 10.

## Known gap
The December 2025 MF/Colilert method-comparison sheet (90 of the 153 Colilert–MF pairs) is not included. Its
paired values are in `validation_panels.json` (`col1`), so the Colilert-vs-MF row of Table 2 recomputes; without
the sheet, `generate_three_way_comparison.R` draws that panel from the 63 pairs in the paired record only.

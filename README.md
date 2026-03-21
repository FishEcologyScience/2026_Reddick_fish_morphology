# Fish Morphology Characterization — Hamilton Harbour, Lake Ontario

**Reddick, D., Rizzuto, M., Bzonek, P., Turner, K., Croft-White, M., Midwood, J.**

This repository contains the R analysis code associated with a study characterizing gross fish morphology in Hamilton Harbour, Lake Ontario. Body width measurements (alongside fork length and mass) are analyzed across multiple species to inform the design and evaluation of structural invasive fish barriers, and to better understand how fence spacing affects fish passage.

---

## Background

Structural fish barriers are increasingly used as a management tool to limit the spread of invasive species. The effectiveness of barrier designs — particularly fence spacing — depends on the body dimensions of target species relative to non-target species. This analysis quantifies width-length and width-mass relationships across the fish community of Hamilton Harbour to provide empirical morphological data for barrier design and assessment.

---

## Getting Started

### Requirements

- R (≥ 4.1)
- RStudio (recommended)
- R packages: `tidyverse`, `patchwork`, `readxl`, `scales`

### Running the analysis

1. Open `2026_Reddick_fish_morphology.Rproj` in RStudio
2. Open and run [`02_scripts/script00-00_user_interface.R`](02_scripts/script00-00_user_interface.R)

That script sources the full pipeline in order. All outputs are written to `03_outputs/01_figures/<YYYY-mm-dd>/` (a dated subfolder created automatically on each run).

---

## Repository Structure

```
2026_Reddick_fish_morphology/
├── 01_data/
│   ├── 01_raw_files/              # Input morphometric data and species lookup
│   ├── 02_processed_files/        # Intermediate outputs (if generated)
│   └── 03_large_files_LFS/        # Large files tracked via Git LFS
├── 02_scripts/
│   ├── 01_functions/              # Helper functions sourced by main scripts
│   ├── 02_scripts_single_use/     # Supplemental one-off analyses
│   ├── 03_scripts_old/            # Archived legacy scripts
│   ├── script00-00_user_interface.R   # START HERE — runs the full pipeline
│   ├── script00-01_load_packages.R    # Package loading and global settings
│   ├── script01-01_import_format_singlefile.R  # Data import and cleaning
│   └── script01-02_morphology_plots.R          # All plots and summary tables
├── 03_outputs/
│   ├── 01_figures/                # Generated plots (dated subfolders)
│   └── 01_tables/                 # Summary statistics and model coefficients
├── 04_notes/                      # Setup guides and documentation
├── .gitignore
├── .gitattributes                 # Git LFS configuration
└── README.md
```

---

## Analysis Pipeline

```
script00-00_user_interface.R
        │
        ├── [1] script00-01_load_packages.R
        │         Loads tidyverse, patchwork; sets global options
        │
        ├── [2] script01-01_import_format_singlefile.R
        │         Reads morphometric Excel file → combined_all, df_all
        │
        └── [3] script01-02_morphology_plots.R
                  Per-species plots (6 types) + patchwork panels
                  Multi-species combined plots (5 types)
                  Summary tables: df_combined_summary, df_combined_models
```

---

## Outputs

### Per-species plots
For each species with sufficient data (n ≥ 10), six plots are produced and assembled into a single patchwork panel:

| Panel position | Plot | Description |
|---|---|---|
| (a) | Histogram | Fork length frequency distribution |
| (b) | Mass ~ Fork Length | Power-law fit: Mass = a · FL^b |
| (c) | Width ~ Fork Length | Linear fit: Width = slope · FL + intercept |
| (d) | log(Width) ~ log(Fork Length) | Log-log allometric scaling |
| (e) | Mass ~ Width | Power-law fit: Mass = a · Width^b |
| (f) | Width ~ Mass^b | Log-log power law |

### Multi-species combined plots
Five plots overlaying all species for direct comparison:
- Width ~ Fork Length
- Width ~ Mass (log x-axis)
- log(Width) ~ log(Fork Length)
- Width ~ Mass^b (log-log)
- Fork Length histograms faceted by species

### Summary tables
- **`df_combined_summary`** — per-species n, mean, min, max for fork length, width, and mass
- **`df_combined_models`** — power-law coefficients (a, b, R²) from Mass ~ Width fits
- **`df_species_counts`** — raw vs. NA/zero-filtered sample sizes per species

---

## Data

Input data are stored in `01_data/01_raw_files/`. The primary file is a consolidated Excel spreadsheet containing morphometric measurements for all species. Required columns:

| Column | Units | Description |
|---|---|---|
| `Species` | — | Species code |
| `ForkLength_mm` | mm | Fork length |
| `Width_mm` | mm | Maximum body width |
| `Mass_g` | g | Wet mass |

---

## Supplemental Analyses

Scripts in `02_scripts/02_scripts_single_use/` contain additional analyses run independently of the main pipeline:

| Script | Purpose |
|---|---|
| `script01_02_elecfish_length_histograms_2012on.R` | Length-frequency distributions from electrofishing surveys (2012–present) |
| `script01_03_gap_analysis_width.R` | Identifies sampling gaps in width measurements relative to the 50 mm threshold |
| `script01_04_rlf_largemouth_bass_bonar_bins.R` | Relative length frequency analysis for Largemouth Bass (Bonar 2002 methodology) |
| `script01_06_plot_height_length_width.R` | Height, length, and width morphology plots |

---

## Citation

A formal citation will be available upon publication. This repository will be archived on Zenodo. In the meantime, please contact the corresponding author for attribution guidance.

---

## Additional Documentation

- [`04_notes/naming_conventions.md`](04_notes/naming_conventions.md) — file, folder, and R object naming standards
- [`04_notes/git_lfs_setup.md`](04_notes/git_lfs_setup.md) — Git LFS installation and troubleshooting

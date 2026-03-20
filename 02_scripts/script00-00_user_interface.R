## --------------------------------------------------------------#
## Script name: script00-00_user_interface.R
##
## Purpose:
##    Master control script — sources all pipeline scripts in
##    dependency order. Run this file to execute the full workflow.
##
##    Script numbering convention:
##      scriptXX-YY  (XX = class, YY = script within class)
##      Higher class numbers depend on lower class numbers.
##      Same-class scripts with letters (a, b) have no dependency.
##
## Author: Paul Bzonek [Claude]
## Date Created: 2026-03-20
##
## --------------------------------------------------------------#
## Modification Notes:
##
## --------------------------------------------------------------#


##### Pipeline ################################################----

# [1] Packages, helper functions, global settings
source("02_scripts/script00-01_load_packages.R")
# Produces: param_seed

# [2] Import and clean morphometric data from single Excel file
source("02_scripts/script01-01_import_format_singlefile.R")
# Produces: combined_all (named list, one data frame per species)
#           df_all       (single combined tibble, all species)

# [3] Per-species and multi-species morphology plots + summary tables
source("02_scripts/script01-02_morphology_plots.R")
# Produces: plots[[species]][[slot]]     (per-species plots + patchwork panels)
#           plots[["combined"]][[slot]]  (multi-species combined plots)
#           df_combined_summary          (descriptive stats per species)
#           df_combined_models           (power-law coefficients per species)
#           df_species_counts            (raw vs filtered sample sizes)

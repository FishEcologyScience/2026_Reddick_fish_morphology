## --------------------------------------------------------------#
## Script name: script00-01_load_packages.R
##
## Purpose:
##    Load packages, source helper functions, and set global options.
##    Only packages actively used by this project are loaded here.
##    Rarely used or context-specific packages should be loaded with
##    package::function() notation directly in the scripts that need them.
##
## Author: Paul Bzonek [Claude]
## Date Created: 2026-03-20
##
## --------------------------------------------------------------#
## Modification Notes:
##   2026-03-20 — Removed ggmap, sf, FESLtelemetry (not used in morphology
##                pipeline). Commented out devtools install (run manually
##                during updates only). Removed theme_set() — theme is
##                managed in script01-02_morphology_plots.R.
## --------------------------------------------------------------#


##### CRAN packages ###########################################----
library(tidyverse)   # core data wrangling and ggplot2
library(patchwork)   # combine multiple ggplots into panels

# Spatial packages — not needed for morphology pipeline; load if extending
# library(ggmap)     # google basemaps
# library(sf)        # shapefiles, projections, distance calculations


##### Lab packages ############################################----
# Install manually during updates only — do NOT leave uncommented during
# routine runs (slow network call on every source).
# devtools::install_github("FishEcologyScience/FESLtelemetry")
# library(FESLtelemetry)


##### Helper functions ########################################----
source("02_scripts/01_functions/function01-01_helper_functions.R")


##### Global settings #########################################----
options(scipen = 999)    # suppress scientific notation
param_seed <- 1987       # random seed for any randomized steps

# Note: ggplot theme is set in script01-02_morphology_plots.R to keep
# theme configuration alongside the plotting code that uses it.

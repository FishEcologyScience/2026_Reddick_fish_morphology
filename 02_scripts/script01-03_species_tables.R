## --------------------------------------------------------------#
## Script name: script01-03_species_tables.R
##
## Purpose:
##   - Use pre-cleaned data (from script01-01_import_format_singlefile.R)
##     to generate species-level summary and results tables.
##   - Outputs include:
##       (1) Descriptive morphology statistics by species
##           (sample size, fork length range, mean fork length,
##            mean width, and SD of width)
##       (2) Species-level morphology model results
##           (coefficients, R², and equation forms for fitted relationships)
##
## Repository layout:
##   01_data/01_raw_files/          # inputs (species lookup)
##   03_outputs/01_tables/          # tables (Excel outputs)
##
## Author: Marcus Rizzuto
## Date Created: 2026-03-27
##
## --------------------------------------------------------------#
##
## TABLE INVENTORY (what this script produces)
##
##   1) species_morphology_summary.xlsx
##        - Per-species descriptive statistics:
##          n, fork length range, mean fork length,
##          mean width, SD of width
##
##   2) morphology_model_results.xlsx
##        - Per-species fitted morphology model results:
##          model type, coefficients, R², and equation form
##
## --------------------------------------------------------------#


##### Setup ##########################################----
#-------------------------------------------------------------#
# This script assumes a completed import/cleaning step.
# The object `combined_all` should already exist in the
# global environment as a named list of cleaned data frames,
# one per species.

if (!exists("combined_all")) {
 stop("ERROR: `combined_all` not found. Run script01-01_import_format_singlefile.R first.")
}

# Tables are written to a fixed output directory so results
# remain consistent across runs and scripts.
path_tables_dir <- file.path("03_outputs", "01_tables")
if (!dir.exists(path_tables_dir)) {
 dir.create(path_tables_dir, recursive = TRUE, showWarnings = FALSE)
}


##### Species name lookup #####################################----
#-------------------------------------------------------------#
# Load the lookup table used to convert internal species
# identifiers (short codes) into human-readable names
# for final reporting outputs.

species_lookup_path <- file.path(
 "01_data", "01_raw_files", "species_lookup.xlsx"
)

species_lookup <- readxl::read_excel(species_lookup_path) |>
 dplyr::rename(
  short_form  = `Short Form`,
  common_name = `Common Name`
 ) |>
 dplyr::mutate(
  short_form  = as.character(short_form),
  common_name = as.character(common_name)
 )


##### Species-level summary table #############################----
#-------------------------------------------------------------#
# Builds a compact per-species descriptive statistics table.
#
# This section:
#   - Iterates once over each species data frame
#   - Computes basic size statistics directly from cleaned data
#   - Does NOT apply minimum-n filters (summary only)
#
# This table is intended for reporting, QA/QC, and overview use.

df_species_summary <- purrr::imap_dfr(
 combined_all,
 function(df_clean, param_species) {
  
  # Guard against blank or placeholder species
  sp_key <- tolower(trimws(param_species))
  if (is.na(param_species) || sp_key %in% c("", "unknown", "unknown_species")) {
   return(NULL)
  }
  
  tibble::tibble(
   species_short = param_species,
   
   # Count based on fork length availability (primary size metric)
   n             = sum(!is.na(df_clean$ForkLength_mm)),
   
   # Range calculations are guarded to avoid Inf when all values are NA
   FL_min_mm     = if (all(is.na(df_clean$ForkLength_mm))) NA_real_
   else suppressWarnings(min(df_clean$ForkLength_mm, na.rm = TRUE)),
   FL_max_mm     = if (all(is.na(df_clean$ForkLength_mm))) NA_real_
   else suppressWarnings(max(df_clean$ForkLength_mm, na.rm = TRUE)),
   
   FL_mean_mm    = mean(df_clean$ForkLength_mm, na.rm = TRUE),
   Width_mean_mm = mean(df_clean$Width_mm,      na.rm = TRUE),
   
   # Width SD retained as a simple indication of spread/variability
   Width_sd_mm   = stats::sd(df_clean$Width_mm, na.rm = TRUE)
  )
 }
)


##### Species-level morphology model results ##################----
#-------------------------------------------------------------#
# Model coefficients are sourced directly from df_combined_models,
# which is built during the plotting step (script01-02_morphology_plots.R).
# This ensures the table coefficients are identical to what appears on the plots.

if (!exists("df_combined_models")) {
 stop("ERROR: `df_combined_models` not found. Run script01-02_morphology_plots.R first.")
}

df_model_results <- df_combined_models


##### Resolve species codes ###################################
#-------------------------------------------------------------#
# Replace short-form species identifiers with common names
# where available for final export tables.

df_model_results <- df_model_results |>
 dplyr::left_join(
  species_lookup |>
   dplyr::select(short_form, common_name),
  by = c("Species name" = "short_form")
 ) |>
 dplyr::mutate(
  `Species name` = dplyr::if_else(
   is.na(common_name) | common_name == "",
   `Species name`,
   common_name
  )
 ) |>
 dplyr::select(-common_name)


##### Formatting and name resolution ##########################----
#-------------------------------------------------------------#

df_species_summary <- df_species_summary |>
 dplyr::left_join(
  species_lookup,
  by = c("species_short" = "short_form")
 ) |>
 dplyr::mutate(
  common_name = dplyr::if_else(
   is.na(common_name) | common_name == "",
   species_short,
   common_name
  ),
  `Fork length range (mm)` = dplyr::if_else(
   is.na(FL_min_mm) | is.na(FL_max_mm),
   NA_character_,
   paste0(round(FL_min_mm), "–", round(FL_max_mm))
  )
 ) |>
 dplyr::select(
  `Species`                = common_name,
  `n`,
  `Fork length range (mm)`,
  `Mean fork length (mm)`  = FL_mean_mm,
  `Mean width (mm)`        = Width_mean_mm,
  `SD width (mm)`          = Width_sd_mm
 ) |>
 dplyr::mutate(
  `Mean fork length (mm)` = round(`Mean fork length (mm)`, 1),
  `Mean width (mm)`       = round(`Mean width (mm)`, 1),
  `SD width (mm)`         = round(`SD width (mm)`, 1)
 ) |>
 dplyr::arrange(dplyr::desc(n))


##### Export summary table ####################################----
#-------------------------------------------------------------#

output_path <- file.path(
 path_tables_dir,
 "species_morphology_summary.xlsx"
)

writexl::write_xlsx(df_species_summary, output_path)

cat("Excel table written to:\n", output_path, "\n")


##### Export model results table ##############################----
#-------------------------------------------------------------#

writexl::write_xlsx(
 df_model_results,
 file.path(path_tables_dir, "morphology_model_results.xlsx")
)

cat("Model results table written to:\n",
    file.path(path_tables_dir, "morphology_model_results.xlsx"), "\n")
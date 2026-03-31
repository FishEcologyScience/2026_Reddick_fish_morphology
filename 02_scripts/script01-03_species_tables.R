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
# Fits morphology relationships at the species level and compiles
# results into a single long-format table.
#
# Key design choices:
#   - All models use the same minimum sample size threshold
#   - Failed fits are silently skipped (no placeholder rows)
#   - R² is always calculated on the scale of the response variable
#
# This table complements the plotting script by providing
# exact numeric results in exportable form.

MIN_N_PER_SPECIES <- 10

df_model_results <- tibble()

for (param_species in names(combined_all)) {
 
 df_clean <- combined_all[[param_species]]
 
 # Skip invalid species identifiers early
 sp_key <- tolower(trimws(param_species))
 if (is.na(param_species) || sp_key %in% c("", "unknown", "unknown_species")) next
 
 
 ### Width ~ Fork Length ######################################
 # Linear fit is used here as a shape diagnostic.
 # Filtering is limited to the variables involved in the model
 # to maximize usable sample size.
 
 df_w_fl <- df_clean |>
  dplyr::filter(!is.na(Width_mm), !is.na(ForkLength_mm))
 
 if (nrow(df_w_fl) >= MIN_N_PER_SPECIES) {
  fit   <- lm(Width_mm ~ ForkLength_mm, data = df_w_fl)
  coefs <- coef(fit)
  
  df_model_results <- dplyr::bind_rows(
   df_model_results,
   tibble(
    `Species name` = param_species,
    Relationship   = "Width ~ Fork Length",
    n              = nrow(df_w_fl),
    intercept_a    = unname(coefs[1]),
    slope_b        = unname(coefs[2]),
    `R^2`          = summary(fit)$r.squared,
    Equation       = paste0(
     "Width = ",
     sprintf("%.4f", coefs[2]),
     " · FL + ",
     sprintf("%.4f", coefs[1])
    )
   )
  )
 }
 
 
 ### Mass ~ Width #############################################
 # Uses nonlinear least squares on the original measurement scale.
 # A preliminary log–log fit is used *only* to obtain stable
 # starting values for the nonlinear solver.
 
 df_m_w <- df_clean |>
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(Mass_g),   Mass_g   > 0)
 
 if (nrow(df_m_w) >= MIN_N_PER_SPECIES) {
  
  lm_start <- lm(log(Mass_g) ~ log(Width_mm), data = df_m_w)
  a0 <- exp(coef(lm_start)[1])
  b0 <- coef(lm_start)[2]
  
  nls_fit <- try(
   nls(Mass_g ~ a * Width_mm^b,
       data = df_m_w,
       start = list(a = a0, b = b0),
       control = nls.control(maxiter = 200, warnOnly = TRUE)),
   silent = TRUE
  )
  
  if (!inherits(nls_fit, "try-error")) {
   a <- coef(nls_fit)["a"]
   b <- coef(nls_fit)["b"]
   
   # Pseudo-R² on the response scale
   y    <- df_m_w$Mass_g
   yhat <- fitted(nls_fit)
   r2   <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
   
   df_model_results <- dplyr::bind_rows(
    df_model_results,
    tibble(
     `Species name` = param_species,
     Relationship   = "Mass ~ Width",
     n              = nrow(df_m_w),
     intercept_a    = a,
     slope_b        = b,
     `R^2`          = r2,
     Equation       = paste0(
      "Mass = ",
      sprintf("%.2e", a),
      " · Width^",
      sprintf("%.4f", b)
     )
    )
   )
  }
 }
 
 
 ### Mass ~ Fork Length #######################################
 # Same fitting strategy as Mass ~ Width, but using length
 # as the predictor. Kept separate to allow direct comparison
 # of model behavior between morphometrics.
 
 df_m_fl <- df_clean |>
  dplyr::filter(!is.na(ForkLength_mm), ForkLength_mm > 0,
                !is.na(Mass_g),       Mass_g       > 0)
 
 if (nrow(df_m_fl) >= MIN_N_PER_SPECIES) {
  
  lm_start <- lm(log(Mass_g) ~ log(ForkLength_mm), data = df_m_fl)
  a0 <- exp(coef(lm_start)[1])
  b0 <- coef(lm_start)[2]
  
  nls_fit <- try(
   nls(Mass_g ~ a * ForkLength_mm^b,
       data = df_m_fl,
       start = list(a = a0, b = b0),
       control = nls.control(maxiter = 200, warnOnly = TRUE)),
   silent = TRUE
  )
  
  if (!inherits(nls_fit, "try-error")) {
   a <- coef(nls_fit)["a"]
   b <- coef(nls_fit)["b"]
   
   y    <- df_m_fl$Mass_g
   yhat <- fitted(nls_fit)
   r2   <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
   
   df_model_results <- dplyr::bind_rows(
    df_model_results,
    tibble(
     `Species name` = param_species,
     Relationship   = "Mass ~ Fork Length",
     n              = nrow(df_m_fl),
     intercept_a    = a,
     slope_b        = b,
     `R^2`          = r2,
     Equation       = paste0(
      "Mass = ",
      sprintf("%.2e", a),
      " · FL^",
      sprintf("%.4f", b)
     )
    )
   )
  }
 }
 
 
 ### log(Width) ~ log(Fork Length) ############################
 # Performed to mirror the log–log diagnostic plots and
 # provide a directly comparable scaling exponent in tabular form.
 
 df_ll <- df_clean |>
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(ForkLength_mm), ForkLength_mm > 0)
 
 if (nrow(df_ll) >= MIN_N_PER_SPECIES) {
  fit   <- lm(log(Width_mm) ~ log(ForkLength_mm), data = df_ll)
  coefs <- coef(fit)
  
  df_model_results <- dplyr::bind_rows(
   df_model_results,
   tibble(
    `Species name` = param_species,
    Relationship   = "log(Width) ~ log(Fork Length)",
    n              = nrow(df_ll),
    intercept_a    = exp(coefs[1]),
    slope_b        = coefs[2],
    `R^2`          = summary(fit)$r.squared,
    Equation       = paste0(
     "Width = ",
     sprintf("%.4f", exp(coefs[1])),
     " · FL^",
     sprintf("%.4f", coefs[2])
    )
   )
  )
 }
 
 
 ### Width ~ Mass^b ###########################################
 # Mirrors the final power-law plot where mass is the driver
 # and width is the response. Fit is performed in log space
 # to prioritize interpretation over prediction.
 
 df_pw <- df_clean |>
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(Mass_g),   Mass_g   > 0)
 
 if (nrow(df_pw) >= MIN_N_PER_SPECIES) {
  fit   <- lm(log(Width_mm) ~ log(Mass_g), data = df_pw)
  coefs <- coef(fit)
  
  df_model_results <- dplyr::bind_rows(
   df_model_results,
   tibble(
    `Species name` = param_species,
    Relationship   = "Width ~ Mass^b",
    n              = nrow(df_pw),
    intercept_a    = exp(coefs[1]),
    slope_b        = coefs[2],
    `R^2`          = summary(fit)$r.squared,
    Equation       = paste0(
     "Width = ",
     sprintf("%.4f", exp(coefs[1])),
     " · Mass^",
     sprintf("%.4f", coefs[2])
    )
   )
  )
 }
}


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
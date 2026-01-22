## --------------------------------------------------------------#
## Script name: script_morphology_plots.R
##
## Purpose of script:
##   - Run morphology plotting workflow for MULTIPLE species (script only)
##   - Avoid overwriting by looping over species and using species-specific paths
##   - Save per-species outputs + combined summary tables
##
## Repository layout (per template):
##   01_data/                       # inputs
##   01_data/processed/<species>/   # cleaned RDS per species
##   03_outputs/01_figures/         # figures
##   03_outputs/01_tables/          # tables (created if missing)
##
## Author: Marcus Rizzuto
## Date Created: 1/22/2026
##
## --------------------------------------------------------------#


## Packages
library(tidyverse)
library(readxl)
library(lubridate)

options(scipen = 999)

## User inputs
param_species_vec    <- c("Rudd", "Goldfish", "Carp")
param_recursive_read <- FALSE

## Paths
path_data_root  <- "01_data"
path_raw_root   <- file.path(path_data_root, "01_raw_files", "Species")
path_figs_dir   <- file.path("03_outputs", "01_figures")
path_tables_dir <- file.path("03_outputs", "01_tables")

## In-memory combined outputs
df_combined_summary <- tibble()
df_combined_models  <- tibble()

## Loop
for (param_species in param_species_vec) {
 
 cat("\n--- Processing:", param_species, "---\n")
 
 path_raw_dir       <- file.path(path_raw_root, param_species)
 path_processed_dir <- file.path(path_data_root, "02_processed_files", param_species)
 
 if (!dir.exists(path_raw_dir)) { cat("WARNING: missing raw dir:", path_raw_dir, "\n"); next }
 
 ## Read & combine
 temp_raw_files <- list.files(path_raw_dir, "\\.xlsx$", ignore.case = TRUE,
                              full.names = TRUE, recursive = param_recursive_read)
 if (length(temp_raw_files) == 0) { cat("WARNING: no .xlsx in:", path_raw_dir, "\n"); next }
 
 cat("Found", length(temp_raw_files), "files; combining…\n")
 temp_df_list    <- lapply(temp_raw_files, readxl::read_excel)
 df_raw_combined <- dplyr::bind_rows(temp_df_list)
 
 ## Clean
 df_clean <- df_raw_combined %>%
  mutate(
   ForkLength_mm = suppressWarnings(as.numeric(ForkLength_mm)),
   Width_mm      = suppressWarnings(as.numeric(Width_mm)),
   Mass_g        = suppressWarnings(as.numeric(Mass_g))
  )
 
 temp_required_cols <- c("ForkLength_mm", "Width_mm", "Mass_g")
 temp_missing <- setdiff(temp_required_cols, names(df_clean))
 if (length(temp_missing)) { cat("WARNING: missing cols:", paste(temp_missing, collapse = ", "), "\n"); next }
 
 ## Summary (kept in memory)
 df_summary <- df_clean %>%
  summarise(
   species       = param_species,
   n_rows        = n(),
   n_FL          = sum(!is.na(ForkLength_mm)),
   n_width       = sum(!is.na(Width_mm)),
   n_weight      = sum(!is.na(Mass_g)),
   FL_mean_mm    = mean(ForkLength_mm, na.rm = TRUE),
   width_mean_mm = mean(Width_mm,      na.rm = TRUE),
   weight_mean_g = mean(Mass_g,        na.rm = TRUE)
  )
 df_combined_summary <- dplyr::bind_rows(df_combined_summary, df_summary)
 
 make_caption <- function(df, text) paste0(param_species, " (n = ", nrow(df), "): ", text)
 
 ## Plot 1 — Width vs Fork Length (linear fit + equation)
 df_scatter_fl <- df_clean %>% filter(!is.na(Width_mm), !is.na(ForkLength_mm))
 
 if (nrow(df_scatter_fl) >= 2) {
  temp_lm_fl  <- lm(Width_mm ~ ForkLength_mm, data = df_scatter_fl)
  temp_coef   <- coef(temp_lm_fl)
  temp_r2_fl  <- summary(temp_lm_fl)$r.squared
  temp_slope  <- unname(temp_coef[["ForkLength_mm"]])
  temp_int    <- unname(temp_coef[["(Intercept)"]])
  
  temp_eq_fl <- paste0(
   "y = ", formatC(temp_slope, format = "f", digits = 2),
   "x", ifelse(temp_int >= 0, " + ", " - "),
   formatC(abs(temp_int), format = "f", digits = 2),
   "\nR^2 = ", formatC(temp_r2_fl, format = "f", digits = 4)
  )
  
  temp_xpos_fl <- quantile(df_scatter_fl$ForkLength_mm, 0.05, na.rm = TRUE)
  temp_ypos_fl <- quantile(df_scatter_fl$Width_mm,       0.95, na.rm = TRUE)
  
  p_scatter_fl <- ggplot(df_scatter_fl, aes(ForkLength_mm, Width_mm)) +
   geom_point(color = "#2c7fb8", alpha = 0.6, size = 2) +
   geom_smooth(method = "lm", se = FALSE, color = "#1f78b4", linewidth = 0.9) +
   annotate("text", x = temp_xpos_fl, y = temp_ypos_fl, label = temp_eq_fl,
            hjust = 0, vjust = 1, size = 3.5) +
   labs(
    title   = paste0(param_species, " - Width by Fork Length"),
    x       = "Fork Length (mm)",
    y       = "Width (mm)",
    caption = make_caption(df_scatter_fl, "Width vs fork length (linear fit).")
   ) +
   theme_minimal()
 } else {
  p_scatter_fl <- ggplot(df_scatter_fl, aes(ForkLength_mm, Width_mm)) +
   geom_point(color = "#2c7fb8", alpha = 0.6, size = 2) +
   labs(
    title   = paste0(param_species, " - Width by Fork Length"),
    x       = "Fork Length (mm)",
    y       = "Width (mm)",
    caption = make_caption(df_scatter_fl, "Width vs fork length.")
   ) +
   theme_minimal()
 }
 p_scatter_fl
 # path_file_fl <- file.path(path_figs_dir, paste0(param_species, "_scatter_width_by_forklength.png"))
 # ggsave(filename = path_file_fl, plot = p_scatter_fl, width = 7, height = 5, dpi = 300)
 
 ## Plot 2 — Width vs Mass (log fit + equation)
 df_scatter_mass <- df_clean %>% filter(!is.na(Width_mm), !is.na(Mass_g), Mass_g > 0)
 
 temp_model_row <- tibble(species = param_species, alpha = NA_real_, beta = NA_real_, r2 = NA_real_)
 
 if (nrow(df_scatter_mass) >= 3) {
  temp_log_fit <- lm(Width_mm ~ log(Mass_g), data = df_scatter_mass)
  
  temp_coefs <- coef(temp_log_fit)
  temp_alpha <- unname(temp_coefs["log(Mass_g)"])
  temp_beta  <- unname(temp_coefs["(Intercept)"])
  temp_r2_m  <- summary(temp_log_fit)$r.squared
  
  temp_model_row <- tibble(species = param_species, alpha = temp_alpha, beta = temp_beta, r2 = temp_r2_m)
  
  temp_eq_m <- paste0(
   "y = ", formatC(temp_alpha, format = "f", digits = 2),
   " ln(x)", ifelse(temp_beta >= 0, " + ", " - "),
   formatC(abs(temp_beta), format = "f", digits = 2),
   "\nR^2 = ", formatC(temp_r2_m, format = "f", digits = 4)
  )
  
  temp_xrng <- range(df_scatter_mass$Mass_g, na.rm = TRUE)
  temp_xseq <- seq(temp_xrng[1], temp_xrng[2], length.out = 200)
  df_pred   <- tibble(Mass_g = temp_xseq) %>%
   mutate(Width_pred = predict(temp_log_fit, newdata = tibble(Mass_g = Mass_g)))
  
  temp_xpos_m <- quantile(df_scatter_mass$Mass_g, 0.05, na.rm = TRUE)
  temp_ypos_m <- quantile(df_scatter_mass$Width_mm, 0.95, na.rm = TRUE)
  
  p_scatter_mass <- ggplot(df_scatter_mass, aes(Mass_g, Width_mm)) +
   geom_point(color = "#f16913", alpha = 0.6, size = 2) +
   geom_line(data = df_pred, aes(Mass_g, Width_pred), color = "#cc4c02", linewidth = 1) +
   annotate("text", x = temp_xpos_m, y = temp_ypos_m, label = temp_eq_m,
            hjust = 0, vjust = 1, size = 3.5) +
   labs(
    title   = paste0(param_species, " - Width by Mass"),
    x       = "Mass (g)",
    y       = "Width (mm)",
    caption = make_caption(df_scatter_mass, "Width vs mass (logarithmic fit).")
   ) +
   theme_minimal()
 } else {
  p_scatter_mass <- ggplot(df_scatter_mass, aes(Mass_g, Width_mm)) +
   geom_point(color = "#f16913", alpha = 0.6, size = 2) +
   labs(
    title   = paste0(param_species, " - Width by Mass"),
    x       = "Mass (g)",
    y       = "Width (mm)",
    caption = make_caption(df_scatter_mass, "Width vs mass.")
   ) +
   theme_minimal()
 }
 p_scatter_mass
 # path_file_mass <- file.path(path_figs_dir, paste0(param_species, "_scatter_width_by_mass.png"))
 # ggsave(filename = path_file_mass, plot = p_scatter_mass, width = 7, height = 5, dpi = 300)
 
 ## (Optional) per-species file saves — commented on purpose
 # path_clean_rds <- file.path(path_processed_dir, paste0(param_species, "_clean.rds"))
 # saveRDS(df_clean, path_clean_rds)
 # path_species_summary_csv <- file.path(path_tables_dir, paste0(param_species, "_summary.csv"))
 # readr::write_csv(df_summary, path_species_summary_csv)
 
 ## Keep combined model row
 df_combined_models <- dplyr::bind_rows(df_combined_models, temp_model_row)
 
 cat("Finished:", param_species, "\n")
}

## (Optional) combined table saves — commented on purpose
# path_combined_summary_csv <- file.path(path_tables_dir, "_combined_summary_by_species.csv")
# readr::write_csv(df_combined_summary, path_combined_summary_csv)
# path_combined_models_csv  <- file.path(path_tables_dir, "_combined_log_model_coefficients.csv")
# readr::write_csv(df_combined_models,  path_combined_models_csv)

## Cleanup temporary objects
temp_objects <- ls(pattern = "^temp_")
if (length(temp_objects)) { rm(list = temp_objects); cat("Removed temporary objects.\n") }

cat("\nAll done.\n")

  
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


## Load Packages
## ----------------------------#
library('tidyverse')
library('readxl')
library('lubridate')


## Adjust Settings (optional)
## ----------------------------#
# theme_set(theme_classic())    # uncomment if you want a global ggplot theme
options(scipen = 999)
param_seed <- 1987
plots <- list()


## User Inputs
## ----------------------------#
# ⚠️ Folder names under '01_data/01_raw_files/Species/' must match exactly.
species_vec <- c("Rudd", "Goldfish", "Carp")   # <- edit this list as needed
recursive_read <- FALSE                        # TRUE if .xlsx are in subfolders per species


## Canonical Paths
## ----------------------------#
data_root  <- "01_data"
raw_root   <- file.path(data_root, "01_raw_files", "Species")
figs_dir   <- file.path("03_outputs", "01_figures")
tables_dir <- file.path("03_outputs", "01_tables")

dir.create(figs_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)


## Holders for combined outputs (across species)
## ----------------------------#
combined_summary_list <- list()
combined_model_list   <- list()


## Loop through species (no custom functions)
## ----------------------------#
for (sp in species_vec) {
 
 message("\n--- Processing species: ", sp, " ---")
 
 # Paths specific to this species
 raw_dir       <- file.path(raw_root, sp)
 processed_dir <- file.path(data_root, "02_processed_files", sp)
 
 # Ensure dirs exist
 if (!dir.exists(raw_dir)) {
  warning("Raw data directory does not exist for ", sp, ": ", raw_dir, ". Skipping this species.")
  next
 }
 dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
 
 # 1) Read & Combine Excel files
 raw_files <- list.files(
  raw_dir,
  pattern = "\\.xlsx$",
  ignore.case = TRUE,
  full.names = TRUE,
  recursive = recursive_read
 )
 
 if (length(raw_files) == 0) {
  warning("No Excel (.xlsx) files found for ", sp, " in: ", raw_dir, ". Skipping this species.")
  next
 }
 
 message("Found ", length(raw_files), " Excel files for ", sp, ". Combining...")
 df_list    <- lapply(raw_files, readxl::read_excel)
 data_field <- bind_rows(df_list)
 
 # 2) Minimal cleaning
 clean_df <- data_field %>%
  mutate(
   ForkLength_mm = suppressWarnings(as.numeric(ForkLength_mm)),
   Width_mm      = suppressWarnings(as.numeric(Width_mm)),
   Mass_g        = suppressWarnings(as.numeric(Mass_g))
  )
 
 # Optional: column check
 required_cols <- c("ForkLength_mm", "Width_mm", "Mass_g")
 missing <- setdiff(required_cols, names(clean_df))
 if (length(missing)) {
  warning("Missing required columns for ", sp, ": ", paste(missing, collapse = ", "),
          ". Skipping this species.")
  next
 }
 
 # 3) Save cleaned data
 clean_rds_path <- file.path(processed_dir, paste0(sp, "_clean.rds"))
 saveRDS(clean_df, clean_rds_path)
 
 # 4) Summary table (per species)
 summary_tbl <- clean_df %>%
  summarise(
   species       = sp,
   n_rows        = n(),
   n_FL          = sum(!is.na(ForkLength_mm)),
   n_width       = sum(!is.na(Width_mm)),
   n_weight      = sum(!is.na(Mass_g)),
   FL_mean_mm    = mean(ForkLength_mm, na.rm = TRUE),
   width_mean_mm = mean(Width_mm,      na.rm = TRUE),
   weight_mean_g = mean(Mass_g,        na.rm = TRUE)
  )
 
 species_summary_csv <- file.path(tables_dir, paste0(sp, "_summary.csv"))
 readr::write_csv(summary_tbl, species_summary_csv)
 combined_summary_list[[sp]] <- summary_tbl
 
 # 5) Plot 1: Width vs Fork Length
 caption_n <- function(df, text) paste0(sp, " (n = ", nrow(df), "): ", text)
 
 df_scatter_fl <- clean_df %>%
  filter(!is.na(Width_mm), !is.na(ForkLength_mm))
 
 p_scatter_fl <- ggplot(df_scatter_fl, aes(x = ForkLength_mm, y = Width_mm)) +
  geom_point(color = "#2c7fb8", alpha = 0.6, size = 2) +
  labs(
   title   = paste0(sp, " - Width by Fork Length"),
   x       = "Fork Length (mm)",
   y       = "Width (mm)",
   caption = caption_n(df_scatter_fl, "Width plotted against fork length.")
  ) +
  theme_minimal()
 
 file_fl <- file.path(figs_dir, paste0(sp, "_scatter_width_by_forklength.png"))
 ggsave(filename = file_fl, plot = p_scatter_fl, width = 7, height = 5, dpi = 300)
 
 # 6) Plot 2: Width vs Mass (logarithmic fit + equation)
 df_scatter_mass <- clean_df %>%
  filter(!is.na(Width_mm), !is.na(Mass_g)) %>%
  filter(Mass_g > 0)  # log requires positive mass
 
 # Default, in case we can't fit a model
 model_row <- tibble(
  species = sp,
  alpha   = NA_real_,
  beta    = NA_real_,
  r2      = NA_real_
 )
 
 if (nrow(df_scatter_mass) >= 3) {
  
  # Fit: Width = alpha * ln(Mass) + beta
  log_fit <- lm(Width_mm ~ log(Mass_g), data = df_scatter_mass)
  
  x_rng <- range(df_scatter_mass$Mass_g, na.rm = TRUE)
  x_seq <- seq(x_rng[1], x_rng[2], length.out = 200)
  
  pred_df <- tibble(Mass_g = x_seq) %>%
   mutate(Width_pred = predict(log_fit, newdata = tibble(Mass_g = Mass_g)))
  
  # Coefficients and R^2
  coefs <- coef(log_fit)
  alpha <- unname(coefs["log(Mass_g)"])
  beta  <- unname(coefs["(Intercept)"])
  r2    <- summary(log_fit)$r.squared
  
  model_row <- tibble(
   species = sp,
   alpha   = alpha,
   beta    = beta,
   r2      = r2
  )
  
  # Equation label: y = 11.97 ln(x) - 31.54   R^2 = 0.9819
  alpha_lab <- formatC(alpha, format = "f", digits = 2)
  beta_lab  <- formatC(abs(beta), format = "f", digits = 2)
  sign_beta <- ifelse(beta >= 0, " + ", " - ")
  r2_lab    <- formatC(r2, format = "f", digits = 4)
  eq_label  <- paste0("y = ", alpha_lab, " ln(x)", sign_beta, beta_lab, "\nR^2 = ", r2_lab)
  
  x_pos <- quantile(df_scatter_mass$Mass_g, 0.05, na.rm = TRUE)
  y_pos <- quantile(df_scatter_mass$Width_mm, 0.95, na.rm = TRUE)
  
  p_scatter_mass <- ggplot(df_scatter_mass, aes(x = Mass_g, y = Width_mm)) +
   geom_point(color = "#f16913", alpha = 0.6, size = 2) +
   geom_line(data = pred_df, aes(x = Mass_g, y = Width_pred),
             color = "#cc4c02", linewidth = 1) +
   annotate("text", x = x_pos, y = y_pos, label = eq_label,
            hjust = 0, vjust = 1, size = 3.5) +
   labs(
    title   = paste0(sp, " - Width by Mass"),
    x       = "Mass (g)",
    y       = "Width (mm)",
    caption = caption_n(df_scatter_mass, "Width vs mass with logarithmic fit.")
   ) +
   theme_minimal()
  
 } else {
  warning("Not enough non-missing Mass/Width pairs to fit a logarithmic model for ", sp)
  p_scatter_mass <- ggplot(df_scatter_mass, aes(x = Mass_g, y = Width_mm)) +
   geom_point(color = "#f16913", alpha = 0.6, size = 2) +
   labs(
    title   = paste0(sp, " - Width by Mass"),
    x       = "Mass (g)",
    y       = "Width (mm)",
    caption = caption_n(df_scatter_mass, "Width plotted against mass.")
   ) +
   theme_minimal()
 }
 
 file_mass <- file.path(figs_dir, paste0(sp, "_scatter_width_by_mass.png"))
 ggsave(filename = file_mass, plot = p_scatter_mass, width = 7, height = 5, dpi = 300)
 
 # Keep combined outputs
 combined_model_list[[sp]] <- model_row
 
 message("Finished: ", sp)
}


## Write combined tables (across species)
## ----------------------------#
if (length(combined_summary_list) > 0) {
 combined_summary <- bind_rows(combined_summary_list)
 combined_summary_csv <- file.path(tables_dir, "_combined_summary_by_species.csv")
 readr::write_csv(combined_summary, combined_summary_csv)
 message("Combined summary written: ", combined_summary_csv)
} else {
 warning("No per-species summaries were produced. Check species and data.")
}

if (length(combined_model_list) > 0) {
 combined_models <- bind_rows(combined_model_list)
 combined_models_csv <- file.path(tables_dir, "_combined_log_model_coefficients.csv")
 readr::write_csv(combined_models, combined_models_csv)
 message("Combined model coefficients written: ", combined_models_csv)
} else {
 warning("No per-species models were produced (model fit likely failed due to insufficient data).")
}

message("\nAll done.")


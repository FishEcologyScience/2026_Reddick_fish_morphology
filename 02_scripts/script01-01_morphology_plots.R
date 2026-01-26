## --------------------------------------------------------------#
## Script name: script01-01_morphology_plots.R
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
## Modification notes:
##   2026-01-26: Added formatting and explanatory comments 
##
## --------------------------------------------------------------#


#####Setup#####################################################----
#-------------------------------------------------------------#

### User Parameters
#----------------------------#
# Vector of species to analyze - script loops through each
# Species without data directories will be skipped with a warning
param_species_vec    <- c("Rudd", "Goldfish", "Carp")
# Whether to search subdirectories for Excel files
param_recursive_read <- FALSE


### Define File Paths
#----------------------------#
# Root data directory
path_data_root  <- "01_data"
# Raw data location (species-specific subdirectories expected)
path_raw_root   <- file.path(path_data_root, "01_raw_files", "Species")
# Output directories (exports commented out per project conventions)
path_figs_dir   <- file.path("03_outputs", "01_figures")
path_tables_dir <- file.path("03_outputs", "01_tables")


### Initialize Output Objects
#----------------------------#
# Accumulator tibbles for cross-species results
# These collect summary stats and model coefficients as each species is processed
df_combined_summary <- tibble()
df_combined_models  <- tibble()


#####Species Processing Loop###################################----
#-------------------------------------------------------------#

# Main analysis loop - iterates through each species in param_species_vec
# For each species:
#   1. Locates and reads raw Excel files
#   2. Combines multiple files into single dataset
#   3. Cleans morphometric columns (coerces to numeric)
#   4. Calculates summary statistics
#   5. Fits regression models (linear and logarithmic)
#   6. Generates scatter plots with fitted lines and equations

for (param_species in param_species_vec) {

 cat("\n--- Processing:", param_species, "---\n")

 ### Define Species-Specific Paths
 #----------------------------#
 path_raw_dir       <- file.path(path_raw_root, param_species)
 path_processed_dir <- file.path(path_data_root, "02_processed_files", param_species)

 # Skip species if raw data directory doesn't exist
 if (!dir.exists(path_raw_dir)) { cat("WARNING: missing raw dir:", path_raw_dir, "\n"); next }


 ### Read and Combine Excel Files
 #----------------------------#
 # Find all .xlsx files in the species directory
 temp_raw_files <- list.files(path_raw_dir, "\\.xlsx$", ignore.case = TRUE,
                              full.names = TRUE, recursive = param_recursive_read)
 # Skip if no Excel files found
 if (length(temp_raw_files) == 0) { cat("WARNING: no .xlsx in:", path_raw_dir, "\n"); next }
 # Read each file and combine into single dataframe
 # Allows multiple years/batches of data to be merged automatically
 cat("Found", length(temp_raw_files), "files; combining…\n")
 temp_df_list    <- suppressMessages(lapply(temp_raw_files, readxl::read_excel))
 df_raw_combined <- dplyr::bind_rows(temp_df_list)


 #####Data Processing##########################################----
 #-------------------------------------------------------------#

 ### Clean and Validate
 #----------------------------#
 # Coerce morphometric columns to numeric type
 # suppressWarnings handles any non-numeric values (converted to NA)
 df_clean <- df_raw_combined %>%
  mutate(
   ForkLength_mm = suppressWarnings(as.numeric(ForkLength_mm)),
   Width_mm      = suppressWarnings(as.numeric(Width_mm)),
   Mass_g        = suppressWarnings(as.numeric(Mass_g))
  )

 # Verify required columns exist before proceeding with analysis
 temp_required_cols <- c("ForkLength_mm", "Width_mm", "Mass_g")
 temp_missing <- setdiff(temp_required_cols, names(df_clean))
 if (length(temp_missing)) { cat("WARNING: missing cols:", paste(temp_missing, collapse = ", "), "\n"); next }
 
 ### Summary Statistics
 #----------------------------#
 # Calculate descriptive statistics for the current species
 # Includes sample sizes (n) and means for each morphometric variable
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

 # Append to combined summary table (accumulates across species)
 df_combined_summary <- dplyr::bind_rows(df_combined_summary, df_summary)

 # Helper function to generate standardized plot captions
 make_caption <- function(df, text) paste0(param_species, " (n = ", nrow(df), "): ", text)


 #####Visualization############################################----
 #-------------------------------------------------------------#

 ### Plot 1: Width vs Fork Length (Linear Fit)
 #----------------------------#
 # Examines relationship between body width and fork length
 # Expectation: Linear relationship reflecting body proportions

 # Filter to complete cases for this analysis
 df_scatter_fl <- df_clean %>% filter(!is.na(Width_mm), !is.na(ForkLength_mm))

 # Base plot - always created
 p_scatter_fl <- ggplot(df_scatter_fl, aes(ForkLength_mm, Width_mm)) +
  geom_point(color = "#2c7fb8", alpha = 0.6, size = 2) +
  labs(
   title   = paste0(param_species, " - Width by Fork Length"),
   x       = "Fork Length (mm)",
   y       = "Width (mm)",
   caption = make_caption(df_scatter_fl, "Width vs fork length.")
  ) +
  theme_minimal()

 # Conditionally add regression line and equation if sufficient data
 if (nrow(df_scatter_fl) >= 2) {

  # Fit linear model: Width = slope * ForkLength + intercept
  temp_lm_fl  <- lm(Width_mm ~ ForkLength_mm, data = df_scatter_fl)
  temp_coef   <- coef(temp_lm_fl)
  temp_r2_fl  <- summary(temp_lm_fl)$r.squared
  temp_slope  <- unname(temp_coef[["ForkLength_mm"]])
  temp_int    <- unname(temp_coef[["(Intercept)"]])

  # Build equation string for plot annotation
  temp_eq_fl <- paste0(
   "y = ", formatC(temp_slope, format = "f", digits = 2),
   "x", ifelse(temp_int >= 0, " + ", " - "),
   formatC(abs(temp_int), format = "f", digits = 2),
   "\nR^2 = ", formatC(temp_r2_fl, format = "f", digits = 4)
  )

  # Position equation in upper-left region of plot (5th/95th percentiles)
  temp_xpos_fl <- quantile(df_scatter_fl$ForkLength_mm, 0.05, na.rm = TRUE)
  temp_ypos_fl <- quantile(df_scatter_fl$Width_mm,       0.95, na.rm = TRUE)

  # Add regression line and equation to base plot
  p_scatter_fl <- p_scatter_fl +
   geom_smooth(method = "lm", se = FALSE, color = "#1f78b4", linewidth = 0.9) +
   annotate("text", x = temp_xpos_fl, y = temp_ypos_fl, label = temp_eq_fl,
            hjust = 0, vjust = 1, size = 3.5)
 }

 # Store plot in species sublist
 # Access later via: plots$Rudd$scatter_fl, plots$Goldfish$scatter_fl, etc.
 plots[[param_species]]$scatter_fl <- p_scatter_fl

 # Optional: Save plot to file (uncomment to enable export)
 # path_file_fl <- file.path(path_figs_dir, paste0(param_species, "_scatter_width_by_forklength.png"))
 # ggsave(filename = path_file_fl, plot = p_scatter_fl, width = 7, height = 5, dpi = 300)


 ### Plot 2: Width vs Mass (Logarithmic Fit)
 #----------------------------#
 # Examines allometric relationship between body width and mass
 # Log transformation captures diminishing width gain as fish grow larger
 # Model form: Width = alpha * ln(Mass) + beta

 # Filter to complete cases with positive mass values (required for log transform)
 df_scatter_mass <- df_clean %>% filter(!is.na(Width_mm), !is.na(Mass_g), Mass_g > 0)

 # Initialize model coefficient row (filled with NA if insufficient data)
 temp_model_row <- tibble(species = param_species, alpha = NA_real_, beta = NA_real_, r2 = NA_real_)

 # Base plot - always created
 p_scatter_mass <- ggplot(df_scatter_mass, aes(Mass_g, Width_mm)) +
  geom_point(color = "#f16913", alpha = 0.6, size = 2) +
  labs(
   title   = paste0(param_species, " - Width by Mass"),
   x       = "Mass (g)",
   y       = "Width (mm)",
   caption = make_caption(df_scatter_mass, "Width vs mass.")
  ) +
  theme_minimal()

 # Conditionally add logarithmic fit if sufficient data
 if (nrow(df_scatter_mass) >= 3) {

  # Fit logarithmic model: Width = alpha * log(Mass) + beta
  temp_log_fit <- lm(Width_mm ~ log(Mass_g), data = df_scatter_mass)

  # Extract model coefficients
  temp_coefs <- coef(temp_log_fit)
  temp_alpha <- unname(temp_coefs["log(Mass_g)"])
  temp_beta  <- unname(temp_coefs["(Intercept)"])
  temp_r2_m  <- summary(temp_log_fit)$r.squared

  # Store coefficients for combined output table
  temp_model_row <- tibble(species = param_species, alpha = temp_alpha, beta = temp_beta, r2 = temp_r2_m)

  # Build equation string for plot annotation
  temp_eq_m <- paste0(
   "y = ", formatC(temp_alpha, format = "f", digits = 2),
   " ln(x)", ifelse(temp_beta >= 0, " + ", " - "),
   formatC(abs(temp_beta), format = "f", digits = 2),
   "\nR^2 = ", formatC(temp_r2_m, format = "f", digits = 4)
  )

  # Generate prediction line across full mass range
  temp_xrng <- range(df_scatter_mass$Mass_g, na.rm = TRUE)
  temp_xseq <- seq(temp_xrng[1], temp_xrng[2], length.out = 200)
  df_pred   <- tibble(Mass_g = temp_xseq) %>%
   mutate(Width_pred = predict(temp_log_fit, newdata = tibble(Mass_g = Mass_g)))

  # Position equation in upper-left region of plot
  temp_xpos_m <- quantile(df_scatter_mass$Mass_g, 0.05, na.rm = TRUE)
  temp_ypos_m <- quantile(df_scatter_mass$Width_mm, 0.95, na.rm = TRUE)

  # Add prediction line and equation to base plot
  p_scatter_mass <- p_scatter_mass +
   geom_line(data = df_pred, aes(Mass_g, Width_pred), color = "#cc4c02", linewidth = 1) +
   annotate("text", x = temp_xpos_m, y = temp_ypos_m, label = temp_eq_m,
            hjust = 0, vjust = 1, size = 3.5)
 }

 # Store plot in species sublist
 # Access later via: plots$Rudd$scatter_mass, plots$Goldfish$scatter_mass, etc.
 plots[[param_species]]$scatter_mass <- p_scatter_mass

 # Optional: Save plot to file (uncomment to enable export)
 # path_file_mass <- file.path(path_figs_dir, paste0(param_species, "_scatter_width_by_mass.png"))
 # ggsave(filename = path_file_mass, plot = p_scatter_mass, width = 7, height = 5, dpi = 300)


 ### Optional Per-Species Exports
 #----------------------------#

 # Uncomment below to save cleaned data and summary stats for each species
 # Exports are disabled by default per project file export policy

 # path_clean_rds <- file.path(path_processed_dir, paste0(param_species, "_clean.rds"))
 # saveRDS(df_clean, path_clean_rds)

 # path_species_summary_csv <- file.path(path_tables_dir, paste0(param_species, "_summary.csv"))
 # readr::write_csv(df_summary, path_species_summary_csv)


 ### Accumulate Model Results
 #----------------------------#

 # Append this species' log model coefficients to combined table
 df_combined_models <- dplyr::bind_rows(df_combined_models, temp_model_row)

 cat("Finished:", param_species, "\n")

} # End species loop

#Display all plots
print(plots)

#####Export Combined Results###################################----
#-------------------------------------------------------------#

# Uncomment below to export cross-species summary tables
# These aggregate results from all processed species into single files

# path_combined_summary_csv <- file.path(path_tables_dir, "_combined_summary_by_species.csv")
# readr::write_csv(df_combined_summary, path_combined_summary_csv)

# path_combined_models_csv  <- file.path(path_tables_dir, "_combined_log_model_coefficients.csv")
# readr::write_csv(df_combined_models,  path_combined_models_csv)


#####Cleanup###################################################----
#-------------------------------------------------------------#

# Remove all temporary objects (prefix: temp_) to keep environment clean
# This is standard practice per project conventions
temp_objects <- ls(pattern = "^temp_")
if (length(temp_objects)) { rm(list = temp_objects); cat("Removed temporary objects.\n") }

cat("\nAll done.\n")

# Objects retained in environment after script completion:
#   - df_combined_summary: Summary statistics for all species
#   - df_combined_models:  Log model coefficients (alpha, beta, R²) for all species
#   - plots:               Nested list of all plots, organized by species
#                          Access via: plots$Rudd$scatter_fl, plots$Goldfish$scatter_mass, etc.

  
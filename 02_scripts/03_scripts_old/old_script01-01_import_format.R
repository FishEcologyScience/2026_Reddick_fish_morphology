## --------------------------------------------------------------#
## Script name: script01-01_import_format.R
##
## Purpose:
##   - Read raw morphometric Excel files for each species
##   - Clean and format data into a consistent structure
##   - Save outputs to:
##         combined_all : list of cleaned data per species
##         df_all       : single combined data frame
##
## Notes:
##   - No modeling or plotting in this script
##   - Downstream scripts rely on combined_all and df_all
##
## Author: Marcus Rizzuto
## Date Created: 2026-01-27
## --------------------------------------------------------------#


##### Setup ####################################################----
#-------------------------------------------------------------#

### User parameters
param_species_vec    <- c("Rudd", "Goldfish", "Carp")
param_recursive_read <- FALSE   # read subfolders? usually FALSE

### Define paths
path_data_root <- "01_data"
path_raw_root  <- file.path(path_data_root, "01_raw_files", "Species")

### Initialize outputs
combined_all <- list()   # cleaned data stored by species


##### Species Loop ############################################----
#-------------------------------------------------------------#

for (param_species in param_species_vec) {
 
 cat("\n--- Processing:", param_species, "---\n")
 
 # Path to raw files for this species
 path_raw_dir <- file.path(path_raw_root, param_species)
 
 # Skip if folder is missing
 if (!dir.exists(path_raw_dir)) {
  cat("WARNING: missing directory:", path_raw_dir, "\n")
  next
 }
 
 ### Locate Excel files
 temp_raw_files <- list.files(
  path_raw_dir,
  pattern = "\\.xlsx$", ignore.case = TRUE,
  full.names = TRUE,
  recursive = param_recursive_read
 )
 
 if (length(temp_raw_files) == 0) {
  cat("WARNING: no .xlsx files found in:", path_raw_dir, "\n")
  next
 }
 
 cat("Found", length(temp_raw_files), "files; combining…\n")
 
 # Read all species files into a list, suppressing readxl messages
 temp_df_list <- suppressMessages(lapply(temp_raw_files, readxl::read_excel))
 
 # Combine into a single data frame
 df_raw_combined <- dplyr::bind_rows(temp_df_list)
 
 
 ##### Cleaning & Formatting ##################################----
 #-------------------------------------------------------------#
 
 # Convert numeric morphometrics; coerce quietly to NA if needed
 df_clean <- df_raw_combined %>%
  dplyr::mutate(
   ForkLength_mm = suppressWarnings(as.numeric(ForkLength_mm)),
   Width_mm      = suppressWarnings(as.numeric(Width_mm)),
   Mass_g        = suppressWarnings(as.numeric(Mass_g))
  )
 
 # Ensure required columns exist
 required_cols <- c("ForkLength_mm", "Width_mm", "Mass_g")
 missing_cols  <- setdiff(required_cols, names(df_clean))
 
 if (length(missing_cols)) {
  cat("WARNING: missing columns:", paste(missing_cols, collapse = ", "), "\n")
  next
 }
 
 # Ensure a Species column exists (overwrite with param_species)
 df_clean$Species <- param_species
 
 # Store cleaned dataset
 combined_all[[param_species]] <- df_clean
 
 cat("Imported rows:", nrow(df_clean), "\n")
}


##### Combine all species ####################################----
#-------------------------------------------------------------#

df_all <- dplyr::bind_rows(combined_all, .id = "species_name")

cat("\nAll species imported. Total rows:", nrow(df_all), "\n")

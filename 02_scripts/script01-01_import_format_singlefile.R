## --------------------------------------------------------------#
## Script name: script01-01_import_format_singlefile.R
##
## Purpose:
##   - Import morphometric data from ONE file containing ALL species
##   - Clean/format the dataset
##   - Split into:
##         combined_all : list of cleaned data per species
##         df_all       : combined tibble of all species
##
## Notes:
##   - Keeps same output objects as the multi-file workflow
##   - No modeling or plotting in this script
##   - Downstream scripts rely on combined_all and df_all
##
## Author: Marcus Rizzuto
## Date Created: 2026-01-27
## --------------------------------------------------------------#


##### Setup ####################################################----
#-------------------------------------------------------------#

# Path to the single Excel file containing all species
path_data_root <- "01_data"
path_singlefile <- file.path(
 path_data_root,
 "01_raw_files",
 "Fish widths_DraftExcelAnalysis_Jan2026.xlsx"
)

# Output objects
combined_all <- list()


##### Import Single Dataset ###################################----
#-------------------------------------------------------------#

cat("\n--- Importing single-table dataset ---\n")

# Read file
df_raw <- suppressMessages(readxl::read_excel(path_singlefile))

cat("Rows imported:", nrow(df_raw), "\n")


##### Cleaning & Formatting ##################################----
#-------------------------------------------------------------#

# Required columns
required_cols <- c("Species", "ForkLength_mm", "Width_mm", "Mass_g")
missing_cols  <- setdiff(required_cols, names(df_raw))

if (length(missing_cols)) {
 stop("ERROR: Missing required columns: ",
      paste(missing_cols, collapse = ", "))
}

# Clean and coerce to numeric
df_all <- df_raw %>%
 dplyr::mutate(
  Species       = as.character(Species),
  ForkLength_mm = suppressWarnings(as.numeric(ForkLength_mm)),
  Width_mm      = suppressWarnings(as.numeric(Width_mm)),
  Mass_g        = suppressWarnings(as.numeric(Mass_g))
 )


##### Split by species #######################################----
#-------------------------------------------------------------#

combined_all <- split(df_all, df_all$Species)

# Sanitize species names if needed
sp_names <- names(combined_all)
bad_idx  <- which(is.na(sp_names) | !nzchar(sp_names))
if (length(bad_idx)) {
 sp_names[bad_idx] <- paste0("UNKNOWN_", bad_idx)
 names(combined_all) <- sp_names
}

cat("Species detected:", paste(names(combined_all), collapse = ", "), "\n")
cat("Total rows:", nrow(df_all), "\n")
cat("\nSingle-table import complete.\n")




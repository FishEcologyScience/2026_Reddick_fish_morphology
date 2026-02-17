## --------------------------------------------------------------#
## Script name: script01_02_elecfish_length_histograms_2012on.R
##
## Purpose:
##   - Read raw ElecFish export (2012+; pre-filtered upstream)
##   - Clean and standardize species + length fields
##   - Generate fork‑length histograms per species
##   - Export QC tables (species counts, min/max lengths)
##
#### Repository layout (per template):
##   01_data/01_raw_files/          # raw ElecFish Excel export
##   03_outputs/01_figures/         # figures
##   03_outputs/01_tables/          # tables
##
## Author: Marcus Rizzuto
## Date Created: 2/4/2026
##
## Notes:
##   - Single-use script (histograms only)
##
## --------------------------------------------------------------#

##### Setup #############################################################
#-----------------------------------------------------------------------#
suppressPackageStartupMessages({
 library(readxl)
 library(dplyr)
 library(stringr)
 library(ggplot2)
 library(tidyr)
 library(readr)
})

##### CONFIG — UPDATE THESE TO MATCH YOUR FILE ##########################
#-----------------------------------------------------------------------#
# Path to the ElecFish Excel export (already filtered to 2012+ upstream)
path_elecfish_xlsx <- file.path("01_data", "01_raw_files", "lw_raw_elecfish_2012_on.xlsx")

# If data are not on the first sheet, set a name or index here (else keep NULL)
excel_sheet_name <- NULL  # e.g., "Export" or 1

# Column names in your Excel export (match EXACTLY to your headers)
COL_SPECIES   <- "Species"        # species/common/scientific name (pick one)
COL_LENGTH_MM <- "Length"         
# Date is not required here since we’re not filtering by year, but keep if present
COL_DATE      <- NULL             # e.g., "SampleDate" or NULL

# Minimum observations required to plot a species histogram
MIN_N_PER_SPECIES <- 5

# Outputs
out_dir_figs   <- file.path("03_outputs", "01_figures")  
out_dir_tables <- file.path("03_outputs", "01_tables")   

# Histogram defaults
DEFAULT_BINS   <- 30
FD_FALLBACK_MM <- 5               # fallback binwidth (mm) if FD binwidth can’t be computed

##### Create output directories ########################################
#-----------------------------------------------------------------------#
dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

##### Read Excel #######################################################
#-----------------------------------------------------------------------#
message("Reading: ", normalizePath(path_elecfish_xlsx, mustWork = FALSE))
df_raw <- if (is.null(excel_sheet_name)) {
 readxl::read_excel(path_elecfish_xlsx)
} else {
 readxl::read_excel(path_elecfish_xlsx, sheet = excel_sheet_name)
}

# Sanity check for required columns
required_cols <- c(COL_SPECIES, COL_LENGTH_MM)
missing_cols  <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
 stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

##### Normalize (no year filter needed) ################################
#-----------------------------------------------------------------------#
df <- df_raw %>%
 # Standardize working column names
 rename(
  Species_raw = !!sym(COL_SPECIES),
  Length_raw  = !!sym(COL_LENGTH_MM)
 ) %>%
 mutate(
  Species   = Species_raw %>% as.character() %>% stringr::str_squish(),
  Length_mm = suppressWarnings(as.numeric(Length_raw))
 ) %>%
 filter(!is.na(Length_mm), Length_mm > 0)

# Optional: if a date column exists and you want it preserved (not used for filtering)
if (!is.null(COL_DATE) && COL_DATE %in% names(df_raw)) {
 df <- df %>% mutate(Date = .data[[COL_DATE]])
}

# Quick QC counts by species
df_counts <- df %>%
 group_by(Species) %>%
 summarise(
  n_all = n(),
  min_length_mm = suppressWarnings(min(Length_mm, na.rm = TRUE)),
  max_length_mm = suppressWarnings(max(Length_mm, na.rm = TRUE)),
  .groups = "drop"
 ) %>%
 arrange(desc(n_all))

out_counts <- file.path(out_dir_tables, "elecfish_counts_by_species_2012plus.csv")
readr::write_csv(df_counts, out_counts)
message("Wrote counts table: ", out_counts)

##### Helper: Freedman–Diaconis binwidth with fallback #################
#-----------------------------------------------------------------------#
fd_binwidth <- function(x, fallback = FD_FALLBACK_MM) {
 x <- x[is.finite(x)]
 n <- length(x)
 if (n <= 1) return(fallback)
 iqr <- IQR(x, na.rm = TRUE)
 bw  <- 2 * iqr / (n)^(1/3)
 if (is.finite(bw) && bw > 0) bw else fallback
}

##### Per-species histograms ###########################################
#-----------------------------------------------------------------------#
species_vec <- sort(unique(df$Species[!is.na(df$Species) & nzchar(df$Species)]))

if (length(species_vec) == 0) {
 warning("No species found with valid positive lengths.")
} else {
 
 ##### PLOTTING LOOP — per-species length histograms #################
 #-------------------------------------------------------------------#
 for (loop_species in species_vec) {
  loop_df <- dplyr::filter(df, Species == loop_species)
  if (nrow(loop_df) < MIN_N_PER_SPECIES) {
   message("Skipping ", loop_species, " (n = ", nrow(loop_df), " < ", MIN_N_PER_SPECIES, ").")
   next
  }
  
  loop_bw       <- fd_binwidth(loop_df$Length_mm)
  loop_fl_mean   <- mean(loop_df$Length_mm, na.rm = TRUE)
  loop_fl_median <- stats::median(loop_df$Length_mm, na.rm = TRUE)
  
  loop_p  <- ggplot(loop_df, aes(x = Length_mm)) +
   geom_histogram(binwidth = loop_bw, boundary = 0,
                  color = "grey30", fill = "#74a9cf", alpha = 0.85) +
   geom_vline(xintercept = loop_fl_mean,   color = "#de2d26", linewidth = 0.7) +
   geom_vline(xintercept = loop_fl_median, color = "#238b45", linetype = "dashed", linewidth = 0.7) +
   labs(
    title   = paste0(loop_species, " — Length histogram (2012+)"),
    x       = "Length (mm)",
    y       = "Count",
    caption = paste0(
     loop_species, " (n = ", nrow(loop_df), "): ",
     "Histogram of lengths. ",
     "Red = mean (", formatC(loop_fl_mean, digits = 1, format = "f"),
     " mm), green dashed = median (", formatC(loop_fl_median, digits = 1, format = "f"), " mm)."
    )
   ) +
   theme_minimal(base_size = 12)
  
  # Show the per-species plot in the RStudio Plots pane
  print(loop_p)
  
  # Optional: briefly pause so plots “tick by” when sourcing
  # Sys.sleep(0.15)
  
  
  # Safe filename for species (saving commented out)
  loop_sp_file <- loop_species %>%
   stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
   stringr::str_replace("^_+|_+$", "")
  
  loop_out_png <- file.path(out_dir_figs, paste0("hist_length_", loop_sp_file, "_2012plus.png"))
  # ggsave(loop_out_png, loop_p, width = 7, height = 5, dpi = 300)
  # message("Saved: ", loop_out_png)
 }
}

##### Optional: combined faceted panel ###############################
#-----------------------------------------------------------------------#
if (nrow(df) > 0) {
 p_all <- ggplot(df, aes(x = Length_mm)) +
  geom_histogram(color = "grey25", fill = "#9ecae1", bins = DEFAULT_BINS, boundary = 0) +
  labs(
   title = "ElecFish lengths by species (2012+)",
   x = "Length (mm)", y = "Count"
  ) +
  facet_wrap(~ Species, scales = "free_y") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
 
 # Show the combined panel in the RStudio Plots pane
 print(p_all)
 
 # Optional: brief pause to view the combined panel
 # Sys.sleep(0.25)
 
 out_png_all <- file.path(out_dir_figs, "hist_length_by_species_2012plus.png")
 # ggsave(out_png_all, p_all, width = 10, height = 8, dpi = 300)
 # message("Saved combined panel: ", out_png_all)
}

message("\nDone. Inspect figures in: ", normalizePath(out_dir_figs, mustWork = FALSE))
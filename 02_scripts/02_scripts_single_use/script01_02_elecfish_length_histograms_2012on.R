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
##   - Goldfish (F0181) were sampled from June 19, 2012 through October 28, 2025 across spring, summer, and fall with summer most sampled (672 records; 54.7%), followed by fall (540; 44.0%) and spring (16; 1.3%). The full dataset contains 1,228 total goldfish records.  Sampling spans 39 unique transects, among rows with valid Time Period codes, effort was predominantly night (code 2: 922; 89.2%) vs. day (code 1: 112; 10.8%), with no crepuscular codes present and 194 rows missing a time‑period code.
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
 library(lubridate)  # robust date parsing
})

##### CONFIG —  #########################################################
#-----------------------------------------------------------------------#
# Path to the ElecFish Excel export (already filtered to 2012+ upstream)
path_elecfish_xlsx <- file.path("01_data", "01_raw_files", "fixed_lw_raw_elecfish_2012_on.xlsx")

# If data are not on the first sheet, set a name or index here (else keep NULL)
excel_sheet_name <- NULL  # e.g., "Export" or 1

# Column names in your Excel export (match EXACTLY to your headers)
COL_SPECIES   <- "Species"        # species/common/scientific name (pick one)
COL_LENGTH_MM <- "Length"
# If NULL, we will auto-detect the first column containing "date" (case-insensitive)
COL_DATE      <- NULL             # e.g., "SampleDate" or NULL

# Path to species lookup (two columns: Species, Common Name)
path_species_lookup <- file.path("01_data", "01_raw_files", "species_lookup.xlsx")
COL_SPECIES_ID      <- "Species"
COL_COMMON_NAME     <- "Common Name"

# Minimum observations required to plot a species histogram
MIN_N_PER_SPECIES <- 5

# Outputs
out_dir_figs   <- file.path("03_outputs", "01_figures")
out_dir_tables <- file.path("03_outputs", "01_tables")

# Histogram defaults
DEFAULT_BINS   <- 30
FD_FALLBACK_MM <- 5               # fallback binwidth (mm) if FD binwidth can’t be computed

# Season palette (clean, no purple/orange/green combo)
PAL_SEASON <- c(Spring = "#A1D99B",  # teal
                Summer = "#FDD0A2",  # blue
                Fall   = "#9ECAE1")  # brown

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
# Auto-detect a date column if COL_DATE is NULL
if (is.null(COL_DATE)) {
 date_candidates <- names(df_raw)[grepl("date", names(df_raw), ignore.case = TRUE)]
 if (length(date_candidates) >= 1) {
  COL_DATE <- date_candidates[1]
  message("Auto-detected date column: ", COL_DATE)
 }
}

# Build df (carry Date before filtering so rows stay aligned)
df <- df_raw %>%
 # Standardize working column names
 rename(
  Species_raw = !!rlang::sym(COL_SPECIES),
  Length_raw  = !!rlang::sym(COL_LENGTH_MM)
 ) %>%
 mutate(
  Species   = Species_raw %>% as.character() %>% stringr::str_squish(),
  Length_mm = suppressWarnings(as.numeric(Length_raw))
 )

# ---- Sanitize Species (prevents NA/blank species from being dropped downstream) ----
df <- df %>%
 mutate(
  Species = dplyr::case_when(
   is.na(Species) ~ "UNKNOWN",
   !nzchar(Species) ~ "UNKNOWN",
   TRUE ~ Species
  )
 )

# Safety checks — should be true after sanitization
stopifnot(!any(is.na(df$Species) | !nzchar(df$Species)))

if (!is.null(COL_DATE) && COL_DATE %in% names(df_raw)) {
 df <- df %>% mutate(Date_raw = .data[[COL_DATE]])
}

# Now filter invalid lengths
df <- df %>% filter(!is.na(Length_mm), Length_mm > 0)

# --- MONTH + SEASON DERIVATION (from Date; fallback to Month) ---------
if ("Date_raw" %in% names(df)) {
 # Parse Date robustly (supports common formats)
 df <- df %>%
  mutate(
   Date_parsed = suppressWarnings(lubridate::parse_date_time(
    as.character(Date_raw),
    orders = c("Ymd","Y-m-d","mdY","m/d/Y","dmy","d/m/Y","Y/m/d","m-d-Y","d-m-Y"),
    quiet  = TRUE
   )),
   Date  = as.Date(coalesce(Date_parsed, suppressWarnings(as.Date(Date_raw)))),
   Month = suppressWarnings(lubridate::month(Date))
  ) %>% select(-Date_parsed)
} else if ("Month" %in% names(df)) {
 df <- df %>% mutate(Month = suppressWarnings(as.integer(.data[["Month"]])))
} else {
 message("No Date or Month column found — will use non-season colored fallback.")
}

# Map Month -> Season per requirement:
#   Spring = May–June (5,6)
#   Summer = July–August (7,8)
#   Fall   = everything else (including NA)
if ("Month" %in% names(df)) {
 df <- df %>%
  mutate(
   Season = dplyr::case_when(
    Month %in% c(5L, 6L) ~ "Spring",
    Month %in% c(7L, 8L) ~ "Summer",
    TRUE                 ~ "Fall"
   ),
   Season = factor(Season, levels = c("Spring","Summer","Fall"))
  )
}

# --- SPECIES LOOKUP + LABEL CREATION ---------------------------------

message("Reading species lookup: ", normalizePath(path_species_lookup, mustWork = FALSE))

lookup_raw <- readxl::read_excel(path_species_lookup)

lookup_missing <- setdiff(c(COL_SPECIES_ID, COL_COMMON_NAME), names(lookup_raw))
if (length(lookup_missing) > 0) {
 stop("Lookup is missing required column(s): ", paste(lookup_missing, collapse = ", "))
}

lookup <- lookup_raw %>%
 transmute(
  Species_id  = !!rlang::sym(COL_SPECIES_ID) %>% as.character() %>% stringr::str_squish(),
  Common_Name = !!rlang::sym(COL_COMMON_NAME) %>% as.character() %>% stringr::str_squish()
 ) %>%
 filter(!is.na(Species_id), nzchar(Species_id)) %>%
 distinct(Species_id, .keep_all = TRUE)

labvec <- setNames(
 paste0(lookup$Species_id, " — ", lookup$Common_Name),
 lookup$Species_id
)

# Also build a named vector for just the common name (names = Species_id)
lab_common <- setNames(lookup$Common_Name, lookup$Species_id)

# Helper to make file-system safe names from labels/common names
make_safe_filename <- function(x, max_chars = 150) {
 x %>%
  stringr::str_to_lower() %>%
  stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
  stringr::str_replace_all("_+", "_") %>%
  stringr::str_replace("^_+|_+$", "") %>%
  stringr::str_sub(1, max_chars)
}

df <- df %>%
 mutate(
  Species_label = dplyr::coalesce(labvec[Species], Species) %>% as.character()
 )

# ----------------------------------------------------------------------

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

# out_counts <- file.path(out_dir_tables, "elecfish_counts_by_species_2012plus.csv")
# readr::write_csv(df_counts, out_counts)
# message("Wrote counts table: ", out_counts)

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
species_vec <- sort(unique(df$Species))

if (length(species_vec) == 0) {
 warning("No species found with valid positive lengths.")
} else {
 
 ##### PLOTTING LOOP — per-species length histograms #################
 #-------------------------------------------------------------------#
 for (loop_species in species_vec) {
  loop_df    <- dplyr::filter(df, Species == loop_species)
  loop_label <- dplyr::first(loop_df$Species_label)
  if (nrow(loop_df) < MIN_N_PER_SPECIES) {
   message("Skipping ", loop_species, " (n = ", nrow(loop_df), " < ", MIN_N_PER_SPECIES, ").")
   next
  }
  
  loop_bw        <- fd_binwidth(loop_df$Length_mm)
  loop_fl_mean   <- mean(loop_df$Length_mm, na.rm = TRUE)
  loop_fl_median <- stats::median(loop_df$Length_mm, na.rm = TRUE)
  
  # Build plot — single panel, clearly color-coded by Season (stacked)
  has_season <- "Season" %in% names(loop_df) && any(!is.na(loop_df$Season))
  
  if (has_season) {
   season_counts <- loop_df %>%
    dplyr::filter(!is.na(Season)) %>%
    dplyr::count(Season, name = "n")
   sc_txt <- if (nrow(season_counts) > 0) {
    paste(paste0(season_counts$Season, "=", season_counts$n), collapse = ", ")
   } else {
    "no seasonal breakdown"
   }
   
   loop_p <- ggplot(loop_df, aes(x = Length_mm, fill = Season)) +
    geom_histogram(
     binwidth  = loop_bw,
     boundary  = 0,
     position  = "stack",   # stacked for clear season contributions
     alpha     = 0.95,
     color     = "grey20"
    ) +
    geom_vline(xintercept = loop_fl_mean,   color = "#de2d26", linewidth = 0.7) +
    geom_vline(xintercept = loop_fl_median, color = "#238b45", linetype = "dashed", linewidth = 0.7) +
    scale_fill_manual(
     values = c(
      Spring = "#66C2A5",   # pastel green
      Summer = "#FC8D62",   # pastel orange
      Fall   = "#8DA0CB"    # pastel light blue
     ),
     drop = FALSE
    ) +
    labs(
     title   = paste0(loop_label, " — Length histogram by season (2012+)"),
     caption = paste0(
      loop_label, " (n = ", nrow(loop_df), "; ", sc_txt, "). ",
      "Red = mean (", formatC(loop_fl_mean, digits = 1, format = "f"),
      " mm), green dashed = median (", formatC(loop_fl_median, digits = 1, format = "f"), " mm)."
     ),
     x = "Fork length (mm)",
     y = "Count",
     fill = "Season"
    ) +
    theme_minimal(base_size = 12)
   
  } else {
   # Fallback: original single-panel histogram (no season available)
   loop_p <- ggplot(loop_df, aes(x = Length_mm)) +
    geom_histogram(
     binwidth = loop_bw,
     boundary = 0,
     color = "grey30",
     fill  = "#74a9cf",
     alpha = 0.85
    ) +
    geom_vline(xintercept = loop_fl_mean,   color = "#de2d26", linewidth = 0.7) +
    geom_vline(xintercept = loop_fl_median, color = "#238b45", linetype = "dashed", linewidth = 0.7) +
    labs(
     title   = paste0(loop_label, " — Length histogram (2012+)"),
     caption = paste0(
      loop_label, " (n = ", nrow(loop_df), "): ",
      "Histogram of lengths. ",
      "Red = mean (", formatC(loop_fl_mean, digits = 1, format = "f"),
      " mm), green dashed = median (", formatC(loop_fl_median, digits = 1, format = "f"), " mm)."
     ),
     x = "Fork length (mm)",
     y = "Count"
    ) +
    theme_minimal(base_size = 12)
  }
  
  # Show the per-species plot in the RStudio Plots pane
  print(loop_p)
  
  # Optional: briefly pause so plots “tick by” when sourcing
  # Sys.sleep(0.15)
  
  # Safe filename for species (saving commented out)
  # Filename based on both ID and common name (avoids collisions)
  loop_common  <- dplyr::coalesce(lab_common[loop_species], loop_species) %>% as.character()
  loop_sp_file <- make_safe_filename(paste(loop_species, loop_common, sep = "_"))
  loop_out_png <- file.path(out_dir_figs, paste0("hist_length_", loop_sp_file, "_2012plus.png"))
  # ggsave(loop_out_png, loop_p, width = 7, height = 5, dpi = 300)
 }
}

message("\nDone. Inspect figures in: ", normalizePath(out_dir_figs, mustWork = FALSE))

## --------------------------------------------------------------#
## Script name: script02_rlf_largemouth_bass_bonar_bins.R
##
## Purpose:
##   - Read ElecFish export
##   - Filter to Largemouth Bass (LMB) and valid positive lengths
##   - Compute Relative Length Frequency (RLF) using Bonar (2002) bins
##     * Exclude <100 mm (age-0) for LMB per Bonar
##     * Bins: 100–200, 201–300, 301–375, 376–450, >450 mm
##   - Export counts + percentages by bin (CSV)
##   - Plot population-only bar chart, or population vs ALF (if provided)
##
#### Repository layout (per template):
##   2026_fish_morphology/
##     01_data/01_raw_files/          # raw ElecFish Excel export
##     01_data/02_derived/            # (optional) ALF CSV if available
##     03_outputs/01_figures/         # figures
##     03_outputs/01_tables/          # tables
##
## Author: Marcus Rizzuto
## Date Created: 2026-02-20
##
## Notes:
##   - Method follows Bonar (2002) RLF: exclude age-0 (>=100 mm for LMB)
##     and use five fixed bins for visual comparison to an Average Length
##     Frequency (ALF) built from many waterbodies with consistent gear/time).
##   - Default LMB species ID set to 'F0317' based on your export.
##   - Add an ALF CSV to show side-by-side RLF vs ALF bars.
##
## References:
##   - Bonar SA (2002) Relative Length Frequency: A Simple, Visual Technique
##     to Evaluate Size Structure in Fish Populations. NAJFM 22:1086–1094.
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

##### CONFIG — ##########################################################
#-----------------------------------------------------------------------#
# Paths relative to the repo root (recommended)
path_elecfish_xlsx <- file.path("01_data", "01_raw_files", "LMB_summary.xlsx")

# Fail fast if the file path is wrong
stopifnot(file.exists(path_elecfish_xlsx))

# If data are not on the first sheet, set a name or index here (else keep NULL)
excel_sheet_name <- NULL  # e.g., "Export" or 1

# Column names in your Excel export 
COL_SPECIES   <- "Species"
COL_LENGTH_MM <- "Length"

# Species filter (Largemouth Bass identifier)
LMB_ID <- "F0317"  # change if your file uses a different code

# Optional: path to an ALF CSV with Bonar bins (columns: Bin, Percent)
# Bin values must match: c("100–200","201–300","301–375","376–450",">450")
path_alf_csv <- NULL  # e.g., file.path("01_data", "02_derived", "alf_lmb_bonar.csv")

# Outputs (project-root-relative)
out_dir_figs   <- file.path("03_outputs", "01_figures")
out_dir_tables <- file.path("03_outputs", "01_tables")

# Bonar-style LMB bins (>= 100 mm)
bonar_breaks <- c(100, 200, 300, 375, 450, Inf)
bonar_labels <- c("100–200", "201–300", "301–375", "376–450", ">450")

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

##### Normalize (filter to LMB, clean length) ##########################
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
 filter(
  Species == LMB_ID,
  !is.na(Length_mm),
  Length_mm > 0
 )

message("LMB rows (Length > 0): ", nrow(df))

##### Exclude <100 mm for RLF (Bonar LMB convention) ###################
#-----------------------------------------------------------------------#
df_ge100 <- df %>% filter(Length_mm >= 100)
df_lt100 <- df %>% filter(Length_mm <  100)

message("LMB rows used for RLF (>=100 mm): ", nrow(df_ge100))
message("LMB rows <100 mm (excluded from RLF): ", nrow(df_lt100))

##### Bin assignment (Bonar bins) ######################################
#-----------------------------------------------------------------------#
df_binned <- df_ge100 %>%
 mutate(
  Bin = cut(
   Length_mm,
   breaks = bonar_breaks,
   include_lowest = TRUE,
   right = TRUE,
   labels = bonar_labels
  ),
  Bin = factor(Bin, levels = bonar_labels) # enforce order
 )

##### Counts & Relative Frequency ######################################
#-----------------------------------------------------------------------#
tab_pop <- df_binned %>%
 count(Bin, name = "Count") %>%
 tidyr::complete(Bin = factor(bonar_labels, levels = bonar_labels), fill = list(Count = 0L)) %>%
 mutate(
  N_total = sum(Count),
  Percent = 100 * Count / N_total
 )

# Export table
out_counts <- file.path(out_dir_tables, "rlf_lmb_counts_percent_bonar.csv")
readr::write_csv(tab_pop, out_counts)
message("Wrote counts table: ", out_counts)

##### Optional: load ALF (for side-by-side comparison) #################
#-----------------------------------------------------------------------#
alf_available <- !is.null(path_alf_csv) && file.exists(path_alf_csv)
if (alf_available) {
 alf <- readr::read_csv(path_alf_csv, show_col_types = FALSE) %>%
  mutate(
   Bin    = factor(Bin, levels = bonar_labels),
   Source = "ALF"
  ) %>%
  select(Bin, Percent, Source)
 
 pop <- tab_pop %>%
  transmute(Bin, Percent, Source = "Your population")
 
 to_plot <- dplyr::bind_rows(alf, pop)
} else {
 to_plot <- tab_pop %>%
  mutate(Source = "Your population") %>%
  select(Bin, Percent, Source)
}

##### Plotting #########################################################
#-----------------------------------------------------------------------#
if (alf_available) {
 p <- ggplot(to_plot, aes(x = Bin, y = Percent, fill = Source)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey25") +
  geom_text(aes(label = sprintf("%.2f%%", Percent)),
            position = position_dodge(width = 0.75), vjust = -0.25, size = 3.1) +
  scale_fill_manual(values = c("Your population" = "#2b8cbe", "ALF" = "#a6bddb")) +
  labs(
   title    = "Largemouth Bass — Relative Length Frequency (Bonar bins, ≥100 mm)",
   subtitle = "Side-by-side with ALF (if provided)",
   x = "Length bin (mm)", y = "Percentage of sample", fill = NULL
  ) +
  theme_minimal(base_size = 12)
 out_png <- file.path(out_dir_figs, "rlf_lmb_vs_alf.png")
} else {
 p <- ggplot(to_plot, aes(x = Bin, y = Percent)) +
  geom_col(fill = "#2b8cbe", colour = "grey25") +
  geom_text(aes(label = sprintf("%.2f%%", Percent)), vjust = -0.25, size = 3.1) +
  labs(
   title    = "Largemouth Bass — Relative Length Frequency (Bonar bins, ≥100 mm)",
   subtitle = "Your population only (add an ALF CSV to compare)",
   x = "Length bin (mm)", y = "Percentage of sample"
  ) +
  theme_minimal(base_size = 12)
 out_png <- file.path(out_dir_figs, "rlf_lmb_population_only.png")
}

print(p)
ggsave(out_png, p, width = 8, height = 5.5, dpi = 300)
message("Wrote figure: ", out_png)

##### (Optional) Skewness g1 (quick diagnostic) ########################
#-----------------------------------------------------------------------#
# Simple moment-based skewness (computed on >=100 mm lengths).
g1 <- function(x) {
 x <- x[is.finite(x)]
 n  <- length(x); if (n < 3) return(NA_real_)
 xb <- mean(x)
 m2 <- mean((x - xb)^2)
 m3 <- mean((x - xb)^3)
 if (m2 == 0) return(NA_real_)
 m3 / (m2^(3/2))
}

g1_pop <- g1(df_ge100$Length_mm)
message(sprintf("Skewness g1 (your >=100 mm lengths): %.3f", g1_pop))

# If ALF provided, simulate an ALF pseudo-sample at bin midpoints to approximate g1
if (alf_available) {
 bin_mid <- c("100–200"=150, "201–300"=250, "301–375"=338, "376–450"=413, ">450"=475)
 Nsim    <- sum(tab_pop$Count)
 alf_sim <- to_plot %>%
  dplyr::filter(Source == "ALF") %>%
  dplyr::mutate(n_bin = round(Percent/100 * Nsim)) %>%
  dplyr::rowwise() %>%
  dplyr::summarise(vals = list(rep(bin_mid[as.character(Bin)], n_bin)), .groups = "drop") %>%
  tidyr::unnest(cols = vals) %>% dplyr::pull(vals)
 g1_alf <- g1(alf_sim)
 message(sprintf("Skewness g1 (ALF pseudo-sample): %.3f", g1_alf))
}

message("\nDone. Inspect tables in: ", normalizePath(out_dir_tables, mustWork = FALSE))
message("Done. Inspect figures in: ", normalizePath(out_dir_figs,   mustWork = FALSE))
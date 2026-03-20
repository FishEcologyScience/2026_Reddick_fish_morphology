## --------------------------------------------------------------#
## Script name: script01_06_plot_height_length_width.R
##
## Purpose:
##   - Import morphometric data from a single CSV:
##       "WL_GF_30May_updated.csv"
##   - Produce two plots for Goldfish:
##       (1) Height vs Fork Length (FL only; no fallback)
##       (2) Height vs Width
##   - Save PNGs only to 03_outputs/01_figures
##
## Author: Marcus Rizzuto
## Date Created: 2026-03-17
## --------------------------------------------------------------#


##### Setup ####################################################----
#-------------------------------------------------------------#

suppressPackageStartupMessages({
 library(readr)
 library(dplyr)
 library(ggplot2)
})

# Paths
path_singlefile <- file.path("01_data", "01_raw_files", "WL_GF_30May_updated.csv")

# Output directory (PNG only)
path_out_dir <- file.path("03_outputs", "01_figures")
if (!dir.exists(path_out_dir)) dir.create(path_out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n--- Plotting: Goldfish Height~FL and Height~Width (CSV input) ---\n")


##### Import ##################################################----
#-------------------------------------------------------------#

df_raw <- readr::read_csv(path_singlefile, show_col_types = FALSE)
cat("Rows imported:", nrow(df_raw), "\n")


##### Cleaning (direct; FL only, no TL) ######################----
#-------------------------------------------------------------#

df <- df_raw %>%
 mutate(
  FL_mm     = suppressWarnings(as.numeric(FL_mm)),
  width_mm  = suppressWarnings(as.numeric(width_mm)),
  height_mm = suppressWarnings(as.numeric(height_mm))
 ) %>%
 select(FL_mm, width_mm, height_mm)

# Availability for each plot
n_hl <- df %>% filter(!is.na(FL_mm),    !is.na(height_mm)) %>% nrow()
n_hw <- df %>% filter(!is.na(width_mm), !is.na(height_mm)) %>% nrow()
cat("Rows available for Height~FL   :", n_hl, "\n")
cat("Rows available for Height~Width:", n_hw, "\n")


##### Helpers: line+equation plot, save PNG ##################----
#-------------------------------------------------------------#

plot_with_lm_eq <- function(data, xvar, yvar, xlab, ylab, title) {
 # Build y ~ x formula safely from strings
 fit <- lm(reformulate(xvar, response = yvar), data = data)
 
 coefs     <- coef(fit)
 slope     <- unname(coefs[2])
 intercept <- unname(coefs[1])
 r2        <- summary(fit)$r.squared
 
 # Equation text
 eq_txt <- paste0(
  "y = ", formatC(slope, format = "f", digits = 2), "x",
  ifelse(intercept >= 0, " + ", " - "), formatC(abs(intercept), format = "f", digits = 2),
  "\nR^2 = ", formatC(r2, format = "f", digits = 4)
 )
 
 # Nice placement (top-left corner within data range)
 xr <- range(data[[xvar]], na.rm = TRUE)
 yr <- range(data[[yvar]], na.rm = TRUE)
 x_anno <- xr[1] + 0.05 * diff(xr)
 y_anno <- yr[2] - 0.05 * diff(yr)
 
 ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
  geom_point(color = "#2C7FB8", alpha = 0.7, size = 2) +
  # Single regression line (no grey band)
  geom_abline(slope = slope, intercept = intercept, color = "#1f78b4", linewidth = 1) +
  annotate("text", x = x_anno, y = y_anno, label = eq_txt,
           hjust = 0, vjust = 1, size = 3.5) +
  labs(title = title, x = xlab, y = ylab) +
  theme_minimal(base_size = 12) +
  theme(
   plot.title = element_text(face = "bold"),
   panel.grid.minor = element_blank()
  )
}

save_png <- function(plot, filename, w = 8, h = 6, dpi = 300) {
 outfile <- file.path(path_out_dir, paste0(filename, ".png"))
 ggsave(outfile, plot, width = w, height = h, dpi = dpi)
 cat("Saved:", outfile, "\n")
}


##### Build EXACTLY TWO plots #################################----
#-------------------------------------------------------------#

# 1) Height vs Fork Length (FL only)
if (n_hl >= 2) {
 d_hl <- df %>% filter(!is.na(FL_mm), !is.na(height_mm))
 p_hl <- plot_with_lm_eq(
  data  = d_hl,
  xvar  = "FL_mm",
  yvar  = "height_mm",
  xlab  = "Fork Length (mm)",
  ylab  = "Height (mm)",
  title = paste0("Goldfish (n = ", nrow(d_hl), "): Height vs Fork Length (FL)")
 )
 save_png(p_hl, "goldfish_height_vs_FL")
} else {
 warning("Not enough points for Height ~ FL (need >= 2).")
}

# 2) Height vs Width
if (n_hw >= 2) {
 d_hw <- df %>% filter(!is.na(width_mm), !is.na(height_mm))
 p_hw <- plot_with_lm_eq(
  data  = d_hw,
  xvar  = "width_mm",
  yvar  = "height_mm",
  xlab  = "Width (mm)",
  ylab  = "Height (mm)",
  title = paste0("Goldfish (n = ", nrow(d_hw), "): Height vs Width")
 )
 save_png(p_hw, "goldfish_height_vs_width")
} else {
 warning("Not enough points for Height ~ Width (need >= 2).")
}

cat("\nDone. PNGs saved to:", path_out_dir, "\n")
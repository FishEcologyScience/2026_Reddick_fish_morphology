## --------------------------------------------------------------#
## Script name: script01-02_morphology_plots.R
##
## Purpose:
##   - Use pre-cleaned data (from script01-02_import_format.R) to
##     generate morphology plots and simple regressions per species,
##     plus multi-species plots for each plot class.
##
#### Repository layout (per template):
##   01_data/                       # inputs
##   01_data/processed/<species>/   # cleaned RDS per species
##   03_outputs/01_figures/         # figures
##   03_outputs/01_tables/          # tables (created if missing)
##
## Author: Marcus Rizzuto
## Date Created: 1/27/2026
##
## --------------------------------------------------------------#
##
## PLOT INVENTORY (what this script produces)
##
##   Per-species plots — stored in plots[[species]][[slot]]:
##   ┌─────┬──────────────────────────┬──────────────────────────────────────────┐
##   │  #  │ Slot name                │ What it shows                            │
##   ├─────┼──────────────────────────┼──────────────────────────────────────────┤
##   │  1  │ scatter_fl               │ Width ~ Fork Length (linear fit)         │
##   │  2  │ scatter_mass             │ Mass ~ Width (power-law, raw scale)      │
##   │  3  │ FL_by_mass               │ Mass ~ Fork Length (power-law, raw scale)│
##   │  4  │ loglog_width_by_FL       │ log(Width) ~ log(Fork Length)            │
##   │  5  │ powerlaw_width_by_mass   │ Width ~ Mass^b (log-log axes)            │
##   │  6  │ hist_FL                  │ Fork Length frequency histogram          │
##   └─────┴──────────────────────────┴──────────────────────────────────────────┘
##
##   Multi-species combined plots — stored in plots[["combined"]][[slot]]:
##     width_by_FL            All species: Width ~ Fork Length overlaid
##     width_by_mass_logx     All species: Width ~ Mass (log x-axis)
##     loglog_width_by_FL     All species: log(Width) ~ log(Fork Length)
##     powerlaw_width_by_mass All species: Width ~ Mass^b (log-log)
##     hist_FL_by_species     Fork Length histograms faceted by species
##
##   Summary tables created in the environment:
##     df_combined_summary    Per-species descriptive stats (n, mean, min, max)
##     df_combined_models     Power-law coefficients (a, b, R²) for Mass~Width
##     df_species_counts      Raw vs. NA/zero-filtered sample sizes per species
##
## --------------------------------------------------------------#


##### Guards & Setup ##########################################----
#-------------------------------------------------------------#
# Require cleaned objects from the import script:
# - combined_all : list of cleaned data frames by species
# - df_all       : single combined data frame of all species

if (!exists("combined_all")) {
 stop("ERROR: `combined_all` not found. Run the import/format script first (e.g., script01-01_import_format_singlefile.R).")
}
if (!exists("df_all")) {
 stop("ERROR: `df_all` not found. Run the import/format script first (e.g., script01-01_import_format_singlefile.R).")
}

# Output directories — figures are saved into a dated subfolder (YYYY-mm-dd)
# so each run's exports are isolated and won't overwrite previous versions.
path_figs_dir   <- file.path("03_outputs", "01_figures", format(Sys.Date(), "%Y-%m-%d"))
path_tables_dir <- file.path("03_outputs", "01_tables")
if (!dir.exists(path_figs_dir))   dir.create(path_figs_dir,   recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(path_tables_dir)) dir.create(path_tables_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize containers for results and plots
df_combined_summary <- tibble()
df_combined_models  <- tibble()  # collects coefficients for Width ~ log(Mass) model
plots <- list()                  # nested list of plots
plots[["combined"]] <- list()    # ensure multi-species plot container exists


# Helper to standardize plot captions per species
make_caption <- function(df, text, sp) paste0(sp, " (n = ", nrow(df), "): ", text)


##### QC Parameters #############################################----
#-------------------------------------------------------------#

# Minimum species-level sample size required to attempt regression fits.
# Effect: species below this threshold get a blank placeholder plot with a
#         note explaining the skip. Raise if you want stricter data requirements;
#         lower (minimum 3) if you need fits on very small samples.
MIN_N_PER_SPECIES <- 10

# Controls how regression trend lines are drawn on COMBINED (multi-species) plots.
# Options:
#   "overall"     - single global trend across all species (default; shows overall signal)
#   "per_species" - one trend line per species (clearer when species differ strongly)
#   "none"        - scatter points only, no regression lines
# Note: this setting does NOT affect per-species plots, only plots[["combined"]].
COMBINED_TREND_MODE <- "overall"  # change to "per_species" or "none" if needed

# ---- Histogram binning parameter (global) -------------------#
# Fixed bin width applied uniformly across all species histograms.
# Effect: smaller values show finer size-class detail but can look noisy on
#         small samples; larger values smooth the distribution.
BIN_WIDTH_CM <- 10   # 10 cm

# ---- Bar-spacing threshold (single reference line) ----------#
# A dashed reference line is drawn at this width value on Width plots,
# but ONLY if the data for that species actually reach this threshold.
# Change this if the biologically meaningful cutoff shifts (e.g., to 40 or 60 mm).
BAR_THRESHOLD_MM <- 50  # 5 cm = 50 mm

## ---- Global typography sizes for all plots ------------------#
# Adjust these to scale all text up or down uniformly across every plot.
TITLE_SIZE        <- 16  # plot title
SUBTITLE_SIZE     <- 13  # plot subtitle (if you use it)
AXIS_TITLE_SIZE   <- 14  # x/y axis titles
AXIS_TEXT_SIZE    <- 12  # x/y tick labels
LEGEND_TITLE_SIZE <- 12
LEGEND_TEXT_SIZE  <- 11
CAPTION_SIZE      <- 10

# Start from minimal but with a slightly larger base font
theme_set(ggplot2::theme_minimal(base_size = 12))

# Override key text elements globally
ggplot2::theme_update(
 plot.title    = ggplot2::element_text(size = TITLE_SIZE, face = "bold"),
 plot.subtitle = ggplot2::element_text(size = SUBTITLE_SIZE),
 plot.caption  = ggplot2::element_text(size = CAPTION_SIZE, colour = "grey30"),

 axis.title.x  = ggplot2::element_text(size = AXIS_TITLE_SIZE),
 axis.title.y  = ggplot2::element_text(size = AXIS_TITLE_SIZE),
 axis.text.x   = ggplot2::element_text(size = AXIS_TEXT_SIZE),
 axis.text.y   = ggplot2::element_text(size = AXIS_TEXT_SIZE),

 legend.title  = ggplot2::element_text(size = LEGEND_TITLE_SIZE),
 legend.text   = ggplot2::element_text(size = LEGEND_TEXT_SIZE)
)

### Core Data Processing
#----------------------------#
# The main per-species processing loop builds:
#   (1) Summary table row (with min/max)
#   (2) 5 plots per species (with guards for min N where applicable)
#   (3) QC flags (negative slope)


##### Species Loop (plots + simple models only) ###############----
#-------------------------------------------------------------#
# For each species:
#   - Compute summary statistics (counts, means, min/max)
#   - Build 5 plot types (with model annotations when data allow)
#   - Enforce minimum n thresholds for regressions (MIN_N_PER_SPECIES)
#   - Record model coefficients and QC flags
for (param_species in names(combined_all)) {
 
 cat("\n--- Plotting:", if (nzchar(param_species)) param_species else "UNKNOWN_SPECIES", "---\n")
 df_clean <- combined_all[[param_species]]
 
 # Skip unknown/blank/NA species names (do not plot/export them)
 sp_raw <- param_species
 sp_key <- tolower(trimws(sp_raw))
 if (is.na(sp_raw) || sp_key %in% c("", "unknown", "unknown_species")) next
 
 # Ensure per-species plot sublist exists before assigning plots
 if (is.null(plots[[param_species]]) || !is.list(plots[[param_species]])) {
  plots[[param_species]] <- list()
 }
 
 ### Species-level summary row (with min/max) ###
 #----------------------------------------------#
 # Summarize counts and central tendency + ranges for key morphometrics.
 df_summary <- df_clean %>%
  dplyr::summarise(
   species        = param_species,
   n_rows         = dplyr::n(),
   n_FL           = sum(!is.na(ForkLength_mm)),
   n_width        = sum(!is.na(Width_mm)),
   n_weight       = sum(!is.na(Mass_g)),
   FL_mean_mm     = mean(ForkLength_mm, na.rm = TRUE),
   FL_min_mm      = suppressWarnings(min(ForkLength_mm, na.rm = TRUE)),
   FL_max_mm      = suppressWarnings(max(ForkLength_mm, na.rm = TRUE)),
   width_mean_mm  = mean(Width_mm,      na.rm = TRUE),
   width_min_mm   = suppressWarnings(min(Width_mm,      na.rm = TRUE)),
   width_max_mm   = suppressWarnings(max(Width_mm,      na.rm = TRUE)),
   weight_mean_g  = mean(Mass_g,        na.rm = TRUE),
   weight_min_g   = suppressWarnings(min(Mass_g,        na.rm = TRUE)),
   weight_max_g   = suppressWarnings(max(Mass_g,        na.rm = TRUE))
  )
 df_combined_summary <- dplyr::bind_rows(df_combined_summary, df_summary)
 
 
 #### Plot 1: Width ~ Fork Length (linear fit) ####----
 # Scatter with optional linear regression; equation and R² annotated on plot.
 
 # Use only rows with both Width and ForkLength present
 loop_scatter_fl <- df_clean %>% dplyr::filter(!is.na(Width_mm), !is.na(ForkLength_mm))
 
 if (nrow(loop_scatter_fl) >= MIN_N_PER_SPECIES) {
  loop_p_scatter_fl <- ggplot(loop_scatter_fl, aes(ForkLength_mm, Width_mm)) +
   geom_point(color = "#2c7fb8", alpha = 0.6, size = 2) +
   labs(
    title   = paste0(param_species, " - Width by Fork Length"),
    x       = "Fork Length (mm)",
    y       = "Width (mm)",
    caption = make_caption(loop_scatter_fl, "Width vs fork length.", param_species)
   ) 
  
  # Fit requires at least 2 points
  if (nrow(loop_scatter_fl) >= 2) {
   loop_lm_fl    <- lm(Width_mm ~ ForkLength_mm, data = loop_scatter_fl)
   loop_coef_fl  <- coef(loop_lm_fl)
   loop_r2_fl    <- summary(loop_lm_fl)$r.squared
   loop_slope_fl <- unname(loop_coef_fl[["ForkLength_mm"]])  # slope
   loop_int_fl   <- unname(loop_coef_fl[["(Intercept)"]])    # intercept
   
   loop_eq_fl <- paste0(
    "Width = ", formatC(loop_slope_fl, format = "f", digits = 2),
    " · FL ", ifelse(loop_int_fl >= 0, "+ ", "- "),
    formatC(abs(loop_int_fl), format = "f", digits = 2), " (intercept)",
    "\nR^2 = ", formatC(loop_r2_fl, format = "f", digits = 3)
   )
   
   loop_p_scatter_fl <- loop_p_scatter_fl +
    geom_smooth(method = "lm", se = FALSE, color = "#1f78b4", linewidth = 0.9) +
    annotate("text",
             x = quantile(loop_scatter_fl$ForkLength_mm, 0.05, na.rm = TRUE),
             y = quantile(loop_scatter_fl$Width_mm,      0.95, na.rm = TRUE),
             label = loop_eq_fl, hjust = 0, vjust = 1, size = 3.5)
  }
  
  # ---- Dashed 50 mm line ONLY if the data reach >= 50 mm ----
  if (is.finite(max(loop_scatter_fl$Width_mm, na.rm = TRUE)) &&
      max(loop_scatter_fl$Width_mm, na.rm = TRUE) >= BAR_THRESHOLD_MM) {
   loop_p_scatter_fl <- loop_p_scatter_fl +
    geom_hline(yintercept = BAR_THRESHOLD_MM,
               linetype = "dashed", color = "firebrick",
               linewidth = 0.5, alpha = 0.8)
  }
 } else {
  # Placeholder when insufficient data for a meaningful plot
  loop_p_scatter_fl <- ggplot() + theme_void() +
   labs(caption = make_caption(loop_scatter_fl,
                               paste0("Insufficient data (n = ", nrow(loop_scatter_fl),
                                      " < ", MIN_N_PER_SPECIES, ") for Width~FL plot."), param_species))
 }
 plots[[param_species]][["scatter_fl"]] <- loop_p_scatter_fl
 
 #### Plot 2: Mass ~ Width (power-law, raw scale) ####----
 # Scatter with Mass on y, Width on x; fits W = a·X^b directly on the raw scale
 # using nls(), falling back to log-log with Duan's smearing correction if nls() fails.
 # Use this plot to predict mass from width in original units.
 # NOTE: Distinct from Plot 5 — Plot 2 uses Width as the predictor on raw axes;
 #       Plot 5 uses Mass as the predictor on log-log axes to visualize scaling.
 
 # Keep only positive, non-missing values
 loop_scatter_mass <- df_clean %>%
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(Mass_g),   Mass_g   > 0)
 
 # Prepare a model-row so we can append even if we don't fit
 model_row <- tibble(
  species = param_species,
  n_used  = nrow(loop_scatter_mass),
  a       = NA_real_,
  b       = NA_real_,
  r2      = NA_real_,
  method  = NA_character_
 )
 
 if (nrow(loop_scatter_mass) >= MIN_N_PER_SPECIES) {
  
  # Base scatter (Mass vs Width)
  loop_p_scatter_mass <- ggplot(loop_scatter_mass, aes(x = Width_mm, y = Mass_g)) +
   geom_point(color = "#f16913", alpha = 0.6, size = 2) +
   labs(
    title   = paste0(param_species, " - Mass by Width "),
    x       = "Width (mm)",
    y       = "Mass (g)",
    caption = make_caption(loop_scatter_mass, "Mass vs width with power-law fit; X = Width_mm.", param_species)
   )
  
  # ----- Fit W = a * X^b on the original scale (preferred) -----------
  # 1) Starting values from log-log regression
  lm_start <- lm(log(Mass_g) ~ log(Width_mm), data = loop_scatter_mass)
  b0       <- unname(coef(lm_start)[2])
  a0       <- exp(unname(coef(lm_start)[1]))
  
  # 2) Try nls(); if it fails, fall back to log–log with Duan smearing
  fit_info <- list(method = NA_character_, a = NA_real_, b = NA_real_, r2 = NA_real_)
  
  nls_fit <- try(
   nls(Mass_g ~ a * (Width_mm ^ b),
       data = loop_scatter_mass,
       start = list(a = a0, b = b0),
       control = nls.control(maxiter = 200, warnOnly = TRUE)),
   silent = TRUE
  )
  
  if (!inherits(nls_fit, "try-error")) {
   # nls success on raw scale
   a_hat <- unname(coef(nls_fit)["a"])
   b_hat <- unname(coef(nls_fit)["b"])
   y     <- loop_scatter_mass$Mass_g
   yhat  <- fitted(nls_fit)
   r2_ml <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
   
   fit_info <- list(method = "nls_raw", a = a_hat, b = b_hat, r2 = r2_ml)
   
  } else {
   # Fallback: log–log with Duan’s smearing correction
   lm_fit <- lm(log(Mass_g) ~ log(Width_mm), data = loop_scatter_mass)
   b_hat  <- unname(coef(lm_fit)[2])
   a_raw  <- exp(unname(coef(lm_fit)[1]))
   smear  <- mean(exp(residuals(lm_fit)))   # Duan smearing factor
   a_hat  <- a_raw * smear
   
   y     <- loop_scatter_mass$Mass_g
   yhat  <- a_hat * (loop_scatter_mass$Width_mm ^ b_hat)
   r2_ml <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
   
   fit_info <- list(method = "loglm_smear", a = a_hat, b = b_hat, r2 = r2_ml)
  }
  
  # ----- Prediction curve on original axes ---------------------------
  x_rng_w <- range(loop_scatter_mass$Width_mm, na.rm = TRUE)
  x_seq_w <- seq(x_rng_w[1], x_rng_w[2], length.out = 200)
  
  pred_curve <- tibble(Width_mm = x_seq_w) |>
   dplyr::mutate(Mass_pred = fit_info$a * (Width_mm ^ fit_info$b))
  
  loop_p_scatter_mass <- loop_p_scatter_mass +
   geom_line(data = pred_curve, aes(x = Width_mm, y = Mass_pred),
             color = "#cc4c02", linewidth = 1)
  
  # ----- Equation annotation (use generic X to avoid implying girth) --
  a_lbl  <- formatC(fit_info$a, format = "e", digits = 2)  # 'a' often small; scientific notation
  b_lbl  <- formatC(fit_info$b, format = "f", digits = 3)
  r2_lbl <- if (is.finite(fit_info$r2)) paste0("\nR^2 = ", formatC(fit_info$r2, format = "f", digits = 3)) else ""
  
  eq_label <- paste0("Mass = ", a_lbl, " · Width^", b_lbl, r2_lbl)
  
  loop_p_scatter_mass <- loop_p_scatter_mass +
   annotate("text",
            x = quantile(loop_scatter_mass$Width_mm, 0.05, na.rm = TRUE),
            y = quantile(loop_scatter_mass$Mass_g,   0.95, na.rm = TRUE),
            label = eq_label,
            hjust = 0, vjust = 1, size = 3.5)
  
  # Optional: keep your 50 mm reference line behavior
  if (is.finite(max(loop_scatter_mass$Width_mm, na.rm = TRUE)) &&
      max(loop_scatter_mass$Width_mm, na.rm = TRUE) >= BAR_THRESHOLD_MM) {
   loop_p_scatter_mass <- loop_p_scatter_mass +
    geom_vline(xintercept = BAR_THRESHOLD_MM,
               linetype = "dashed", color = "firebrick",
               linewidth = 0.5, alpha = 0.8)
  }
  
  # Record fitted parameters
  model_row <- tibble(
   species = param_species,
   n_used  = nrow(loop_scatter_mass),
   a       = fit_info$a,
   b       = fit_info$b,
   r2      = fit_info$r2,
   method  = fit_info$method
  )
  
 } else {
  loop_p_scatter_mass <- ggplot() + theme_void() +
   labs(caption = make_caption(loop_scatter_mass,
                               paste0("Insufficient data (n = ", nrow(loop_scatter_mass),
                                      " < ", MIN_N_PER_SPECIES, ") for Mass~Width (power-law)."),
                               param_species))
 }
 
 # Store the plot and the model row (same destinations as before)
 plots[[param_species]][["scatter_mass"]] <- loop_p_scatter_mass
 df_combined_models <- dplyr::bind_rows(df_combined_models, model_row)
 
 #### Plot 3: Mass ~ Fork Length (power-law, raw scale) ####----
 # Scatter with Mass on y, Fork Length on x; fits W = a·L^b on the raw scale.
 # Standard fisheries length-weight relationship plot.
 
 # Use only rows with positive W and L
 loop_fl_mass <- df_clean %>%
  dplyr::filter(!is.na(ForkLength_mm), ForkLength_mm > 0,
                !is.na(Mass_g),       Mass_g       > 0)
 
 if (nrow(loop_fl_mass) >= MIN_N_PER_SPECIES) {
  
  # Base scatter
  loop_p_fl_mass <- ggplot(loop_fl_mass, aes(x = ForkLength_mm, y = Mass_g)) +
   geom_point(color = "#d95f02", alpha = 0.6, size = 2) +
   labs(
    title   = paste0(param_species, " - Mass by Fork Length "),
    x       = "Fork Length (mm)",
    y       = "Mass (g)",
    caption = make_caption(loop_fl_mass, "Mass vs fork length with power-law fit (W = a·L^b).", param_species)
   )
  
  # ----- Fit W = a * L^b ----------------------------------------------
  # 1) Get robust starting values from log-log linear fit
  lm_start   <- lm(log(Mass_g) ~ log(ForkLength_mm), data = loop_fl_mass)
  b0         <- unname(coef(lm_start)[2])
  a0         <- exp(unname(coef(lm_start)[1]))
  
  # 2) Try nls() on the original scale; if it fails, fall back to log-log with smearing
  fit_info <- list(method = NA_character_, a = NA_real_, b = NA_real_, r2 = NA_real_)
  
  nls_ok <- TRUE
  nls_fit <- try(
   nls(Mass_g ~ a * (ForkLength_mm ^ b),
       data = loop_fl_mass,
       start = list(a = a0, b = b0),
       control = nls.control(maxiter = 200, warnOnly = TRUE)),
   silent = TRUE
  )
  
  if (inherits(nls_fit, "try-error")) nls_ok <- FALSE
  
  if (nls_ok) {
   # Extract coefficients from nls
   coef_nls <- coef(nls_fit)
   a_hat    <- unname(coef_nls["a"])
   b_hat    <- unname(coef_nls["b"])
   
   # Raw-scale R^2 (pseudo-R^2)
   y     <- loop_fl_mass$Mass_g
   yhat  <- fitted(nls_fit)
   r2_ml <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
   
   fit_info <- list(method = "nls_raw", a = a_hat, b = b_hat, r2 = r2_ml)
   
  } else {
   # ----- Fallback: log-log fit with Duan's smearing correction -----
   # log(W) = log(a) + b log(L)  =>  W = a * L^b, with a corrected by smearing factor
   lm_fit <- lm(log(Mass_g) ~ log(ForkLength_mm), data = loop_fl_mass)
   b_hat  <- unname(coef(lm_fit)[2])
   a_raw  <- exp(unname(coef(lm_fit)[1]))
   smear  <- mean(exp(residuals(lm_fit)))   # Duan smearing
   a_hat  <- a_raw * smear
   
   # R^2 on raw scale using back-transformed fitted values
   y     <- loop_fl_mass$Mass_g
   yhat  <- a_hat * (loop_fl_mass$ForkLength_mm ^ b_hat)
   r2_ml <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
   
   fit_info <- list(method = "loglm_smear", a = a_hat, b = b_hat, r2 = r2_ml)
  }
  
  # ----- Prediction curve on original axes ---------------------------
  x_rng_fl <- range(loop_fl_mass$ForkLength_mm, na.rm = TRUE)
  x_seq_fl <- seq(x_rng_fl[1], x_rng_fl[2], length.out = 200)
  
  pred_curve <- tibble(ForkLength_mm = x_seq_fl) |>
   dplyr::mutate(Mass_pred = fit_info$a * (ForkLength_mm ^ fit_info$b))
  
  loop_p_fl_mass <- loop_p_fl_mass +
   geom_line(data = pred_curve, aes(x = ForkLength_mm, y = Mass_pred),
             color = "#bf5b17", linewidth = 1)
  
  # ----- Equation annotation: W = a · L^b ----------------------------
  # 'a' can be small; show in engineering/scientific format
  a_lbl <- formatC(fit_info$a, format = "e", digits = 2)
  b_lbl <- formatC(fit_info$b, format = "f", digits = 3)
  r2_lbl <- if (is.finite(fit_info$r2)) paste0("\nR^2 = ", formatC(fit_info$r2, format = "f", digits = 3)) else ""
  
  eq_label <- paste0("Mass = ", a_lbl, " · FL^", b_lbl, r2_lbl)
  
  loop_p_fl_mass <- loop_p_fl_mass +
   annotate(
    "text",
    x = quantile(loop_fl_mass$ForkLength_mm, 0.05, na.rm = TRUE),
    y = quantile(loop_fl_mass$Mass_g,        0.95, na.rm = TRUE),
    label = eq_label,
    hjust = 0, vjust = 1, size = 3.6
   )
  
 } else {
  loop_p_fl_mass <- ggplot() + theme_void() +
   labs(
    caption = make_caption(
     loop_fl_mass,
     paste0("Insufficient data (n = ", nrow(loop_fl_mass),
            " < ", MIN_N_PER_SPECIES, ") for Mass~Fork Length (power-law) plot."),
     param_species
    )
   )
 }
 
 # Store the plot (same slot name as before)
 plots[[param_species]][["FL_by_mass"]] <- loop_p_fl_mass

 #### Plot 4: log(Width) ~ log(Fork Length) (log-log axes) ####----
 # Both axes on log10 scale; linearizes the power-law relationship.
 # Good for checking whether Width and Fork Length scale allometrically.
 
 # Use positive values only for log transformation
 loop_loglog_fl <- df_clean %>%
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(ForkLength_mm), ForkLength_mm > 0)
 
 if (nrow(loop_loglog_fl) >= MIN_N_PER_SPECIES) {
  loop_p_loglog_fl <- ggplot(loop_loglog_fl, aes(ForkLength_mm, Width_mm)) +
   geom_point(color = "#4daf4a", alpha = 0.6, size = 2) +
   scale_x_log10(breaks = scales::log_breaks(n = 10), labels = scales::label_comma()) +
   scale_y_log10(breaks = scales::log_breaks(n = 10), labels = scales::label_comma()) +
   annotation_logticks(sides = "bl") +
   labs(
    title   = paste0(param_species, " - log(Width) ~ log(Fork Length)"),
    x       = "Fork Length (mm, log10 scale)",
    y       = "Width (mm, log10 scale)",
    caption = make_caption(loop_loglog_fl, "log(Width) vs log(ForkLength).", param_species)
   ) 
  
  if (nrow(loop_loglog_fl) >= 3) {
   loop_fit_ll   <- lm(log(Width_mm) ~ log(ForkLength_mm), data = loop_loglog_fl)
   loop_coefs_ll <- coef(loop_fit_ll)
   loop_r2_ll    <- summary(loop_fit_ll)$r.squared
   loop_slope_ll <- unname(loop_coefs_ll[["log(ForkLength_mm)"]])  # slope in log space
   loop_int_ll   <- unname(loop_coefs_ll[["(Intercept)"]])
   
   # Predict on original scale across observed ForkLength range
   loop_rng_x_fl <- range(loop_loglog_fl$ForkLength_mm, na.rm = TRUE)
   loop_xseq_fl  <- seq(loop_rng_x_fl[1], loop_rng_x_fl[2], length.out = 200)
   loop_pred_ll  <- tibble(ForkLength_mm = loop_xseq_fl) %>%
    dplyr::mutate(Width_pred = exp(loop_int_ll + loop_slope_ll * log(ForkLength_mm)))
   
   loop_p_loglog_fl <- loop_p_loglog_fl +
    geom_line(data = loop_pred_ll, aes(ForkLength_mm, Width_pred),
              color = "#238b45", linewidth = 1) +
    annotate("text",
             x = quantile(loop_loglog_fl$ForkLength_mm, 0.05, na.rm = TRUE),
             y = quantile(loop_loglog_fl$Width_mm,       0.95, na.rm = TRUE),
             label = paste0("log(Width) = slope·log(FL) ", ifelse(loop_int_ll >= 0, "+ ", "- "),
                            formatC(abs(loop_int_ll), format = "f", digits = 2),
                            "\nslope = ", formatC(loop_slope_ll, format = "f", digits = 3),
                            "  R^2 = ", formatC(loop_r2_ll, format = "f", digits = 3)),
             hjust = 0, vjust = 1, size = 3.5)
  }
 } else {
  loop_p_loglog_fl <- ggplot() + theme_void() +
   labs(caption = make_caption(loop_loglog_fl,
                               paste0("Insufficient data (n = ", nrow(loop_loglog_fl),
                                      " < ", MIN_N_PER_SPECIES, ") for log-log Width~FL."),
                               param_species))
 }
 plots[[param_species]][["loglog_width_by_FL"]] <- loop_p_loglog_fl
 
 
 #### Plot 5: Width ~ Mass^b (power-law, log-log axes) ####----
 # Both axes on log10 scale; Mass is the predictor, Width is the response.
 # Fit via log-log linear regression: log(Width) = b·log(Mass) + log(a).
 # NOTE: Distinct from Plot 2 — here Mass drives Width (reverse direction),
 #       and the fit is on log-log axes rather than the raw scale.
 #       Use this plot to visualize the allometric scaling exponent b.
 
 # Use positive values only for log transformation
 loop_power_mass <- df_clean %>%
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(Mass_g),   Mass_g > 0)
 
 if (nrow(loop_power_mass) >= MIN_N_PER_SPECIES) {
  loop_p_power_mass <- ggplot(loop_power_mass, aes(Mass_g, Width_mm)) +
   geom_point(color = "#984ea3", alpha = 0.6, size = 2) +
   scale_x_log10(breaks = scales::log_breaks(n = 10), labels = scales::label_comma()) +
   scale_y_log10(breaks = scales::log_breaks(n = 10), labels = scales::label_comma()) +
   annotation_logticks(sides = "bl") +
   labs(
    title   = paste0(param_species, " - Power law: Width ~ Mass^b"),
    x       = "Mass (g, log10 scale)",
    y       = "Width (mm, log10 scale)",
    caption = make_caption(loop_power_mass, "Power-law width vs mass.", param_species)
   ) 
  
  if (nrow(loop_power_mass) >= 3) {
   loop_fit_pw   <- lm(log(Width_mm) ~ log(Mass_g), data = loop_power_mass)
   loop_coefs_pw <- coef(loop_fit_pw)
   loop_r2_pw    <- summary(loop_fit_pw)$r.squared
   loop_b_pw     <- unname(loop_coefs_pw[["log(Mass_g)"]])  # exponent b
   loop_log_a_pw <- unname(loop_coefs_pw[["(Intercept)"]])
   loop_a_pw     <- exp(loop_log_a_pw)                      # coefficient a
   
   # Predict on original scale across observed Mass range
   loop_rng_x_mass <- range(loop_power_mass$Mass_g, na.rm = TRUE)
   loop_xseq_mass2 <- seq(loop_rng_x_mass[1], loop_rng_x_mass[2], length.out = 200)
   loop_pred_pw    <- tibble(Mass_g = loop_xseq_mass2) %>%
    dplyr::mutate(Width_pred = loop_a_pw * (Mass_g ^ loop_b_pw))
   
   loop_p_power_mass <- loop_p_power_mass +
    geom_line(data = loop_pred_pw, aes(Mass_g, Width_pred),
              color = "#6a51a3", linewidth = 1) +
    annotate("text",
             x = quantile(loop_power_mass$Mass_g, 0.05, na.rm = TRUE),
             y = quantile(loop_power_mass$Width_mm, 0.95, na.rm = TRUE),
             label = paste0("Width = ", formatC(loop_a_pw, format = "f", digits = 2),
                            " · Mass^", formatC(loop_b_pw, format = "f", digits = 3),
                            "\nR^2 = ", formatC(loop_r2_pw, format = "f", digits = 3)),
             hjust = 0, vjust = 1, size = 3.5)
  }
  
  # ---- Dashed 50 mm line ONLY if the data reach >= 50 mm (log scale OK) ----
  if (is.finite(max(loop_power_mass$Width_mm, na.rm = TRUE)) &&
      max(loop_power_mass$Width_mm, na.rm = TRUE) >= BAR_THRESHOLD_MM) {
   loop_p_power_mass <- loop_p_power_mass +
    geom_hline(yintercept = BAR_THRESHOLD_MM,
               linetype = "dashed", color = "firebrick",
               linewidth = 0.5, alpha = 0.8)
  }
  
 } else {
  loop_p_power_mass <- ggplot() + theme_void() +
   labs(caption = make_caption(loop_power_mass,
                               paste0("Insufficient data (n = ", nrow(loop_power_mass),
                                      " < ", MIN_N_PER_SPECIES, ") for power-law Width~Mass."),
                               param_species))
 }
 plots[[param_species]][["powerlaw_width_by_mass"]] <- loop_p_power_mass
 
 
 #### Plot 6: Fork Length histogram ####----
 # Fixed bin width (BIN_WIDTH_CM) across all species for direct comparison.
 # Red solid line = mean; green dashed line = median.
 
 # Data for histogram
 loop_hist_df <- df_clean %>%
  dplyr::filter(!is.na(ForkLength_mm))
 
 # Only build if we have *some* data points (>= 1).
 loop_hist_caption <- NA_character_  # populated below; used in patchwork panel caption
 if (nrow(loop_hist_df) >= 1) {

  # Mean & median markers
  loop_fl_mean   <- mean(loop_hist_df$ForkLength_mm, na.rm = TRUE)
  loop_fl_median <- stats::median(loop_hist_df$ForkLength_mm, na.rm = TRUE)

  # Caption string for the patchwork panel footer
  loop_mean_mass <- mean(df_clean$Mass_g, na.rm = TRUE)
  loop_hist_caption <- paste0(
   "n = ", nrow(loop_hist_df), "  |  ",
   "Mean fork length: ", formatC(loop_fl_mean,   digits = 1, format = "f"), " mm  |  ",
   "Mean mass: ",        formatC(loop_mean_mass, digits = 1, format = "f"), " g"
  )
  
  p_hist_fl <- ggplot(loop_hist_df, aes(x = ForkLength_mm)) +
   geom_histogram(
    binwidth = BIN_WIDTH_CM,        # <--- fixed width (10 cm)
    boundary = 0,
    color = "grey30",
    fill = "#74a9cf",
    alpha = 0.8
   ) +
   annotate("segment",
            x = loop_fl_mean, xend = loop_fl_mean,
            y = -Inf, yend = Inf,
            colour = "#de2d26", linetype = "solid", linewidth = 1.0, lineend = "butt") +
   annotate("segment",
            x = loop_fl_median, xend = loop_fl_median,
            y = -Inf, yend = Inf,
            colour = "#238b45", linetype = "dashed", linewidth = 1.0, lineend = "butt") +
   labs(
    title   = paste0(param_species, " - Histogram of Fork Length"),
    x       = "Fork Length (mm)",
    y       = "Count",
    caption = make_caption(
     loop_hist_df,
     paste0(
      "Histogram of fork lengths. Fixed bin width = ",
      formatC(BIN_WIDTH_CM, digits = 0, format = "f"), " cm. ",
      "Red = mean (", formatC(loop_fl_mean, digits = 1, format = "f"),
      " mm), green dashed = median (", formatC(loop_fl_median, digits = 1, format = "f"), " mm)."
     ),
     param_species
    )
   ) +
 scale_x_continuous(
  breaks = scales::breaks_width(50),
  minor_breaks = scales::breaks_width(25),
  expand = c(0.01, 0)
 )
  
  plots[[param_species]][["hist_FL"]] <- p_hist_fl
 } else {
  plots[[param_species]][["hist_FL"]] <- ggplot() + theme_void() +
   labs(caption = make_caption(loop_hist_df, "No fork length data available.", param_species))
 }
 #### Patchwork panel ####----
 # Assembled here while all 6 plot objects are in scope.
 # Only built if all regression plots had sufficient data (nrow >= MIN_N_PER_SPECIES).
 # Layout: Histogram | Mass~FL  /  Width~FL | log(Width)~log(FL)  /  Mass~Width | Width~Mass^b
 if (nrow(dplyr::filter(df_clean, !is.na(ForkLength_mm), !is.na(Width_mm)))     >= MIN_N_PER_SPECIES &&
     nrow(dplyr::filter(df_clean, !is.na(Mass_g), Mass_g > 0, !is.na(Width_mm))) >= MIN_N_PER_SPECIES &&
     nrow(dplyr::filter(df_clean, !is.na(ForkLength_mm), ForkLength_mm > 0, !is.na(Mass_g), Mass_g > 0)) >= MIN_N_PER_SPECIES) {

  # Step 1: assemble layout and apply sub-plot theme adjustments.
  # Two-step assignment ensures & theme() only affects sub-plots, not the panel title.
  loop_pw <-
   (plots[[param_species]][["hist_FL"]]      | plots[[param_species]][["FL_by_mass"]])          /
   (plots[[param_species]][["scatter_fl"]]   | plots[[param_species]][["loglog_width_by_FL"]])  /
   (plots[[param_species]][["scatter_mass"]] | plots[[param_species]][["powerlaw_width_by_mass"]]) &
   ggplot2::theme(
    plot.title   = ggplot2::element_text(size = 9, face = "plain"),  # short sub-plot titles
    plot.caption = ggplot2::element_blank()                          # redundant at panel level
   )

  # Step 2: add panel-level title, letter tags (a)-(f), and histogram caption footer.
  plots[[param_species]][["patchwork"]] <- loop_pw +
   patchwork::plot_annotation(
    title      = param_species,
    caption    = if (!is.na(loop_hist_caption)) loop_hist_caption else NULL,
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
    theme      = ggplot2::theme(
     plot.title   = ggplot2::element_text(size = TITLE_SIZE + 2, face = "bold"),
     plot.caption = ggplot2::element_text(size = CAPTION_SIZE, colour = "grey30")
    )
   )

 } else {
  cat("Skipping patchwork for", param_species, "— insufficient data for one or more regression plots.\n")
 }

 cat("Finished:", if (nzchar(param_species)) param_species else "UNKNOWN_SPECIES", "\n")
} # end species loop


##### Multi-species plots #####################################----
#-------------------------------------------------------------#
# Build cross-species views using df_all (all species combined).
# Trend line behaviour on these plots is controlled by COMBINED_TREND_MODE
# (set in QC Parameters above): "overall", "per_species", or "none".

#### Combined 1: Width ~ Fork Length (linear) ####----
# Shows whether species share similar width-to-length slopes or diverge,
# suggesting morphological differences relevant to gear selectivity or condition.
df_combined_fl <- df_all %>% dplyr::filter(!is.na(Width_mm), !is.na(ForkLength_mm))
plots[["combined"]][["width_by_FL"]] <- ggplot(df_combined_fl, aes(ForkLength_mm, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 geom_smooth(method = "lm", se = FALSE) +
 labs(title = "All species - Width by Fork Length (linear)",
      x = "Fork Length (mm)", y = "Width (mm)") +
 theme(legend.position = "bottom")

#### Combined 2: Width ~ Mass (log x-axis) ####----
# Log-transforms mass to linearize the relationship; useful for comparing
# how body width tracks mass gain across species of very different sizes.
df_combined_mass_logx <- df_all %>% dplyr::filter(!is.na(Width_mm), !is.na(Mass_g), Mass_g > 0)
plots[["combined"]][["width_by_mass_logx"]] <- ggplot(df_combined_mass_logx, aes(Mass_g, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 stat_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
 labs(title = "All species - Width by Mass (log x)",
      x = "Mass (g)", y = "Width (mm)") +
 theme(legend.position = "bottom")

#### Combined 3: log(Width) ~ log(Fork Length) (log-log) ####----
# Log-log plot linearizes the power-law; parallel lines = same scaling exponent,
# offset lines = same shape but different size at a given length.
# Trend controlled by COMBINED_TREND_MODE.
df_combined_loglog_fl <- df_all %>%
 dplyr::filter(!is.na(Width_mm), Width_mm > 0,
               !is.na(ForkLength_mm), ForkLength_mm > 0)

plots[["combined"]][["loglog_width_by_FL"]] <-
 ggplot(df_combined_loglog_fl, aes(ForkLength_mm, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 scale_x_log10() + scale_y_log10() +
 labs(title = "All species - log(Width) ~ log(Fork Length)",
      x = "Fork Length (mm, log10 scale)", y = "Width (mm, log10 scale)") +
 theme(legend.position = "bottom")

if (COMBINED_TREND_MODE == "overall") {
 # One global trend line across all species
 plots[["combined"]][["loglog_width_by_FL"]] <- plots[["combined"]][["loglog_width_by_FL"]] +
  stat_smooth(method = "lm",
              formula = y ~ x, se = FALSE,
              mapping = aes(x = log(ForkLength_mm), y = log(Width_mm)),
              inherit.aes = FALSE, color = "black")
} else if (COMBINED_TREND_MODE == "per_species") {
 # One trend line per species (clearer when species differ strongly)
 plots[["combined"]][["loglog_width_by_FL"]] <- plots[["combined"]][["loglog_width_by_FL"]] +
  geom_smooth(method = "lm", se = FALSE)
} # else "none": scatter points only

#### Combined 4: Width ~ Mass^b (power-law, log-log) ####----
# Visualizes allometric scaling of width with mass across species.
# Steeper slopes indicate width grows faster relative to mass.
# Trend controlled by COMBINED_TREND_MODE.
df_combined_power_mass <- df_all %>%
 dplyr::filter(!is.na(Width_mm), Width_mm > 0,
               !is.na(Mass_g),   Mass_g > 0)

plots[["combined"]][["powerlaw_width_by_mass"]] <-
 ggplot(df_combined_power_mass, aes(Mass_g, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 scale_x_log10() + scale_y_log10() +
 labs(title = "All species - Power law: Width ~ Mass^b",
      x = "Mass (g, log10 scale)", y = "Width (mm, log10 scale)") +
 theme(legend.position = "bottom")

if (COMBINED_TREND_MODE == "overall") {
 plots[["combined"]][["powerlaw_width_by_mass"]] <- plots[["combined"]][["powerlaw_width_by_mass"]] +
  stat_smooth(method = "lm",
              formula = y ~ x, se = FALSE,
              mapping = aes(x = log(Mass_g), y = log(Width_mm)),
              inherit.aes = FALSE, color = "black")
} else if (COMBINED_TREND_MODE == "per_species") {
 plots[["combined"]][["powerlaw_width_by_mass"]] <- plots[["combined"]][["powerlaw_width_by_mass"]] +
  geom_smooth(method = "lm", se = FALSE)
} # else "none": no trend lines

#### Combined 5: Fork Length histograms faceted by species ####----
# Side-by-side size structure comparison; y-axes are free (scales = "free_y")
# so rare species are still legible alongside abundant ones.
df_all_hist <- df_all %>% dplyr::filter(!is.na(ForkLength_mm))
plots[["combined"]][["hist_FL_by_species"]] <- ggplot(df_all_hist, aes(x = ForkLength_mm)) +
 geom_histogram(color = "grey30", fill = "#9ecae1", alpha = 0.85,
                binwidth = BIN_WIDTH_CM, boundary = 0) +
 labs(
  title = "Histogram of Fork Lengths by Species",
  x = "Fork Length (mm)",
  y = "Count"
 ) +
 facet_wrap(~ Species, scales = "free_y") +
 theme_minimal() +
 theme(legend.position = "none")

# Send all plots to the RStudio Plots/Viewer pane.
# To view a single plot without printing all, use e.g.:
#   plots[["Goldfish"]][["scatter_fl"]]
#   plots[["combined"]][["hist_FL_by_species"]]
#   print(plots)
for (sp in names(combined_all)) print(plots[[sp]][["patchwork"]])

##### Species-level counts (raw vs filtered) ##################----
# Helps diagnose if species are genuinely low-sample vs heavily filtered by NA/zero.
df_species_counts <- df_all %>%
 dplyr::group_by(Species) %>%
 dplyr::summarise(
  n_all           = dplyr::n(),
  n_FL_Width      = sum(!is.na(ForkLength_mm) & !is.na(Width_mm)),
  n_Mass_Width_pos= sum(!is.na(Mass_g) & Mass_g > 0 & !is.na(Width_mm)),
  .groups = "drop"
 ) %>%
 dplyr::arrange(n_all)


##### Optional exports (kept commented) #######################----
#-------------------------------------------------------------#
# readr::write_csv(df_combined_summary, file.path(path_tables_dir, "_combined_summary_by_species.csv"))
# readr::write_csv(df_combined_models,  file.path(path_tables_dir, "_combined_log_model_coefficients.csv"))
# ggsave(file.path(path_figs_dir, "combined_width_by_FL.png"),        plots[["combined"]][["width_by_FL"]],        width = 7, height = 5, dpi = 300)
# ggsave(file.path(path_figs_dir, "combined_width_by_mass_logx.png"), plots[["combined"]][["width_by_mass_logx"]], width = 7, height = 5, dpi = 300)
# ggsave(file.path(path_figs_dir, "combined_loglog_width_by_FL.png"), plots[["combined"]][["loglog_width_by_FL"]], width = 7, height = 5, dpi = 300)
# ggsave(file.path(path_figs_dir, "combined_powerlaw_width_by_mass.png"), plots[["combined"]][["powerlaw_width_by_mass"]], width = 7, height = 5, dpi = 300)
# ggsave(file.path(path_figs_dir, "combined_hist_FL_by_species.png"),
    # plots[["combined"]][["hist_FL_by_species"]], width = 9, height = 7, dpi = 300)
# readr::write_csv(df_species_counts,  file.path(path_tables_dir, "_species_counts_raw_vs_filtered.csv"))
# readr::write_csv(df_combined_summary, file.path(path_tables_dir, "_combined_summary_by_species.csv"))  # already present
# readr::write_csv(df_combined_models,  file.path(path_tables_dir, "_combined_log_model_coefficients.csv"))  # already present

 # Per-species plot export (compact, all individual plots)
# for (sp in names(combined_all)) for (nm in names(plots[[sp]]))
# ggsave(file.path(path_figs_dir, paste0(gsub("[^A-Za-z0-9_\\-]", "_", sp), "_", nm, ".png")),
# plots[[sp]][[nm]], width = 7, height = 5, dpi = 300)

# Per-species patchwork panel export — letter size (8.5x11), all species at once
for (sp in names(combined_all))
 ggsave(file.path(path_figs_dir, paste0(gsub("[^A-Za-z0-9_\\-]", "_", sp), "_patchwork.png")),
        plots[[sp]][["patchwork"]], width = 8.5, height = 11, dpi = 300)

# Single-species patchwork export (swap species name as needed)
# ggsave(file.path(path_figs_dir, "Goldfish_patchwork.png"), plots[["Goldfish"]][["patchwork"]], width = 8.5, height = 11, dpi = 300)
# ggsave(file.path(path_figs_dir, "Rudd_patchwork.png"),     plots[["Rudd"]][["patchwork"]],     width = 8.5, height = 11, dpi = 300)

##### Cleanup (remove temporary loop/temp objects) #############----
#-------------------------------------------------------------#
# Remove transient objects created inside the loop or as temp_* helpers.
loop_temp_vars <- ls(pattern = "^(loop_|temp_)")
if (length(loop_temp_vars)) {
 rm(list = loop_temp_vars)
}

cat("\nAll plots complete.\n")
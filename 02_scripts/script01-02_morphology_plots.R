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

# Output directories (created if missing; exports remain commented out below)
path_figs_dir   <- file.path("03_outputs", "01_figures")
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
# Minimum species-level sample size required to draw/fits
MIN_N_PER_SPECIES <- 10

# Controls how combined regression lines are drawn:
#   - "overall"      : single global trend (original behavior for log-log/power)
#   - "per_species"  : one trend line per species (usually clearer)
#   - "none"         : no regression lines on combined plots
COMBINED_TREND_MODE <- "overall"  # change to "per_species" or "none" if needed

# ---- Histogram binning parameter (global) -------------------#
# single fixed width across all lengths.
BIN_WIDTH_CM <- 10   # 10 cm 

## ---- Global typography sizes for all plots ------------------#
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

# ---- Bar-spacing threshold (single reference line) ----------#
BAR_THRESHOLD_MM <- 50  # 5 cm = 50 mm

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
 
 
 #### Plot 1: Width ~ Fork Length (linear) ------------------#
 ### Minor: Scatter + optional linear fit with equation ###
 
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
    "y = ", formatC(loop_slope_fl, format = "f", digits = 2),
    "x", ifelse(loop_int_fl >= 0, " + ", " - "),
    formatC(abs(loop_int_fl), format = "f", digits = 2),
    "\nR^2 = ", formatC(loop_r2_fl, format = "f", digits = 4)
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
 
 #### Plot 2: Width ~ log(Mass) (Mass on Y) -------------------#
 ### Scatter with Mass on y, model still Width ~ ln(Mass); inverse curve + fitted-form equation ###
 
 # Use only rows with Width and positive Mass
 loop_scatter_mass <- df_clean %>%
  dplyr::filter(!is.na(Width_mm), !is.na(Mass_g), Mass_g > 0)
 
 # Initialize model row (remains NA if too few points)
 temp_model_row <- tibble(
  species = param_species,
  n_used  = nrow(loop_scatter_mass),
  alpha   = NA_real_,   # slope for ln(Mass)
  beta    = NA_real_,   # intercept
  r2      = NA_real_,
  flag_negative_slope = NA  # TRUE if alpha < 0
 )
 
 if (nrow(loop_scatter_mass) >= MIN_N_PER_SPECIES) {
  # MASS on Y, WIDTH on X
  loop_p_scatter_mass <- ggplot(loop_scatter_mass, aes(x = Width_mm, y = Mass_g)) +
   geom_point(color = "#f16913", alpha = 0.6, size = 2) +
   labs(
    title   = paste0(param_species, " - Mass by Width"),
    x       = "Width (mm)",
    y       = "Mass (g)",
    caption = make_caption(loop_scatter_mass, "Mass vs width (model on ln(Mass)).", param_species)
   ) 
  
  if (nrow(loop_scatter_mass) >= 3) {
   # Fit original model (predicting Width): Width = alpha * ln(Mass) + beta
   loop_log_fit_mass <- lm(Width_mm ~ log(Mass_g), data = loop_scatter_mass)
   loop_coefs_mass   <- coef(loop_log_fit_mass)
   loop_alpha_mass   <- unname(loop_coefs_mass["log(Mass_g)"])     # slope alpha
   loop_beta_mass    <- unname(loop_coefs_mass["(Intercept)"])     # intercept beta
   loop_r2_mass      <- summary(loop_log_fit_mass)$r.squared
   
   temp_model_row <- tibble(
    species = param_species,
    n_used  = nrow(loop_scatter_mass),
    alpha   = loop_alpha_mass,
    beta    = loop_beta_mass,
    r2      = loop_r2_mass,
    flag_negative_slope = is.finite(loop_alpha_mass) && loop_alpha_mass < 0
   )
   
   # Invert model to draw curve on flipped axes:
   # Width = alpha*ln(Mass) + beta  ==>  Mass = exp((Width - beta)/alpha)
   x_rng_w  <- range(loop_scatter_mass$Width_mm, na.rm = TRUE)
   x_seq_w  <- seq(x_rng_w[1], x_rng_w[2], length.out = 200)
   pred_inv <- tibble(Width_mm = x_seq_w) %>%
    dplyr::mutate(Mass_pred = exp((Width_mm - loop_beta_mass) / loop_alpha_mass))
   
   loop_p_scatter_mass <- loop_p_scatter_mass +
    geom_line(data = pred_inv, aes(x = Width_mm, y = Mass_pred), color = "#cc4c02", linewidth = 1)
   
# ---- Equation annotation  ----
   loop_eq_mass <- paste0(
    "y = ", formatC(loop_alpha_mass, format = "f", digits = 2), " ln(Mass)",
    ifelse(loop_beta_mass >= 0, " + ", " - "),
    formatC(abs(loop_beta_mass), format = "f", digits = 2),
    "\nR^2 = ", formatC(loop_r2_mass, format = "f", digits = 4)
   )
   
   loop_p_scatter_mass <- loop_p_scatter_mass +
    annotate(
     "text",
     x = quantile(loop_scatter_mass$Width_mm, 0.05, na.rm = TRUE),
     y = quantile(loop_scatter_mass$Mass_g,   0.95, na.rm = TRUE),
     label = loop_eq_mass,
     hjust = 0, vjust = 1, size = 3.5
    )
  }
  
  # ---- Dashed 50 mm line ONLY if the data reach >= 50 mm (Width on x) ----
  if (is.finite(max(loop_scatter_mass$Width_mm, na.rm = TRUE)) &&
      max(loop_scatter_mass$Width_mm, na.rm = TRUE) >= BAR_THRESHOLD_MM) {
   loop_p_scatter_mass <- loop_p_scatter_mass +
    geom_vline(xintercept = BAR_THRESHOLD_MM,
               linetype = "dashed", color = "firebrick",
               linewidth = 0.5, alpha = 0.8)
  }
  
 } else {
  loop_p_scatter_mass <- ggplot() + theme_void() +
   labs(caption = make_caption(loop_scatter_mass,
                               paste0("Insufficient data (n = ", nrow(loop_scatter_mass),
                                      " < ", MIN_N_PER_SPECIES, ") for Mass~Width plot."),
                               param_species))
 }
 
 # Store the plot and model row
 plots[[param_species]][["scatter_mass"]] <- loop_p_scatter_mass
 df_combined_models <- dplyr::bind_rows(df_combined_models, temp_model_row)
 
 #### Plot 3: Fork Length ~ Mass (Mass on Y) ####################
 ### Scatter with Mass on y, model still FL ~ ln(Mass); inverse curve + fitted-form equation ###
 
 # Use only rows with both FL and Mass present and Mass > 0
 loop_fl_mass <- df_clean %>%
  dplyr::filter(!is.na(ForkLength_mm), !is.na(Mass_g), Mass_g > 0)
 
 if (nrow(loop_fl_mass) >= MIN_N_PER_SPECIES) {
  # MASS on Y, FORK LENGTH on X
  loop_p_fl_mass <- ggplot(loop_fl_mass, aes(x = ForkLength_mm, y = Mass_g)) +
   geom_point(color = "#d95f02", alpha = 0.6, size = 2) +
   labs(
    title   = paste0(param_species, " - Mass by Fork Length"),
    x       = "Fork Length (mm)",
    y       = "Mass (g)",
    caption = make_caption(loop_fl_mass, "Mass vs fork length (model on ln(Mass)).", param_species)
   ) 
  
  # Fit requires ≥3 points
  if (nrow(loop_fl_mass) >= 3) {
   # Fit original model (predicting FL): FL = slope_lnx * ln(Mass) + int_lnx
   loop_lm_flm <- lm(ForkLength_mm ~ log(Mass_g), data = loop_fl_mass)
   loop_coef   <- coef(loop_lm_flm)
   slope_lnx   <- unname(loop_coef["log(Mass_g)"])  # slope on ln(Mass)
   int_lnx     <- unname(loop_coef["(Intercept)"])  # intercept
   r2_lnx      <- summary(loop_lm_flm)$r.squared
   
   # Invert for plotting on flipped axes:
   # FL = slope*ln(Mass) + intercept  ==>  Mass = exp((FL - intercept)/slope)
   x_rng_fl <- range(loop_fl_mass$ForkLength_mm, na.rm = TRUE)
   x_seq_fl <- seq(x_rng_fl[1], x_rng_fl[2], length.out = 200)
   pred_inv <- tibble(ForkLength_mm = x_seq_fl) %>%
    dplyr::mutate(Mass_pred = exp((ForkLength_mm - int_lnx) / slope_lnx))
   
   loop_p_fl_mass <- loop_p_fl_mass +
    geom_line(data = pred_inv, aes(x = ForkLength_mm, y = Mass_pred), color = "#bf5b17", linewidth = 1)
   
   # ---- Equation annotation (display with y = ..., to match other graphs) ----
   loop_eq <- paste0(
    "y = ", formatC(slope_lnx, format = "f", digits = 2), " ln(Mass)",
    ifelse(int_lnx >= 0, " + ", " - "),
    formatC(abs(int_lnx), format = "f", digits = 2),
    "\nR^2 = ", formatC(r2_lnx, format = "f", digits = 4)
   )
   
   loop_p_fl_mass <- loop_p_fl_mass +
    annotate(
     "text",
     x = quantile(loop_fl_mass$ForkLength_mm, 0.05, na.rm = TRUE),
     y = quantile(loop_fl_mass$Mass_g,        0.95, na.rm = TRUE),
     label = loop_eq,
     hjust = 0, vjust = 1, size = 3.5
    )
  }
 } else {
  loop_p_fl_mass <- ggplot() + theme_void() +
   labs(
    caption = make_caption(
     loop_fl_mass,
     paste0("Insufficient data (n = ", nrow(loop_fl_mass),
            " < ", MIN_N_PER_SPECIES, ") for Mass~Fork Length plot."),
     param_species
    )
   )
 }
 
 # Store the plot
 plots[[param_species]][["FL_by_mass"]] <- loop_p_fl_mass


 #### Plot 4: log(Width) ~ log(Fork Length) (log-log) --------#
 ### Minor: Scatter on original scale + log-log fit, line drawn on original scale ###
 
 # Use positive values only for log transformation
 loop_loglog_fl <- df_clean %>%
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(ForkLength_mm), ForkLength_mm > 0)
 
 if (nrow(loop_loglog_fl) >= MIN_N_PER_SPECIES) {
  loop_p_loglog_fl <- ggplot(loop_loglog_fl, aes(ForkLength_mm, Width_mm)) +
   geom_point(color = "#4daf4a", alpha = 0.6, size = 2) +
   scale_x_log10() + scale_y_log10() +
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
             label = paste0("log(y) = ", formatC(loop_slope_ll, format = "f", digits = 2),
                            " log(x) ", ifelse(loop_int_ll >= 0, "+ ", "- "),
                            formatC(abs(loop_int_ll), format = "f", digits = 2),
                            "\nR^2 = ", formatC(loop_r2_ll, format = "f", digits = 4)),
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
 
 
 #### Plot 5: Power law (Width = a * Mass^b) -----------------#
 ### Minor: Power-law fit via log-log regression (Width = a * Mass^b) ###
 
 # Use positive values only for log transformation
 loop_power_mass <- df_clean %>%
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(Mass_g),   Mass_g > 0)
 
 if (nrow(loop_power_mass) >= MIN_N_PER_SPECIES) {
  loop_p_power_mass <- ggplot(loop_power_mass, aes(Mass_g, Width_mm)) +
   geom_point(color = "#984ea3", alpha = 0.6, size = 2) +
   scale_x_log10() + scale_y_log10() +
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
             label = paste0("y = ", formatC(loop_a_pw, format = "f", digits = 2),
                            " x^", formatC(loop_b_pw, format = "f", digits = 2),
                            "\nR^2 = ", formatC(loop_r2_pw, format = "f", digits = 4)),
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
 
 
 #### Plot 6: Histogram of Fork Lengths per species -----------#
 ### Fixed binwidth (single width across all lengths) ###
 
 # Data for histogram
 loop_hist_df <- df_clean %>%
  dplyr::filter(!is.na(ForkLength_mm))
 
 # Only build if we have *some* data points (>= 1).
 if (nrow(loop_hist_df) >= 1) {
  
  # Mean & median markers
  loop_fl_mean   <- mean(loop_hist_df$ForkLength_mm, na.rm = TRUE)
  loop_fl_median <- stats::median(loop_hist_df$ForkLength_mm, na.rm = TRUE)
  
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
 cat("Finished:", if (nzchar(param_species)) param_species else "UNKNOWN_SPECIES", "\n")
} # end species loop

##### Multi-species plots #####################################----
#-------------------------------------------------------------#
# Build cross-species views using the combined dataset (df_all)

### Combined: Width ~ Fork Length (linear) ###
df_combined_fl <- df_all %>% dplyr::filter(!is.na(Width_mm), !is.na(ForkLength_mm))
plots[["combined"]][["width_by_FL"]] <- ggplot(df_combined_fl, aes(ForkLength_mm, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 geom_smooth(method = "lm", se = FALSE) +
 labs(title = "All species - Width by Fork Length (linear)",
      x = "Fork Length (mm)", y = "Width (mm)") +
 theme(legend.position = "bottom")

### Combined: Width ~ log(Mass) (log x) ###
df_combined_mass_logx <- df_all %>% dplyr::filter(!is.na(Width_mm), !is.na(Mass_g), Mass_g > 0)
plots[["combined"]][["width_by_mass_logx"]] <- ggplot(df_combined_mass_logx, aes(Mass_g, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 stat_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
 labs(title = "All species - Width by Mass (log x)",
      x = "Mass (g)", y = "Width (mm)") +
 theme(legend.position = "bottom")

### Combined: log(Width) ~ log(Fork Length) ###
df_combined_loglog_fl <- df_all %>%
 dplyr::filter(!is.na(Width_mm), Width_mm > 0,
               !is.na(ForkLength_mm), ForkLength_mm > 0)

# Base scatter + log scales
plots[["combined"]][["loglog_width_by_FL"]] <-
 ggplot(df_combined_loglog_fl, aes(ForkLength_mm, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 scale_x_log10() + scale_y_log10() +
 labs(title = "All species - log(Width) ~ log(Fork Length)",
      x = "Fork Length (mm, log10 scale)", y = "Width (mm, log10 scale)") +
 theme(legend.position = "bottom")

# Add trend(s) per toggle
if (COMBINED_TREND_MODE == "overall") {
 # One global line across all species on log-log (original behavior)
 plots[["combined"]][["loglog_width_by_FL"]] <- plots[["combined"]][["loglog_width_by_FL"]] +
  stat_smooth(method = "lm",
              formula = y ~ x, se = FALSE,
              mapping = aes(x = log(ForkLength_mm), y = log(Width_mm)),
              inherit.aes = FALSE, color = "black")
} else if (COMBINED_TREND_MODE == "per_species") {
 # One line per species (clearer)
 plots[["combined"]][["loglog_width_by_FL"]] <- plots[["combined"]][["loglog_width_by_FL"]] +
  geom_smooth(method = "lm", se = FALSE)
} # else "none": points only

### Combined: Power law (Width ~ Mass^b) ###
df_combined_power_mass <- df_all %>%
 dplyr::filter(!is.na(Width_mm), Width_mm > 0,
               !is.na(Mass_g),   Mass_g > 0)

# Base scatter + log scales
plots[["combined"]][["powerlaw_width_by_mass"]] <-
 ggplot(df_combined_power_mass, aes(Mass_g, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 scale_x_log10() + scale_y_log10() +
 labs(title = "All species - Power law: Width ~ Mass^b",
      x = "Mass (g, log10 scale)", y = "Width (mm, log10 scale)") +
 theme(legend.position = "bottom")

# Add trend(s) per toggle
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

### Combined: Histogram of Fork Lengths by Species ###
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

# Show plots in RStudio viewer
print(plots)


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

# Per-species plot export (compact)
# for (sp in names(combined_all)) for (nm in names(plots[[sp]]))
# ggsave(file.path(path_figs_dir, paste0(gsub("[^A-Za-z0-9_\\-]", "_", sp), "_", nm, ".png")),
# plots[[sp]][[nm]], width = 7, height = 5, dpi = 300)

##### Cleanup (remove temporary loop/temp objects) #############----
#-------------------------------------------------------------#
# Remove transient objects created inside the loop or as temp_* helpers.
loop_temp_vars <- ls(pattern = "^(loop_|temp_)")
if (length(loop_temp_vars)) {
 rm(list = loop_temp_vars)
}

cat("\nAll plots complete.\n")
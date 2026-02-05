## --------------------------------------------------------------#
## Script name: script01-02_morphology_plots.R
##
## Purpose:
##   - Use pre-cleaned data (from script01-02_import_format.R) to
##     generate morphology plots and simple regressions per species,
##     plus multi-species plots for each plot class.
##
##
#### Repository layout (per template):
##   01_data/                       # inputs
##   01_data/processed/<species>/   # cleaned RDS per species
##   03_outputs/01_figures/         # figures
##   03_outputs/01_tables/          # tables (created if missing)
##
##
## Author: Marcus Rizzuto
## Date Created: 1/27/2026
## --------------------------------------------------------------#


##### Guards & Setup ##########################################----
#-------------------------------------------------------------#
# Require cleaned objects from the import script:
# - combined_all : list of cleaned data frames by species
# - df_all       : single combined data frame of all species
if (!exists("combined_all")) {
 stop("ERROR: `combined_all` not found. Run script01-02_import_format.R first.")
}
if (!exists("df_all")) {
 stop("ERROR: `df_all` not found. Run script01-02_import_format.R first.")
}

# Output directories (created if missing; exports remain commented out below)
path_figs_dir   <- file.path("03_outputs", "01_figures")
path_tables_dir <- file.path("03_outputs", "01_tables")
if (!dir.exists(path_figs_dir))   dir.create(path_figs_dir,   recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(path_tables_dir)) dir.create(path_tables_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize containers for results and plots
df_combined_summary <- tibble()
df_combined_models  <- tibble()  # collects coefficients for log(Mass) model (for continuity)
plots <- list()   # nested list of plots
plots[["combined"]] <- list()    # ensure multi-species plot container exists

# Sanitize species names in case of empty keys (defensive against zero-length names)
sp_names <- names(combined_all)
bad_idx <- which(is.na(sp_names) | !nzchar(sp_names))
if (length(bad_idx)) {
 warning("Empty species names detected at positions: ",
         paste(bad_idx, collapse = ", "),
         ". Renaming to 'UNKNOWN_<idx>'.")
 sp_names[bad_idx] <- paste0("UNKNOWN_", bad_idx)
 names(combined_all) <- sp_names
}

# Helper to standardize plot captions per species
make_caption <- function(df, text, sp) paste0(sp, " (n = ", nrow(df), "): ", text)


##### Species Loop (plots + simple models only) ###############----
#-------------------------------------------------------------#
# For each species:
#   - Compute summary statistics
#   - Build 4 plot types (with model annotations when data allow)
#   - Store per-species plots in `plots[[<species>]]`
for (param_species in names(combined_all)) {
 
 cat("\n--- Plotting:", if (nzchar(param_species)) param_species else "UNKNOWN_SPECIES", "---\n")
 df_clean <- combined_all[[param_species]]
 
 # Ensure per-species plot sublist exists before assigning plots
 if (is.null(plots[[param_species]]) || !is.list(plots[[param_species]])) {
  plots[[param_species]] <- list()
 }
 
 # Summary stats (row counts and means for key morphometrics)
 df_summary <- df_clean %>%
  dplyr::summarise(
   species       = param_species,
   n_rows        = dplyr::n(),
   n_FL          = sum(!is.na(ForkLength_mm)),
   n_width       = sum(!is.na(Width_mm)),
   n_weight      = sum(!is.na(Mass_g)),
   FL_mean_mm    = mean(ForkLength_mm, na.rm = TRUE),
   width_mean_mm = mean(Width_mm,      na.rm = TRUE),
   weight_mean_g = mean(Mass_g,        na.rm = TRUE)
  )
 df_combined_summary <- dplyr::bind_rows(df_combined_summary, df_summary)
 
 #### Plot 1: Width ~ Fork Length (linear) ------------------#
 # Scatter with optional linear fit and equation annotation
 df_scatter_fl <- df_clean %>% dplyr::filter(!is.na(Width_mm), !is.na(ForkLength_mm))
 
 p_scatter_fl <- ggplot(df_scatter_fl, aes(ForkLength_mm, Width_mm)) +
  geom_point(color = "#2c7fb8", alpha = 0.6, size = 2) +
  labs(
   title   = paste0(param_species, " - Width by Fork Length"),
   x       = "Fork Length (mm)",
   y       = "Width (mm)",
   caption = make_caption(df_scatter_fl, "Width vs fork length.", param_species)
  ) +
  theme_minimal()
 
 if (nrow(df_scatter_fl) >= 2) {
  temp_lm_fl  <- lm(Width_mm ~ ForkLength_mm, data = df_scatter_fl)
  temp_coef   <- coef(temp_lm_fl)
  temp_r2_fl  <- summary(temp_lm_fl)$r.squared
  temp_slope  <- unname(temp_coef[["ForkLength_mm"]])
  temp_int    <- unname(temp_coef[["(Intercept)"]])
  
  temp_eq_fl <- paste0(
   "y = ", formatC(temp_slope, format = "f", digits = 2),
   "x", ifelse(temp_int >= 0, " + ", " - "),
   formatC(abs(temp_int), format = "f", digits = 2),
   "\nR^2 = ", formatC(temp_r2_fl, format = "f", digits = 4)
  )
  
  p_scatter_fl <- p_scatter_fl +
   geom_smooth(method = "lm", se = FALSE, color = "#1f78b4", linewidth = 0.9) +
   annotate("text",
            x = quantile(df_scatter_fl$ForkLength_mm, 0.05, na.rm = TRUE),
            y = quantile(df_scatter_fl$Width_mm,      0.95, na.rm = TRUE),
            label = temp_eq_fl, hjust = 0, vjust = 1, size = 3.5)
 }
 
 plots[[param_species]][["scatter_fl"]] <- p_scatter_fl
 
 #### Plot 2: Width ~ log(Mass) (logarithmic x) --------------#
 # Scatter with log(Mass) linear model; store alpha/beta/R² per species
 df_scatter_mass <- df_clean %>% dplyr::filter(!is.na(Width_mm), !is.na(Mass_g), Mass_g > 0)
 
 # Initialize model row (remains NA if too few points)
 temp_model_row <- tibble(species = param_species, alpha = NA_real_, beta = NA_real_, r2 = NA_real_)
 
 p_scatter_mass <- ggplot(df_scatter_mass, aes(Mass_g, Width_mm)) +
  geom_point(color = "#f16913", alpha = 0.6, size = 2) +
  labs(
   title   = paste0(param_species, " - Width by Mass"),
   x       = "Mass (g)",
   y       = "Width (mm)",
   caption = make_caption(df_scatter_mass, "Width vs mass.", param_species)
  ) +
  theme_minimal()
 
 if (nrow(df_scatter_mass) >= 3) {
  temp_log_fit <- lm(Width_mm ~ log(Mass_g), data = df_scatter_mass)
  temp_coefs   <- coef(temp_log_fit)
  temp_alpha   <- unname(temp_coefs["log(Mass_g)"])
  temp_beta    <- unname(temp_coefs["(Intercept)"])
  temp_r2_m    <- summary(temp_log_fit)$r.squared
  
  temp_model_row <- tibble(species = param_species, alpha = temp_alpha, beta = temp_beta, r2 = temp_r2_m)
  
  temp_eq_m <- paste0(
   "y = ", formatC(temp_alpha, format = "f", digits = 2),
   " ln(x)", ifelse(temp_beta >= 0, " + ", " - "),
   formatC(abs(temp_beta), format = "f", digits = 2),
   "\nR^2 = ", formatC(temp_r2_m, format = "f", digits = 4)
  )
  
  temp_xrng <- range(df_scatter_mass$Mass_g, na.rm = TRUE)
  temp_xseq <- seq(temp_xrng[1], temp_xrng[2], length.out = 200)
  df_pred   <- tibble(Mass_g = temp_xseq) %>%
   dplyr::mutate(Width_pred = predict(temp_log_fit, newdata = tibble(Mass_g = Mass_g)))
  
  p_scatter_mass <- p_scatter_mass +
   geom_line(data = df_pred, aes(Mass_g, Width_pred), color = "#cc4c02", linewidth = 1) +
   annotate("text",
            x = quantile(df_scatter_mass$Mass_g, 0.05, na.rm = TRUE),
            y = quantile(df_scatter_mass$Width_mm, 0.95, na.rm = TRUE),
            label = temp_eq_m, hjust = 0, vjust = 1, size = 3.5)
 }
 
 plots[[param_species]][["scatter_mass"]] <- p_scatter_mass
 df_combined_models <- dplyr::bind_rows(df_combined_models, temp_model_row)
 
 #### Plot 3: log(Width) ~ log(Fork Length) (log-log) --------#
 # Scatter on original scale with log-log fit annotated (line drawn on original scale)
 df_loglog_fl <- df_clean %>%
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(ForkLength_mm), ForkLength_mm > 0)
 
 p_loglog_fl <- ggplot(df_loglog_fl, aes(ForkLength_mm, Width_mm)) +
  geom_point(color = "#4daf4a", alpha = 0.6, size = 2) +
  scale_x_log10() + scale_y_log10() +
  labs(
   title   = paste0(param_species, " - log(Width) ~ log(Fork Length)"),
   x       = "Fork Length (mm, log10 scale)",
   y       = "Width (mm, log10 scale)",
   caption = make_caption(df_loglog_fl, "log(Width) vs log(ForkLength).", param_species)
  ) +
  theme_minimal()
 
 if (nrow(df_loglog_fl) >= 3) {
  temp_fit_ll <- lm(log(Width_mm) ~ log(ForkLength_mm), data = df_loglog_fl)
  temp_coefs  <- coef(temp_fit_ll)
  temp_r2_ll  <- summary(temp_fit_ll)$r.squared
  temp_slope  <- unname(temp_coefs[["log(ForkLength_mm)"]])
  temp_int    <- unname(temp_coefs[["(Intercept)"]])
  
  rng_x <- range(df_loglog_fl$ForkLength_mm, na.rm = TRUE)
  xseq  <- seq(rng_x[1], rng_x[2], length.out = 200)
  df_pred_ll <- tibble(ForkLength_mm = xseq) %>%
   dplyr::mutate(Width_pred = exp(temp_int + temp_slope * log(ForkLength_mm)))
  
  p_loglog_fl <- p_loglog_fl +
   geom_line(data = df_pred_ll, aes(ForkLength_mm, Width_pred),
             color = "#238b45", linewidth = 1) +
   annotate("text",
            x = quantile(df_loglog_fl$ForkLength_mm, 0.05, na.rm = TRUE),
            y = quantile(df_loglog_fl$Width_mm,       0.95, na.rm = TRUE),
            label = paste0("log(y) = ", formatC(temp_slope, format = "f", digits = 2),
                           " log(x) ", ifelse(temp_int >= 0, "+ ", "- "),
                           formatC(abs(temp_int), format = "f", digits = 2),
                           "\nR^2 = ", formatC(temp_r2_ll, format = "f", digits = 4)),
            hjust = 0, vjust = 1, size = 3.5)
 }
 
 plots[[param_species]][["loglog_width_by_FL"]] <- p_loglog_fl
 
 #### Plot 4: Power law (Width = a * Mass^b) -----------------#
 # Power-law curve via log-log regression: Width = a * Mass^b
 df_power_mass <- df_clean %>%
  dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                !is.na(Mass_g),   Mass_g > 0)
 
 p_power_mass <- ggplot(df_power_mass, aes(Mass_g, Width_mm)) +
  geom_point(color = "#984ea3", alpha = 0.6, size = 2) +
  scale_x_log10() + scale_y_log10() +
  labs(
   title   = paste0(param_species, " - Power law: Width ~ Mass^b"),
   x       = "Mass (g, log10 scale)",
   y       = "Width (mm, log10 scale)",
   caption = make_caption(df_power_mass, "Power-law width vs mass.", param_species)
  ) +
  theme_minimal()
 
 if (nrow(df_power_mass) >= 3) {
  temp_fit_pw <- lm(log(Width_mm) ~ log(Mass_g), data = df_power_mass)
  temp_coefs  <- coef(temp_fit_pw)
  temp_r2_pw  <- summary(temp_fit_pw)$r.squared
  temp_b      <- unname(temp_coefs[["log(Mass_g)"]])  # exponent b
  temp_log_a  <- unname(temp_coefs[["(Intercept)"]])
  temp_a      <- exp(temp_log_a)
  
  rng_x <- range(df_power_mass$Mass_g, na.rm = TRUE)
  xseq  <- seq(rng_x[1], rng_x[2], length.out = 200)
  df_pred_pw <- tibble(Mass_g = xseq) %>%
   dplyr::mutate(Width_pred = temp_a * (Mass_g ^ temp_b))
  
  p_power_mass <- p_power_mass +
   geom_line(data = df_pred_pw, aes(Mass_g, Width_pred),
             color = "#6a51a3", linewidth = 1) +
   annotate("text",
            x = quantile(df_power_mass$Mass_g, 0.05, na.rm = TRUE),
            y = quantile(df_power_mass$Width_mm, 0.95, na.rm = TRUE),
            label = paste0("y = ", formatC(temp_a, format = "f", digits = 2),
                           " x^", formatC(temp_b, format = "f", digits = 2),
                           "\nR^2 = ", formatC(temp_r2_pw, format = "f", digits = 4)),
            hjust = 0, vjust = 1, size = 3.5)
 }
 
 plots[[param_species]][["powerlaw_width_by_mass"]] <- p_power_mass
 
 cat("Finished:", if (nzchar(param_species)) param_species else "UNKNOWN_SPECIES", "\n")
} # end species loop


##### Multi-species plots #####################################----
#-------------------------------------------------------------#
# Build cross-species views using the combined dataset (df_all)

# 1) Width ~ Fork Length (linear)
df_all1 <- df_all %>% dplyr::filter(!is.na(Width_mm), !is.na(ForkLength_mm))
plots[["combined"]][["width_by_FL"]] <- ggplot(df_all1, aes(ForkLength_mm, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 geom_smooth(method = "lm", se = FALSE) +
 labs(title = "All species - Width by Fork Length (linear)",
      x = "Fork Length (mm)", y = "Width (mm)") +
 theme(legend.position = "bottom")

# 2) Width ~ log(Mass)
df_all2 <- df_all %>% dplyr::filter(!is.na(Width_mm), !is.na(Mass_g), Mass_g > 0)
plots[["combined"]][["width_by_mass_logx"]] <- ggplot(df_all2, aes(Mass_g, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 stat_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
 labs(title = "All species - Width by Mass (log x)",
      x = "Mass (g)", y = "Width (mm)") +
 theme(legend.position = "bottom")

# 3) log(Width) ~ log(Fork Length)
df_all3 <- df_all %>% dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                                    !is.na(ForkLength_mm), ForkLength_mm > 0)
plots[["combined"]][["loglog_width_by_FL"]] <- ggplot(df_all3, aes(ForkLength_mm, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 scale_x_log10() + scale_y_log10() +
 stat_smooth(method = "lm",
             formula = y ~ x, se = FALSE,
             mapping = aes(x = log(ForkLength_mm), y = log(Width_mm)),
             inherit.aes = FALSE, color = "black") +
 labs(title = "All species - log(Width) ~ log(Fork Length)",
      x = "Fork Length (mm, log10 scale)", y = "Width (mm, log10 scale)") +
 theme(legend.position = "bottom")

# 4) Power law: Width ~ Mass^b
df_all4 <- df_all %>% dplyr::filter(!is.na(Width_mm), Width_mm > 0,
                                    !is.na(Mass_g),   Mass_g > 0)
plots[["combined"]][["powerlaw_width_by_mass"]] <- ggplot(df_all4, aes(Mass_g, Width_mm, color = Species)) +
 geom_point(alpha = 0.5, size = 1.8) +
 scale_x_log10() + scale_y_log10() +
 stat_smooth(method = "lm",
             formula = y ~ x, se = FALSE,
             mapping = aes(x = log(Mass_g), y = log(Width_mm)),
             inherit.aes = FALSE, color = "black") +
 labs(title = "All species - Power law: Width ~ Mass^b",
      x = "Mass (g, log10 scale)", y = "Width (mm, log10 scale)") +
 theme(legend.position = "bottom")

# Show plots in RStudio viewer
print(plots)

##### Optional exports (kept commented) #######################----
#-------------------------------------------------------------#
# readr::write_csv(df_combined_summary, file.path(path_tables_dir, "_combined_summary_by_species.csv"))
# readr::write_csv(df_combined_models,  file.path(path_tables_dir, "_combined_log_model_coefficients.csv"))
# ggsave(file.path(path_figs_dir, "combined_width_by_FL.png"),        plots[["combined"]][["width_by_FL"]],        width = 7, height = 5, dpi = 300)
# ggsave(file.path(path_figs_dir, "combined_width_by_mass_logx.png"), plots[["combined"]][["width_by_mass_logx"]], width = 7, height = 5, dpi = 300)
# ggsave(file.path(path_figs_dir, "combined_loglog_width_by_FL.png"), plots[["combined"]][["loglog_width_by_FL"]], width = 7, height = 5, dpi = 300)
# ggsave(file.path(path_figs_dir, "combined_powerlaw_width_by_mass.png"), plots[["combined"]][["powerlaw_width_by_mass"]], width = 7, height = 5, dpi = 300)

cat("\nAll plots complete.\n")

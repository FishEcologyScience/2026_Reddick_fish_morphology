## --------------------------------------------------------------#
## Script name: script03_binning_FD_legacy.R
##
## Purpose:
##   - Preserve the original adaptive (FD‑style) binning logic used
##     in the morphology workflow prior to fixed‑bin standardization.
##   - Compute and print the legacy FD binwidth and binned Fork
##     Length counts for each species (no plotting, no saving).
##
## Author: Marcus Rizzuto
## Date Created: 3/5/2026
##
## Notes:
##   - This script is intentionally limited to binning only.
##   - Plotting and visualization occur exclusively in 
##     script01-02_morphology_plots.R.
## --------------------------------------------------------------#

library(dplyr)

# df_all MUST already exist (loaded from the import/format script)
if (!exists("df_all")) {
 stop("ERROR: `df_all` not found. Run the import/format script first.")
}

# Function: original safe FD binwidth
fd_binwidth_safe <- function(x, fallback = 5) {
 x <- x[is.finite(x)]
 n <- length(x)
 if (n <= 1) return(fallback)
 iqr <- IQR(x, na.rm = TRUE)
 bw  <- 2 * iqr / (n)^(1/3)
 if (is.finite(bw) && bw > 0) bw else fallback
}

# Compute and PRINT legacy FD binning for each species
for (sp in unique(df_all$Species)) {
 
 df_sp <- df_all %>%
  filter(Species == sp, !is.na(ForkLength_mm))
 
 if (nrow(df_sp) == 0) {
  cat("\nSpecies:", sp, "- no FL data.\n")
  next
 }
 
 x <- df_sp$ForkLength_mm
 bw <- fd_binwidth_safe(x)
 
 rng <- range(x)
 breaks <- seq(
  floor(rng[1] / bw) * bw,
  ceiling(rng[2] / bw) * bw + bw,
  by = bw
 )
 
 bins <- cut(x, breaks = breaks, include.lowest = TRUE, right = FALSE)
 
 cat("\n-----------------------------\n")
 cat("Species:", sp, "\n")
 cat("FD binwidth =", format(bw, digits = 3), "mm\n")
 print(table(bins))
}


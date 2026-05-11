library(tidyverse)

poor_naming_scheme_range <- df_all %>%
  dplyr::filter(!is.na(ForkLength_mm)) %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    min_fl = min(ForkLength_mm),
    max_fl = max(ForkLength_mm)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Species = forcats::fct_reorder(Species, max_fl))

poor_naming_scheme_lollipop_size_range <- ggplot(poor_naming_scheme_range, aes(y = Species)) +
  geom_segment(aes(x = min_fl, xend = max_fl, yend = Species),
               linewidth = 1, colour = "grey50") +
  geom_point(aes(x = min_fl), size = 3, colour = "#2166ac") +
  geom_point(aes(x = max_fl), size = 3, colour = "#d6604d") +
  labs(
    x = "Fork Length (mm)",
    y = NULL
  ) +
  theme_classic()

poor_naming_scheme_lollipop_size_range

# ================================================================
# VISUALISE BIRTH YEARS OF FOUNDERS AND ANCESTORS
# SOUTH AFRICAN BOERPERD
#
# Requires:
# founder_list
# ancestor_list
#
# Produced by the PurgeR founder/ancestor script
# ================================================================


# ================================================================
# 1. PACKAGES
# ================================================================

# Install once if required:
# install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


# ================================================================
# 2. CHECK THAT THE DATA OBJECTS EXIST
# ================================================================

sum(is.na(founder_list$BirthYear))
sum(is.na(ancestor_list$BirthYear))

range(founder_list$BirthYear, na.rm = TRUE)
range(ancestor_list$BirthYear, na.rm = TRUE)

if (!exists("founder_list")) {
  
  founder_list <- read.csv(
    "PurgeR_REG_founders_birthdates.csv",
    stringsAsFactors = FALSE
  )
}


if (!exists("ancestor_list")) {
  
  ancestor_list <- read.csv(
    "PurgeR_REG_ancestors_birthdates.csv",
    stringsAsFactors = FALSE
  )
}


# Make sure BirthYear is numeric

founder_list$BirthYear <- as.numeric(
  founder_list$BirthYear
)

ancestor_list$BirthYear <- as.numeric(
  ancestor_list$BirthYear
)


# ================================================================
# 3. BASIC COUNTS
# ================================================================

n_founders_total <- nrow(founder_list)

n_founders_known <- sum(
  !is.na(founder_list$BirthYear)
)

n_founders_unknown <- sum(
  is.na(founder_list$BirthYear)
)


n_ancestors_total <- nrow(ancestor_list)

n_ancestors_known <- sum(
  !is.na(ancestor_list$BirthYear)
)

n_ancestors_unknown <- sum(
  is.na(ancestor_list$BirthYear)
)


founder_known_pct <- 100 *
  n_founders_known /
  n_founders_total


ancestor_known_pct <- 100 *
  n_ancestors_known /
  n_ancestors_total


cat(
  "\nFounders:",
  n_founders_total,
  "\nKnown birth year:",
  n_founders_known,
  paste0(" (", round(founder_known_pct, 1), "%)"),
  "\nUnknown birth year:",
  n_founders_unknown,
  "\n"
)


cat(
  "\nAncestors:",
  n_ancestors_total,
  "\nKnown birth year:",
  n_ancestors_known,
  paste0(" (", round(ancestor_known_pct, 1), "%)"),
  "\nUnknown birth year:",
  n_ancestors_unknown,
  "\n"
)


# ================================================================
# 4. COMBINE FOUNDER AND ANCESTOR BIRTH YEARS
# ================================================================

birthyear_data <- bind_rows(
  
  founder_list %>%
    transmute(
      ANI_ID = ANI_ID,
      BirthYear = BirthYear,
      Group = "Founders"
    ),
  
  ancestor_list %>%
    transmute(
      ANI_ID = ANI_ID,
      BirthYear = BirthYear,
      Group = "Ancestors"
    )
)


birthyear_data$Group <- factor(
  birthyear_data$Group,
  levels = c(
    "Founders",
    "Ancestors"
  )
)


# ================================================================
# 5. LIMIT TO ANIMALS WITH KNOWN BIRTH YEARS
# ================================================================

birthyear_known <- birthyear_data %>%
  filter(
    !is.na(BirthYear)
  )


# ================================================================
# 6. CREATE 5-YEAR BINS
# ================================================================
#
# Example:
#
# 1948 -> 1945-1949
# 1950 -> 1950-1954
# 1973 -> 1970-1974
#
# BinStart is kept numeric so historical reference years
# can be positioned correctly on the x-axis.
# ================================================================

birthyear_5yr <- birthyear_known %>%
  
  mutate(
    BinStart = floor(BirthYear / 5) * 5
  ) %>%
  
  count(
    Group,
    BinStart,
    name = "Number"
  ) %>%
  
  group_by(
    Group
  ) %>%
  
  mutate(
    Percentage = 100 * Number / sum(Number)
  ) %>%
  
  ungroup()


# Inspect the table

print(
  birthyear_5yr
)


# Save the table used to make the figure

write.csv(
  birthyear_5yr,
  "PurgeR_founder_ancestor_birthyear_5year_distribution.csv",
  row.names = FALSE
)


# ================================================================
# 7. SET COMMON X-AXIS RANGE
# ================================================================

min_year <- floor(
  min(
    birthyear_known$BirthYear,
    na.rm = TRUE
  ) / 5
) * 5


max_year <- ceiling(
  max(
    birthyear_known$BirthYear,
    na.rm = TRUE
  ) / 5
) * 5


cat(
  "\nPlot year range:",
  min_year,
  "-",
  max_year,
  "\n"
)


# ================================================================
# 8. SET COMMON Y-AXIS FOR PANELS A AND B
# ================================================================

distribution_ymax <- ceiling(
  max(
    birthyear_5yr$Percentage,
    na.rm = TRUE
  ) / 5
) * 5


# Add a little space for visual clarity
distribution_ymax <- distribution_ymax + 5


# ================================================================
# 9. PANEL A — FOUNDER BIRTH-YEAR DISTRIBUTION
# ================================================================

founder_plot_data <- birthyear_5yr %>%
  filter(
    Group == "Founders"
  )


plot_A <- ggplot(
  founder_plot_data,
  aes(
    x = BinStart + 2.5,
    y = Percentage
  )
) +
  
  geom_col(
    width = 4.5,
    fill = "#377EB8"
  ) +
  
  # Historical markers
  geom_vline(
    xintercept = 1973,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_vline(
    xintercept = 1998,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  scale_x_continuous(
    limits = c(
      min_year,
      max_year
    ),
    breaks = seq(
      1950,
      max_year,
      by = 10
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      distribution_ymax
    ),
    breaks = seq(
      0,
      distribution_ymax,
      by = 5
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  labs(
    title = "A. Birth-year distribution of founders",
    subtitle = paste0(
      "Known birth year: ",
      n_founders_known,
      " of ",
      n_founders_total,
      " founders (",
      round(founder_known_pct, 1),
      "%)"
    ),
    x = "Birth year",
    y = "Percentage of founders\nwith known birth year"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.title = element_text(
      face = "bold"
    )
  )


# ================================================================
# 10. PANEL B — ANCESTOR BIRTH-YEAR DISTRIBUTION
# ================================================================

ancestor_plot_data <- birthyear_5yr %>%
  filter(
    Group == "Ancestors"
  )


plot_B <- ggplot(
  ancestor_plot_data,
  aes(
    x = BinStart + 2.5,
    y = Percentage
  )
) +
  
  geom_col(
    width = 4.5,
    fill = "#D55E00"
  ) +
  
  # Historical markers
  geom_vline(
    xintercept = 1973,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_vline(
    xintercept = 1998,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  scale_x_continuous(
    limits = c(
      min_year,
      max_year
    ),
    breaks = seq(
      1950,
      max_year,
      by = 10
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      distribution_ymax
    ),
    breaks = seq(
      0,
      distribution_ymax,
      by = 5
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  labs(
    title = "B. Birth-year distribution of ancestors",
    subtitle = paste0(
      "Known birth year: ",
      n_ancestors_known,
      " of ",
      n_ancestors_total,
      " ancestors (",
      round(ancestor_known_pct, 1),
      "%)"
    ),
    x = "Birth year",
    y = "Percentage of ancestors\nwith known birth year"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.title = element_text(
      face = "bold"
    )
  )


# ================================================================
# 11. CREATE CUMULATIVE BIRTH-YEAR DATA
# ================================================================

cumulative_data <- birthyear_known %>%
  
  count(
    Group,
    BirthYear,
    name = "Number"
  ) %>%
  
  group_by(
    Group
  ) %>%
  
  complete(
    BirthYear = seq(
      min_year,
      max_year,
      by = 1
    ),
    fill = list(
      Number = 0
    )
  ) %>%
  
  arrange(
    BirthYear,
    .by_group = TRUE
  ) %>%
  
  mutate(
    CumulativeNumber = cumsum(Number),
    CumulativePercentage =
      100 * CumulativeNumber / sum(Number)
  ) %>%
  
  ungroup()


# Save cumulative data

write.csv(
  cumulative_data,
  "PurgeR_founder_ancestor_cumulative_birthyears.csv",
  row.names = FALSE
)


# ================================================================
# 12. PANEL C — CUMULATIVE DISTRIBUTION
# ================================================================

plot_C <- ggplot(
  cumulative_data,
  aes(
    x = BirthYear,
    y = CumulativePercentage,
    colour = Group
  )
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  # Historical markers
  geom_vline(
    xintercept = 1973,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_vline(
    xintercept = 1998,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  annotate(
    "text",
    x = 1973,
    y = 98,
    label = "1973",
    angle = 90,
    vjust = -0.4,
    hjust = 1,
    size = 3.2
  ) +
  
  annotate(
    "text",
    x = 1998,
    y = 98,
    label = "1998",
    angle = 90,
    vjust = -0.4,
    hjust = 1,
    size = 3.2
  ) +
  
  scale_colour_manual(
    values = c(
      "Founders" = "#377EB8",
      "Ancestors" = "#D55E00"
    )
  ) +
  
  scale_x_continuous(
    limits = c(
      min_year,
      max_year
    ),
    breaks = seq(
      1950,
      max_year,
      by = 10
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  labs(
    title = "C. Cumulative distribution of founder and ancestor birth years",
    x = "Birth year",
    y = "Cumulative percentage",
    colour = NULL
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )


# ================================================================
# 13. COMBINE ALL THREE PANELS
# ================================================================

birthyear_figure <- (
  plot_A /
    plot_B /
    plot_C
) +
  
  plot_layout(
    heights = c(
      1,
      1,
      1.15
    )
  ) +
  
  plot_annotation(
    
    title = "Birth-year distribution of founders and ancestors in the South African Boerperd pedigree",
    
    caption = paste0(
      "Birth years were available for ",
      n_founders_known,
      " of ",
      n_founders_total,
      " founders (",
      round(founder_known_pct, 1),
      "%) and ",
      n_ancestors_known,
      " of ",
      n_ancestors_total,
      " ancestors (",
      round(ancestor_known_pct, 1),
      "%). ",
      "Individuals with unknown birth years were excluded from the temporal distributions. ",
      "Dashed vertical lines indicate 1973 and 1998."
    ),
    
    theme = theme(
      
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      ),
      
      plot.caption = element_text(
        size = 9,
        hjust = 0
      )
    )
  )


# ================================================================
# 14. DISPLAY FIGURE
# ================================================================

print(
  birthyear_figure
)


# ================================================================
# 15. SAVE FIGURE
# ================================================================

ggsave(
  filename =
    "PurgeR_founder_ancestor_birthyear_figure.png",
  plot = birthyear_figure,
  width = 10,
  height = 13,
  units = "in",
  dpi = 600
)


ggsave(
  filename =
    "PurgeR_founder_ancestor_birthyear_figure.pdf",
  plot = birthyear_figure,
  width = 10,
  height = 13,
  units = "in"
)


# ================================================================
# 16. CALCULATE HISTORICALLY USEFUL PERCENTAGES
# ================================================================
#
# These values will help us interpret Panel C.
#
# We calculate:
# percentage born by 1973
# percentage born by 1998
#
# IMPORTANT:
# These percentages apply ONLY to animals with known birth years.
# ================================================================


historical_summary <- birthyear_known %>%
  
  group_by(
    Group
  ) %>%
  
  summarise(
    
    Known_birth_years = n(),
    
    Born_by_1973 =
      sum(
        BirthYear <= 1973
      ),
    
    Percent_born_by_1973 =
      100 *
      mean(
        BirthYear <= 1973
      ),
    
    Born_by_1998 =
      sum(
        BirthYear <= 1998
      ),
    
    Percent_born_by_1998 =
      100 *
      mean(
        BirthYear <= 1998
      ),
    
    Median_birth_year =
      median(
        BirthYear
      ),
    
    .groups = "drop"
  )


print(
  historical_summary
)


write.csv(
  historical_summary,
  "PurgeR_founder_ancestor_historical_birthyear_summary.csv",
  row.names = FALSE
)


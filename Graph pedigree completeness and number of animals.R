# ================================================================
# PEDIGREE COMPLETENESS AND NUMBER OF ANIMALS BY YEAR
# ================================================================


# Install packages once if needed:
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "scales"))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)


# ================================================================
# 1. SET WORKING DIRECTORY
# ================================================================

setwd(
  "C:/2025/Paper drafts/Pedigree paper/1973 cohort analysis check"
)


# ================================================================
# READ AND CLEAN THE DATA
# ================================================================

ped_year <- read_excel(
  "POPREP inbreeding TABLE 1 year number animals pedigree completeness.xlsx",
  skip = 1,
  col_types = "text"
)


# Replace the unusual middle dot with a decimal point
# and remove spaces in all completeness columns

ped_year <- ped_year %>%
  mutate(
    across(
      starts_with("Complete gen"),
      ~ gsub(
        "·",
        ".",
        gsub(
          " ",
          "",
          .
        )
      )
    )
  )


# Convert everything to numeric

ped_year <- ped_year %>%
  mutate(
    Year = as.numeric(Year),
    `No of animals` = as.numeric(`No of animals`),
    across(
      starts_with("Complete gen"),
      as.numeric
    )
  )


# Keep years up to 2024

ped_year <- ped_year %>%
  filter(
    Year <= 2024
  )


# Check the data

head(ped_year)

tail(ped_year)

colSums(
  is.na(
    ped_year[
      ,
      grep(
        "^Complete gen",
        names(ped_year)
      )
    ]
  )
)

# ================================================================
# PEDIGREE COMPLETENESS AND NUMBER OF ANIMALS BY YEAR
# ================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)


# ================================================================
# 1. CONVERT COMPLETENESS DATA TO LONG FORMAT
# ================================================================

complete_long <- ped_year %>%
  pivot_longer(
    cols = starts_with("Complete gen"),
    names_to = "Generation",
    values_to = "Completeness"
  )


# Keep generations in the correct order in the legend
complete_long$Generation <- factor(
  complete_long$Generation,
  levels = paste0(
    "Complete gen ",
    1:6
  )
)


# ================================================================
# 2. CALCULATE SCALING FOR NUMBER OF ANIMALS
# ================================================================

# Pedigree completeness is plotted from 0 to 100%.
# Number of animals is temporarily rescaled onto this same
# plotting range, then converted back to real animal numbers
# on the right-hand y-axis.

max_animals <- max(
  ped_year$`No of animals`,
  na.rm = TRUE
)

animal_scale <- 100 / max_animals


# ================================================================
# CREATE THE GRAPH
# ================================================================

pedigree_plot <- ggplot() +
  
  # --------------------------------------------------------------
# SIX PEDIGREE COMPLETENESS LINES
# --------------------------------------------------------------

geom_line(
  data = complete_long,
  aes(
    x = Year,
    y = Completeness,
    colour = Generation,
    group = Generation
  ),
  linewidth = 1
) +
  
  # --------------------------------------------------------------
# NUMBER OF ANIMALS
# Solid black line using the right-hand y-axis
# --------------------------------------------------------------

geom_line(
  data = ped_year,
  aes(
    x = Year,
    y = `No of animals` * animal_scale,
    colour = "Number of animals"
  ),
  linewidth = 1.3
) +
  
  # --------------------------------------------------------------
# COLOURS
# --------------------------------------------------------------

scale_colour_manual(
  values = c(
    "Complete gen 1" = "#E41A1C",
    "Complete gen 2" = "#377EB8",
    "Complete gen 3" = "#4DAF4A",
    "Complete gen 4" = "#984EA3",
    "Complete gen 5" = "#FF7F00",
    "Complete gen 6" = "#A65628",
    "Number of animals" = "black"
  ),
  
  breaks = c(
    "Complete gen 1",
    "Complete gen 2",
    "Complete gen 3",
    "Complete gen 4",
    "Complete gen 5",
    "Complete gen 6",
    "Number of animals"
  )
) +
  
  # --------------------------------------------------------------
# X AXIS
# --------------------------------------------------------------

scale_x_continuous(
  name = "Year",
  
  limits = c(
    min(
      ped_year$Year,
      na.rm = TRUE
    ),
    2024
  ),
  
  breaks = c(
    seq(
      ceiling(
        min(
          ped_year$Year,
          na.rm = TRUE
        ) / 10
      ) * 10,
      2020,
      by = 10
    ),
    2024
  ),
  
  expand = expansion(
    mult = c(
      0.01,
      0.01
    )
  )
) +
  
  # --------------------------------------------------------------
# LEFT Y AXIS = PEDIGREE COMPLETENESS
# RIGHT Y AXIS = NUMBER OF ANIMALS
# --------------------------------------------------------------

scale_y_continuous(
  name = "Pedigree completeness (%)",
  
  limits = c(
    0,
    100
  ),
  
  breaks = seq(
    0,
    100,
    by = 10
  ),
  
  labels = function(x) {
    paste0(
      x,
      "%"
    )
  },
  
  sec.axis = sec_axis(
    ~ . / animal_scale,
    
    name = "Number of animals",
    
    labels = comma
  ),
  
  expand = expansion(
    mult = c(
      0,
      0
    )
  )
) +
  
  # --------------------------------------------------------------
# TITLES AND LEGEND
# --------------------------------------------------------------

labs(
  title =
    "Pedigree completeness and number of animals by year",
  
  colour =
    NULL
) +
  
  # --------------------------------------------------------------
# STYLE
# --------------------------------------------------------------

theme_classic(
  base_size = 13
) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    
    axis.title.x = element_text(
      face = "bold"
    ),
    
    axis.title.y.left = element_text(
      face = "bold"
    ),
    
    axis.title.y.right = element_text(
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    legend.position = "bottom",
    
    legend.direction = "horizontal",
    
    legend.box = "horizontal"
  )


# ================================================================
# DISPLAY GRAPH
# ================================================================

print(
  pedigree_plot
)


# ================================================================
# 5. SAVE GRAPH
# ================================================================

ggsave(
  "Pedigree_completeness_and_number_animals_by_year.png",
  pedigree_plot,
  width = 12,
  height = 7,
  dpi = 300
)
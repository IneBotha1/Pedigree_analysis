# ================================================================
# GRAPH THE NUMBER AND PERCENTAGE OF INDIVIDUALS BY LAP
# ================================================================

# Install once if necessary:
# install.packages(c("readxl", "ggplot2"))

library(readxl)
library(ggplot2)


# 1. READ THE EXCEL FILE ------------------------------------------

# The file is in the current working directory.
lap_data <- read_excel(
  "LAP_table.xlsx",
  sheet = 1
)


# Check what was imported.
print(lap_data)
names(lap_data)


# 2. STANDARDISE THE COLUMN NAMES ---------------------------------

# Retain the first three columns in case Excel contains empty columns.
lap_data <- lap_data[, 1:3]

names(lap_data) <- c(
  "LAP",
  "Number",
  "Percentage"
)


# Convert all three columns to numeric.
lap_data$LAP <- as.numeric(
  lap_data$LAP
)

lap_data$Number <- as.numeric(
  lap_data$Number
)

lap_data$Percentage <- as.numeric(
  lap_data$Percentage
)


# Remove completely empty rows.
lap_data <- lap_data[
  complete.cases(
    lap_data[, c(
      "LAP",
      "Number",
      "Percentage"
    )]
  ),
]


# If Excel stored percentages as proportions, for example 0.069
# instead of 6.9, convert them to percentage-point values.
if (
  max(
    lap_data$Percentage,
    na.rm = TRUE
  ) <= 1
) {
  lap_data$Percentage <-
    lap_data$Percentage * 100
}


# Order the rows by LAP.
lap_data <- lap_data[
  order(lap_data$LAP),
]


# 3. CHECK THE IMPORTED VALUES ------------------------------------

total_individuals <- sum(
  lap_data$Number,
  na.rm = TRUE
)


cat(
  "Total number of individuals:",
  total_individuals,
  "\n"
)

cat(
  "Sum of rounded percentages:",
  sum(
    lap_data$Percentage,
    na.rm = TRUE
  ),
  "%\n"
)


# Compare the supplied percentage with the percentage calculated
# directly from the number of individuals.
lap_data$Calculated_percentage <-
  100 *
  lap_data$Number /
  total_individuals


print(
  lap_data
)

# ================================================================
# SINGLE-LINE LAP GRAPH WITH TWO Y-AXES
# ================================================================

library(ggplot2)

# Total pedigree population represented in the table.
total_individuals <- sum(
  lap_data$Number,
  na.rm = TRUE
)

# Conversion factor:
# 1 percentage point corresponds to this number of individuals.
scale_factor <- total_individuals / 100


lap_line_plot <- ggplot(
  lap_data,
  aes(
    x = LAP,
    y = Number
  )
) +
  
  # Only one line is required because number and percentage
  # describe exactly the same distribution.
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 2.8
  ) +
  
  scale_x_continuous(
    breaks = sort(
      unique(lap_data$LAP)
    )
  ) +
  
  scale_y_continuous(
    name = "Number of individuals",
    
    labels = scales::comma,
    
    # Convert number of individuals to percentage on the right axis.
    sec.axis = sec_axis(
      trans = ~ . / scale_factor,
      
      name = "Percentage of population (%)",
      
      labels = function(x) {
        paste0(
          round(x, 1),
          "%"
        )
      }
    ),
    
    expand = expansion(
      mult = c(
        0,
        0.06
      )
    )
  ) +
  
  labs(
    title = "Distribution of longest ancestral path length",
    
    subtitle = paste0(
      "Full pedigree: ",
      format(
        total_individuals,
        big.mark = ","
      ),
      " individuals"
    ),
    
    x = "Longest ancestral path (generations)"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    axis.title.y.left = element_text(
      face = "bold"
    ),
    
    axis.title.y.right = element_text(
      face = "bold"
    )
  )


print(
  lap_line_plot
)

##Make the line Blue
lap_line_plot <- ggplot(
  lap_data,
  aes(
    x = LAP,
    y = Number
  )
) +
  
  geom_line(
    linewidth = 1,
    colour = "blue"
  ) +
  
  geom_point(
    size = 2.8,
    colour = "blue"
  ) +
  
  scale_x_continuous(
    breaks = sort(
      unique(lap_data$LAP)
    )
  ) +
  
  scale_y_continuous(
    name = "Number of individuals",
    
    labels = scales::comma,
    
    sec.axis = sec_axis(
      trans = ~ . / scale_factor,
      name = "Percentage of population (%)",
      labels = function(x) {
        paste0(round(x, 1), "%")
      }
    ),
    
    expand = expansion(
      mult = c(0, 0.06)
    )
  ) +
  
  labs(
    title = "Distribution of longest ancestral path length",
    
    subtitle = paste0(
      "Full pedigree: ",
      format(
        total_individuals,
        big.mark = ","
      ),
      " individuals"
    ),
    
    x = "Longest ancestral path (generations)"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    axis.title.y.left = element_text(
      face = "bold"
    ),
    
    axis.title.y.right = element_text(
      face = "bold"
    )
  )

print(lap_line_plot)

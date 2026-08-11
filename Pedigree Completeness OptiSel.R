#Pedigree Completeness OptiSel

library(optiSel)

# Read the complete, loop-resolved pedigree.
pedigree <- read.csv(
  "SAH_basic_ped_input_ungrouped_loopsres.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

# Prepare the three pedigree columns required by optiSel.
ped_input <- pedigree[, c(
  "ANI_ID",
  "SIRE_ID",
  "DAM_ID"
)]

names(ped_input) <- c(
  "Indiv",
  "Sire",
  "Dam"
)

# Prepare and validate the pedigree.
Pedig <- prePed(ped_input)

# Calculate individual pedigree-completeness measures.
pedigree_summary <- summary.Pedig(
  Pedig,
  maxd = 50,
  d = 4
)

summary(
  pedigree_summary[, c(
    "equiGen",
    "fullGen",
    "maxGen",
    "PCI",
    "Inbreeding"
  )]
)

c(
  Total_animals = nrow(pedigree_summary),
  
  Animals_with_ancestors =
    sum(pedigree_summary$maxGen > 0),
  
  Animals_without_ancestors =
    sum(pedigree_summary$maxGen == 0),
  
  Maximum_depth =
    max(pedigree_summary$maxGen),
  
  Maximum_equivalent_generations =
    max(pedigree_summary$equiGen),
  
  Mean_equivalent_generations =
    mean(pedigree_summary$equiGen),
  
  Mean_PCI =
    mean(pedigree_summary$PCI)
)

# ================================================================
# LINE GRAPH OF PEDIGREE COMPLETENESS BY GENERATION
# FULL PEDIGREE
# ================================================================


# Determine the deepest traced generation observed in the pedigree.
plot_maxd <- max(
  pedigree_summary$maxGen,
  na.rm = TRUE
)


cat(
  "Maximum observed traced generation:",
  plot_maxd,
  "\n"
)


# Calculate the proportion of known ancestors for every individual
# in every ancestral generation.
generation_completeness <- optiSel::completeness(
  Pedig,
  
  keep = as.character(
    pedigree_summary$Indiv
  ),
  
  maxd = plot_maxd,
  
  by = "Indiv"
)


# Check the returned columns.
head(
  generation_completeness
)

# Calculate mean completeness by ancestral generation.
mean_generation_completeness <- aggregate(
  Completeness ~ Generation,
  
  data = generation_completeness,
  
  FUN = mean,
  na.rm = TRUE
)


# Generation zero represents the individual itself and is always
# complete, so remove it from the graph.
mean_generation_completeness <-
  mean_generation_completeness[
    mean_generation_completeness$Generation >= 1,
  ]


# Express completeness as a percentage.
mean_generation_completeness$Completeness_percent <-
  100 *
  mean_generation_completeness$Completeness


print(
  mean_generation_completeness
)

# Create the line graph
pedigree_completeness_lineplot <- ggplot(
  mean_generation_completeness,
  
  aes(
    x = Generation,
    y = Completeness_percent
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 2.5
  ) +
  
  scale_x_continuous(
    breaks = seq(
      1,
      plot_maxd,
      by = 1
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
      by = 10
    ),
    
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  
  labs(
    x =
      "Ancestral generation",
    
    y =
      "Known ancestors (%)"
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
    )
  )


print(
  pedigree_completeness_lineplot
)

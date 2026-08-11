# SOUTH AFRICAN BOERPERD
# HISTORICAL MARGINAL GENETIC CONTRIBUTION ANALYSIS

# Reference populations:
#   Group 1 = animals born in 1973
#   Group 3 = animals born from 1974 to 1999
#   Group 4 = animals born from 2000 to 2024
#
# Historical pools:
#   S1 = 29 confirmed pre-1973 ancestors of the 1973 cohort
#   S2 = S1 plus 10 ancestors with unknown birth years
#   S3 = all 343 animals with known birth years before 1973

# 1. SET THE WORKING DIRECTORY ------------------------------------

setwd(
  "C:/2025/Paper drafts/Pedigree paper/1973 cohort analysis check"
)


# 2. LOAD PACKAGES ------------------------------------------------

install.packages(
  "rlang",
  repos = "https://cloud.r-project.org"
)

packageVersion("rlang")




# Run these lines once if the packages are not installed:
# install.packages("visPedigree")
install.packages("ggplot2")

library(visPedigree)
library(ggplot2)


# 3. SPECIFY THE THREE UNCUT CFC FILES ----------------------------

S1_filename <- "CFC_S1_crosscheck_input.txt"
S2_filename <- "CFC_S2_crosscheck_input.txt"
S3_filename <- "CFC_S3_crosscheck_input.txt"


# 4. FUNCTION TO READ A CFC FILE ---------------------------------

read_cfc <- function(filename) {
  
  # Read IDs as character values so that long IDs are not rounded.
  x <- read.table(
    filename,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    check.names = FALSE
  )
  
  # Remove accidental spaces from column headings.
  names(x) <- trimws(names(x))
  
  # Standardise the spelling of the group column.
  names(x)[tolower(names(x)) == "group"] <- "Group"
  
  # Confirm that all four required columns are present.
  required <- c(
    "Progeny",
    "Sire",
    "Dam",
    "Group"
  )
  
  if (!all(required %in% names(x))) {
    stop(
      filename,
      " is missing the following columns: ",
      paste(
        setdiff(required, names(x)),
        collapse = ", "
      )
    )
  }
  
  # Keep only the required columns.
  x <- x[, required]
  
  # Remove spaces from the IDs.
  x$Progeny <- trimws(x$Progeny)
  x$Sire <- trimws(x$Sire)
  x$Dam <- trimws(x$Dam)
  
  # Convert unknown-parent codes to missing values.
  unknown_parent <- c(
    "",
    "0",
    "-9",
    "*",
    "NA",
    "N/A",
    "UNKNOWN"
  )
  
  x$Sire[x$Sire %in% unknown_parent] <- NA_character_
  x$Dam[x$Dam %in% unknown_parent] <- NA_character_
  
  # Convert the CFC group to an integer.
  x$Group <- as.integer(x$Group)
  
  # Stop if invalid group values were found.
  if (anyNA(x$Group)) {
    stop(
      "Missing or invalid group values were found in ",
      filename
    )
  }
  
  # Each animal should occur only once in the pedigree.
  if (anyDuplicated(x$Progeny) > 0) {
    stop(
      "Duplicated Progeny IDs were found in ",
      filename
    )
  }
  
  return(x)
}


# 5. READ THE THREE FILES -----------------------------------------

S1_file <- read_cfc(S1_filename)
S2_file <- read_cfc(S2_filename)
S3_file <- read_cfc(S3_filename)


# 6. CONFIRM THAT THE PEDIGREE IS IDENTICAL IN ALL FILES ----------

# The group assignments may differ, but Progeny, Sire and Dam
# should be identical in all three scenario files.

get_pedigree_core <- function(x) {
  
  core <- x[
    order(x$Progeny),
    c("Progeny", "Sire", "Dam")
  ]
  
  rownames(core) <- NULL
  
  return(core)
}


if (
  !identical(
    get_pedigree_core(S1_file),
    get_pedigree_core(S2_file)
  ) ||
  !identical(
    get_pedigree_core(S1_file),
    get_pedigree_core(S3_file)
  )
) {
  stop(
    "The pedigree structure differs among the three files. ",
    "Use files created from the same complete, uncut pedigree."
  )
}


# 7. EXTRACT THE HISTORICAL-POOL IDS ------------------------------

# Group 2 contains the historical pool in each scenario.
# Use S1_ids, S2_ids and S3_ids rather than S1, S2 and S3
# so these vectors cannot be confused with complete data frames.

S1_ids <- unique(
  S1_file$Progeny[S1_file$Group == 2]
)

S2_ids <- unique(
  S2_file$Progeny[S2_file$Group == 2]
)

S3_ids <- unique(
  S3_file$Progeny[S3_file$Group == 2]
)


# Check the historical-pool sizes.

pool_check <- data.frame(
  Scenario = c("S1", "S2", "S3"),
  Expected = c(29, 39, 343),
  Observed = c(
    length(S1_ids),
    length(S2_ids),
    length(S3_ids)
  )
)

print(pool_check)


# Stop if the historical pools are incorrect.

stopifnot(
  length(S1_ids) == 29,
  length(S2_ids) == 39,
  length(S3_ids) == 343
)


# Confirm that all S1 animals also occur in S2 and S3.

stopifnot(
  all(S1_ids %in% S2_ids),
  all(S1_ids %in% S3_ids)
)


# Show how many animals are added by S2 and S3.

cat(
  "\nAnimals added in S2:",
  length(setdiff(S2_ids, S1_ids)),
  "\n"
)

cat(
  "Additional animals in S3:",
  length(setdiff(S3_ids, S1_ids)),
  "\n"
)


# 8. DEFINE THE REFERENCE POPULATIONS -----------------------------

# The reference populations are the same in all three scenarios,
# so their IDs can be obtained from the S1 file.

ref_1973_ids <- unique(
  S1_file$Progeny[S1_file$Group == 1]
)

ref_1974_1999_ids <- unique(
  S1_file$Progeny[S1_file$Group == 3]
)

ref_2000_2024_ids <- unique(
  S1_file$Progeny[S1_file$Group == 4]
)


# Check the reference-population sizes.

reference_check <- data.frame(
  Reference_population = c(
    "1973 cohort",
    "1974-1999",
    "2000-2024"
  ),
  Expected = c(
    48,
    3114,
    12216
  ),
  Observed = c(
    length(ref_1973_ids),
    length(ref_1974_1999_ids),
    length(ref_2000_2024_ids)
  )
)

print(reference_check)


# Stop if the reference populations are incorrect.

stopifnot(
  length(ref_1973_ids) == 48,
  length(ref_1974_1999_ids) == 3114,
  length(ref_2000_2024_ids) == 12216
)


# 9. CREATE ONE COMPLETE PEDIGREE FOR visPedigree ----------------

# Only one pedigree is needed because the pedigree structure is
# identical in S1, S2 and S3.

ped_input <- S1_file[, c(
  "Progeny",
  "Sire",
  "Dam"
)]


# tidyped() expects the first three columns to represent:
# individual, sire and dam.

names(ped_input) <- c(
  "Ind",
  "Sire",
  "Dam"
)


# Create the tidyped object from the complete, uncut pedigree.
# The historical animals retain their original parent records.

ped_tidy <- tidyped(
  ped_input,
  addgen = TRUE,
  addnum = TRUE
)


# Confirm that the conversion worked.

class(ped_tidy)

head(ped_tidy)

cat(
  "\nOriginal pedigree records:",
  nrow(ped_input),
  "\n"
)

cat(
  "Prepared pedigree records:",
  nrow(ped_tidy),
  "\n"
)


# 10. CONFIRM THAT ALL REFERENCE ANIMALS ARE PRESENT --------------

tidy_ids <- as.character(ped_tidy$Ind)

stopifnot(
  all(ref_1973_ids %in% tidy_ids),
  all(ref_1974_1999_ids %in% tidy_ids),
  all(ref_2000_2024_ids %in% tidy_ids)
)


# 11. RUN THE MARGINAL ANCESTOR ANALYSIS --------------------------

# mode = "ancestor" requests marginal ancestor contributions.
#
# top = nrow(ped_tidy) requests all contributing ancestors,
# rather than only the default top contributors.


# Reference population 1: 1973 cohort.

contrib_1973 <- pedcontrib(
  ped = ped_tidy,
  reference = ref_1973_ids,
  mode = "ancestor",
  top = nrow(ped_tidy)
)


# Reference population 2: animals born from 1974 to 1999.

contrib_1974_1999 <- pedcontrib(
  ped = ped_tidy,
  reference = ref_1974_1999_ids,
  mode = "ancestor",
  top = nrow(ped_tidy)
)


# Reference population 3: animals born from 2000 to 2024.

contrib_2000_2024 <- pedcontrib(
  ped = ped_tidy,
  reference = ref_2000_2024_ids,
  mode = "ancestor",
  top = nrow(ped_tidy)
)


# 12. EXAMINE THE SUMMARY OF EACH ANALYSIS ------------------------

# n_ref = number of animals in the reference population
# n_ancestor = number of contributing ancestors
# f_a = effective number of ancestors
#
# f_a is useful supporting information, but it is not the summed
# historical-pool contribution used for the final graph.

contrib_1973$summary

contrib_1974_1999$summary

contrib_2000_2024$summary


# 13. EXTRACT THE ANCESTOR TABLES ---------------------------------

# Each ancestor table contains:
#
# Ind        = ancestor ID
# Contrib    = marginal genetic contribution
# CumContrib = cumulative contribution
# Rank       = contribution rank

anc_1973 <- as.data.frame(
  contrib_1973$ancestors
)

anc_1974_1999 <- as.data.frame(
  contrib_1974_1999$ancestors
)

anc_2000_2024 <- as.data.frame(
  contrib_2000_2024$ancestors
)


# Ensure the IDs and contributions use the correct data types.

anc_1973$Ind <- as.character(anc_1973$Ind)
anc_1973$Contrib <- as.numeric(anc_1973$Contrib)

anc_1974_1999$Ind <- as.character(
  anc_1974_1999$Ind
)
anc_1974_1999$Contrib <- as.numeric(
  anc_1974_1999$Contrib
)

anc_2000_2024$Ind <- as.character(
  anc_2000_2024$Ind
)
anc_2000_2024$Contrib <- as.numeric(
  anc_2000_2024$Contrib
)


# Open these tables manually when needed.

# View(anc_1973)
# View(anc_1974_1999)
# View(anc_2000_2024)


# 14. CHECK THE TOTAL CONTRIBUTIONS -------------------------------

# Because all contributing ancestors were requested, the complete
# marginal-contribution table should sum to approximately 1.

contribution_total_check <- data.frame(
  Reference_population = c(
    "1973 cohort",
    "1974-1999",
    "2000-2024"
  ),
  Total_contribution = c(
    sum(anc_1973$Contrib, na.rm = TRUE),
    sum(anc_1974_1999$Contrib, na.rm = TRUE),
    sum(anc_2000_2024$Contrib, na.rm = TRUE)
  )
)

print(contribution_total_check)


# Warn if a table does not sum to approximately 1.

if (
  any(
    abs(
      contribution_total_check$Total_contribution - 1
    ) > 0.000001
  )
) {
  warning(
    "At least one contribution table does not sum to approximately 1."
  )
}


# 15. FUNCTION TO SUM ONE HISTORICAL POOL -------------------------

summarise_historical_pool <- function(
    ancestor_table,
    historical_ids,
    scenario_name,
    reference_name
) {
  
  # Select ancestors that belong to the specified historical pool.
  matched <- ancestor_table[
    ancestor_table$Ind %in% historical_ids,
    ,
    drop = FALSE
  ]
  
  # Sum their marginal genetic contributions.
  contribution <- sum(
    matched$Contrib,
    na.rm = TRUE
  )
  
  # Create one summary row.
  summary <- data.frame(
    Reference_population = reference_name,
    Scenario = scenario_name,
    Historical_pool_size = length(historical_ids),
    Historical_ancestors_matched = nrow(matched),
    Historical_ancestors_not_matched =
      length(historical_ids) - nrow(matched),
    Contribution_proportion = contribution,
    Contribution_percent = 100 * contribution
  )
  
  # Add labels to the detailed matched ancestor table.
  matched$Reference_population <- reference_name
  matched$Scenario <- scenario_name
  
  return(
    list(
      summary = summary,
      matched = matched
    )
  )
}


# 16. STORE POOLS AND REFERENCE TABLES IN LISTS -------------------

historical_pools <- list(
  S1 = S1_ids,
  S2 = S2_ids,
  S3 = S3_ids
)

ancestor_tables <- list(
  "1973 cohort" = anc_1973,
  "1974-1999" = anc_1974_1999,
  "2000-2024" = anc_2000_2024
)


# 17. CALCULATE ALL NINE RESULTS ---------------------------------

summary_list <- list()
detail_list <- list()

result_number <- 1


for (reference_name in names(ancestor_tables)) {
  
  for (scenario_name in names(historical_pools)) {
    
    current_result <- summarise_historical_pool(
      ancestor_table =
        ancestor_tables[[reference_name]],
      
      historical_ids =
        historical_pools[[scenario_name]],
      
      scenario_name =
        scenario_name,
      
      reference_name =
        reference_name
    )
    
    summary_list[[result_number]] <-
      current_result$summary
    
    detail_list[[result_number]] <-
      current_result$matched
    
    result_number <- result_number + 1
  }
}


# Combine the nine summary rows.

results <- do.call(
  rbind,
  summary_list
)

rownames(results) <- NULL


# Combine the detailed matched ancestor records.

matched_ancestor_details <- do.call(
  rbind,
  detail_list
)

rownames(matched_ancestor_details) <- NULL


# 18. ORDER AND DISPLAY THE FINAL RESULTS -------------------------

results$Reference_population <- factor(
  results$Reference_population,
  levels = c(
    "1973 cohort",
    "1974-1999",
    "2000-2024"
  )
)

results$Scenario <- factor(
  results$Scenario,
  levels = c(
    "S1",
    "S2",
    "S3"
  )
)

results <- results[
  order(
    results$Reference_population,
    results$Scenario
  ),
]

rownames(results) <- NULL


# This is the main results table.

print(results)


# Display only the variables most important for interpretation.

print(
  results[, c(
    "Reference_population",
    "Scenario",
    "Historical_pool_size",
    "Historical_ancestors_matched",
    "Contribution_proportion",
    "Contribution_percent"
  )]
)


# 19. CREATE A SIMPLE WIDE PERCENTAGE TABLE -----------------------

percentage_table <- reshape(
  results[, c(
    "Reference_population",
    "Scenario",
    "Contribution_percent"
  )],
  idvar = "Reference_population",
  timevar = "Scenario",
  direction = "wide"
)

names(percentage_table) <- sub(
  "Contribution_percent\\.",
  "",
  names(percentage_table)
)

print(percentage_table)


# 20. CREATE AN OUTPUT FOLDER -------------------------------------

if (!dir.exists("marginal_contribution_results")) {
  dir.create("marginal_contribution_results")
}


# 21. SAVE THE RESULTS TABLES -------------------------------------

write.csv(
  pool_check,
  paste0(
    "marginal_contribution_results/",
    "historical_pool_check.csv"
  ),
  row.names = FALSE
)

write.csv(
  reference_check,
  paste0(
    "marginal_contribution_results/",
    "reference_population_check.csv"
  ),
  row.names = FALSE
)

write.csv(
  contribution_total_check,
  paste0(
    "marginal_contribution_results/",
    "contribution_total_check.csv"
  ),
  row.names = FALSE
)

write.csv(
  results,
  paste0(
    "marginal_contribution_results/",
    "marginal_contribution_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  percentage_table,
  paste0(
    "marginal_contribution_results/",
    "marginal_contribution_percentage_table.csv"
  ),
  row.names = FALSE
)

write.csv(
  matched_ancestor_details,
  paste0(
    "marginal_contribution_results/",
    "matched_historical_ancestor_details.csv"
  ),
  row.names = FALSE
)


# Save the full ancestor-contribution tables.

write.csv(
  anc_1973,
  paste0(
    "marginal_contribution_results/",
    "all_ancestors_1973.csv"
  ),
  row.names = FALSE
)

write.csv(
  anc_1974_1999,
  paste0(
    "marginal_contribution_results/",
    "all_ancestors_1974_1999.csv"
  ),
  row.names = FALSE
)

write.csv(
  anc_2000_2024,
  paste0(
    "marginal_contribution_results/",
    "all_ancestors_2000_2024.csv"
  ),
  row.names = FALSE
)


# 22. CREATE THE FINAL GRAPH --------------------------------------

contribution_plot <- ggplot(
  results,
  aes(
    x = Reference_population,
    y = Contribution_percent,
    fill = Scenario
  )
) +
  
  # Draw the three scenarios beside each other.
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "black"
  ) +
  
  # Add the percentage above each bar.
  geom_text(
    aes(
      label = paste0(
        round(Contribution_percent, 1),
        "%"
      )
    ),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 3
  ) +
  
  # Use specific colours (obtain codes)
  scale_fill_manual(
    values = c(
      "S1" = "#67BFA5",
      "S2" = "#FC8D62",
      "S3" = "#8DA0CB"
    ),
    labels = c(
      "S1: Known pre-1973 ancestors",
      "S2: Ancestors including unknown years",
      "S3: All known pre-1973 horses"
    )
  ) +
  
  # Leave room above the tallest bar.
  scale_y_continuous(
    limits = c(
      0,
      max(
        results$Contribution_percent,
        na.rm = TRUE
      ) * 1.15
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  labs(
    title =
      "Persistence of historical genetic contribution",
    
    subtitle =
      paste(
        "Marginal contributions across the selected",
        "South African Boerperd reference populations"
      ),
    
    x = "Reference population",
    
    y = "Marginal genetic contribution (%)",
    
    fill = "Historical pool definition"
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
      size = 12
    ),
    
    axis.text.x = element_text(
      angle = 15,
      hjust = 1
    ),
    
    legend.position = "right"
  )


#Displaggplot2# Display the graph.

print(contribution_plot)


# 23. SAVE THE GRAPH ----------------------------------------------

ggsave(
  filename = paste0(
    "marginal_contribution_results/",
    "historical_marginal_contribution.png"
  ),
  plot = contribution_plot,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  filename = paste0(
    "marginal_contribution_results/",
    "historical_marginal_contribution.pdf"
  ),
  plot = contribution_plot,
  width = 11,
  height = 7
)


# 24. SAVE PACKAGE AND R VERSION INFORMATION ----------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  paste0(
    "marginal_contribution_results/",
    "sessionInfo.txt"
  )
)


# 25. FINAL MESSAGE -----------------------------------------------

cat(
  "\nAnalysis complete.\n",
  "Results were saved in the folder:\n",
  "marginal_contribution_results\n",
  sep = ""
)

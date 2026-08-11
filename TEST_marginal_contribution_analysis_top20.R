# ================================================================
# SOUTH AFRICAN BOERPERD
# HISTORICAL MARGINAL GENETIC CONTRIBUTION ANALYSIS
# SUPERVISOR REFERENCE POPULATIONS — TOP-20 REPLICATION TEST
# ================================================================
#
# This script skips the CFC workflow and works directly from the
# complete, loop-resolved and uncut pedigree.
#
# Input file:
#   SAH_basic_ped_input_ungrouped_loopsres.csv
#
# Historical pools:
#   S1 = 29 confirmed ancestors of the 1973 cohort with known
#        birth years before 1973
#   S2 = S1 plus 10 ancestors of the 1973 cohort with unknown
#        birth years, giving 39 animals
#   S3 = all 343 animals in the complete pedigree with known
#        birth years before 1973
#
# Supervisor reference populations:
#   1. 1973-1989
#   2. 1990-2004
#   3. 2005-2014
#   4. 2015-2024
#
# IMPORTANT:
#   Animals born in 2025 or later are explicitly excluded.
#
# This version uses top = 20 in pedcontrib() as a replication test.
# If the supervisor confirms that all contributing ancestors were
# returned, change TOP_N from 20L to nrow(ped_tidy).
# ================================================================

# 1. SET THE WORKING DIRECTORY ------------------------------------

setwd(
  "C:/2025/Paper drafts/Pedigree paper/1973 cohort analysis check"
)

# 2. LOAD PACKAGES ------------------------------------------------

# Run only if the packages are not already installed:
# install.packages("visPedigree")
# install.packages("ggplot2")

library(visPedigree)
library(ggplot2)

# 3. SET THE NUMBER OF ANCESTORS RETURNED -------------------------

TOP_N <- 20L

# 4. READ THE COMPLETE LOOP-RESOLVED PEDIGREE ---------------------

pedigree <- read.csv(
  "SAH_basic_ped_input_ungrouped_loopsres.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

# 5. CONFIRM THAT REQUIRED COLUMNS ARE PRESENT --------------------

required_columns <- c(
  "ANI_ID",
  "SIRE_ID",
  "DAM_ID",
  "BIRTH_DTM"
)

if (!all(required_columns %in% names(pedigree))) {
  stop(
    "Missing required columns: ",
    paste(
      setdiff(required_columns, names(pedigree)),
      collapse = ", "
    )
  )
}

# 6. STANDARDISE ANIMAL AND PARENT IDS ----------------------------

pedigree$ANI_ID <- trimws(as.character(pedigree$ANI_ID))
pedigree$SIRE_ID <- trimws(as.character(pedigree$SIRE_ID))
pedigree$DAM_ID <- trimws(as.character(pedigree$DAM_ID))

unknown_parent_codes <- c(
  "",
  "0",
  "-9",
  "*",
  "NA",
  "N/A",
  "UNKNOWN"
)

pedigree$SIRE_ID[
  pedigree$SIRE_ID %in% unknown_parent_codes
] <- NA_character_

pedigree$DAM_ID[
  pedigree$DAM_ID %in% unknown_parent_codes
] <- NA_character_

if (
  any(is.na(pedigree$ANI_ID)) ||
  any(pedigree$ANI_ID == "")
) {
  stop("Missing or blank ANI_ID values were found.")
}

if (anyDuplicated(pedigree$ANI_ID) > 0) {
  duplicated_ids <- unique(
    pedigree$ANI_ID[
      duplicated(pedigree$ANI_ID)
    ]
  )
  
  stop(
    "Duplicated ANI_ID values were found. Examples: ",
    paste(head(duplicated_ids, 10), collapse = ", ")
  )
}

# 7. PARSE BIRTH DATES AND EXTRACT BIRTH YEARS --------------------

birth_date <- as.Date(
  pedigree$BIRTH_DTM,
  format = "%Y/%m/%d"
)

needs_second_format <- (
  is.na(birth_date) &
    !is.na(pedigree$BIRTH_DTM) &
    trimws(pedigree$BIRTH_DTM) != ""
)

birth_date[needs_second_format] <- as.Date(
  pedigree$BIRTH_DTM[needs_second_format],
  format = "%Y-%m-%d"
)

birth_year <- as.integer(
  format(birth_date, "%Y")
)

pedigree$BIRTH_YEAR_ANALYSIS <- birth_year

cat(
  "\nAnimals with a known birth year:",
  sum(!is.na(birth_year)),
  "\n"
)

cat(
  "Animals with an unknown birth year:",
  sum(is.na(birth_year)),
  "\n"
)

cat(
  "Earliest known birth year:",
  min(birth_year, na.rm = TRUE),
  "\n"
)

cat(
  "Latest known birth year in the pedigree:",
  max(birth_year, na.rm = TRUE),
  "\n"
)

# Show the number of animals assigned to each recent birth year.

recent_year_counts <- as.data.frame(
  table(
    birth_year[
      !is.na(birth_year) &
        birth_year >= 2010
    ]
  )
)

names(recent_year_counts) <- c(
  "Birth_year",
  "Number_of_animals"
)

recent_year_counts$Birth_year <- as.integer(
  as.character(recent_year_counts$Birth_year)
)

print(recent_year_counts)

failed_date_rows <- pedigree[
  is.na(birth_year) &
    !is.na(pedigree$BIRTH_DTM) &
    trimws(pedigree$BIRTH_DTM) != "",
  c(
    "ANI_ID",
    "BIRTH_DTM"
  ),
  drop = FALSE
]

cat(
  "\nNumber of nonblank dates that failed to parse:",
  nrow(failed_date_rows),
  "\n"
)

print(
  head(
    failed_date_rows,
    30
  )
)

animals_2024 <- pedigree[
  !is.na(birth_year) &
    birth_year == 2024,
  c(
    "ANI_ID",
    "SIRE_ID",
    "DAM_ID",
    "BIRTH_DTM"
  ),
  drop = FALSE
]

animals_2025 <- pedigree[
  !is.na(birth_year) &
    birth_year == 2025,
  c(
    "ANI_ID",
    "SIRE_ID",
    "DAM_ID",
    "BIRTH_DTM"
  ),
  drop = FALSE
]

cat(
  "\nAnimals born in 2024:",
  nrow(animals_2024),
  "\n"
)

cat(
  "Animals born in 2025:",
  nrow(animals_2025),
  "\n"
)

print(animals_2024)
print(animals_2025)

# 8. DEFINE THE 1973 COHORT ---------------------------------------

cohort_1973_ids <- unique(
  pedigree$ANI_ID[
    !is.na(birth_year) &
      birth_year == 1973
  ]
)

cat(
  "\nNumber of animals born in 1973:",
  length(cohort_1973_ids),
  "\n"
)

stopifnot(
  length(cohort_1973_ids) == 48
)

# 9. TRACE ALL RECORDED ANCESTORS OF THE 1973 COHORT --------------

get_ancestors <- function(start_ids, ped) {
  
  current_ids <- unique(as.character(start_ids))
  examined_ids <- current_ids
  ancestor_ids <- character()
  
  repeat {
    
    current_rows <- ped[
      ped$ANI_ID %in% current_ids,
      c("SIRE_ID", "DAM_ID"),
      drop = FALSE
    ]
    
    parent_ids <- unique(
      c(
        current_rows$SIRE_ID,
        current_rows$DAM_ID
      )
    )
    
    parent_ids <- parent_ids[
      !is.na(parent_ids) &
        parent_ids != ""
    ]
    
    new_parent_ids <- setdiff(
      parent_ids,
      examined_ids
    )
    
    if (length(new_parent_ids) == 0) {
      break
    }
    
    ancestor_ids <- union(
      ancestor_ids,
      new_parent_ids
    )
    
    examined_ids <- union(
      examined_ids,
      new_parent_ids
    )
    
    current_ids <- new_parent_ids
  }
  
  return(
    intersect(
      ancestor_ids,
      ped$ANI_ID
    )
  )
}

ancestors_1973_ids <- get_ancestors(
  start_ids = cohort_1973_ids,
  ped = pedigree
)

cat(
  "Number of recorded ancestors of the 1973 cohort:",
  length(ancestors_1973_ids),
  "\n"
)

# 10. DEFINE THE THREE HISTORICAL POOLS ---------------------------

S1_ids <- unique(
  pedigree$ANI_ID[
    pedigree$ANI_ID %in% ancestors_1973_ids &
      !is.na(birth_year) &
      birth_year < 1973
  ]
)

S2_ids <- unique(
  pedigree$ANI_ID[
    pedigree$ANI_ID %in% ancestors_1973_ids &
      (
        is.na(birth_year) |
          birth_year < 1973
      )
  ]
)

S3_ids <- unique(
  pedigree$ANI_ID[
    !is.na(birth_year) &
      birth_year < 1973
  ]
)

historical_pools <- list(
  S1 = S1_ids,
  S2 = S2_ids,
  S3 = S3_ids
)

historical_pool_check <- data.frame(
  Scenario = c("S1", "S2", "S3"),
  Expected_size = c(29, 39, 343),
  Observed_size = c(
    length(S1_ids),
    length(S2_ids),
    length(S3_ids)
  )
)

cat("\nHistorical-pool size check:\n")
print(historical_pool_check)

stopifnot(
  length(S1_ids) == 29,
  length(S2_ids) == 39,
  length(S3_ids) == 343
)

stopifnot(
  all(S1_ids %in% S2_ids),
  all(S1_ids %in% S3_ids)
)

cat(
  "\nAnimals added to S2:",
  length(setdiff(S2_ids, S1_ids)),
  "\n"
)

cat(
  "Additional animals in S3 relative to S1:",
  length(setdiff(S3_ids, S1_ids)),
  "\n"
)

# 11. DEFINE THE FOUR SUPERVISOR REFERENCE POPULATIONS ------------

ref_1973_1989_ids <- unique(
  pedigree$ANI_ID[
    !is.na(birth_year) &
      birth_year >= 1973 &
      birth_year <= 1989
  ]
)

ref_1990_2004_ids <- unique(
  pedigree$ANI_ID[
    !is.na(birth_year) &
      birth_year >= 1990 &
      birth_year <= 2004
  ]
)

ref_2005_2014_ids <- unique(
  pedigree$ANI_ID[
    !is.na(birth_year) &
      birth_year >= 2005 &
      birth_year <= 2014
  ]
)

# Explicitly restricted to 2015-2024.
# Animals born in 2025 or later are excluded.

ref_2015_2024_ids <- unique(
  pedigree$ANI_ID[
    !is.na(birth_year) &
      birth_year >= 2015 &
      birth_year <= 2024
  ]
)

reference_populations <- list(
  "1973-1989" = ref_1973_1989_ids,
  "1990-2004" = ref_1990_2004_ids,
  "2005-2014" = ref_2005_2014_ids,
  "2015-2024" = ref_2015_2024_ids
)

reference_population_check <- data.frame(
  Reference_population = names(reference_populations),
  Start_year = c(1973, 1990, 2005, 2015),
  End_year = c(1989, 2004, 2014, 2024),
  Number_of_animals = lengths(reference_populations)
)

cat("\nSupervisor reference-population sizes:\n")
print(reference_population_check)

if (any(reference_population_check$Number_of_animals == 0)) {
  stop(
    "At least one reference population contains no animals."
  )
}

all_reference_ids <- unlist(
  reference_populations,
  use.names = FALSE
)

if (anyDuplicated(all_reference_ids) > 0) {
  stop(
    "The four reference populations overlap unexpectedly."
  )
}

post_2024_ids <- pedigree$ANI_ID[
  !is.na(birth_year) &
    birth_year > 2024
]

stopifnot(
  length(
    intersect(
      ref_2015_2024_ids,
      post_2024_ids
    )
  ) == 0
)

# 12. SAVE HISTORICAL AND REFERENCE ID LISTS ----------------------

write.csv(
  data.frame(ANI_ID = S1_ids),
  "historical_pool_S1_IDs_supervisor_periods.csv",
  row.names = FALSE
)

write.csv(
  data.frame(ANI_ID = S2_ids),
  "historical_pool_S2_IDs_supervisor_periods.csv",
  row.names = FALSE
)

write.csv(
  data.frame(ANI_ID = S3_ids),
  "historical_pool_S3_IDs_supervisor_periods.csv",
  row.names = FALSE
)

write.csv(
  data.frame(ANI_ID = ref_1973_1989_ids),
  "reference_IDs_1973_1989.csv",
  row.names = FALSE
)

write.csv(
  data.frame(ANI_ID = ref_1990_2004_ids),
  "reference_IDs_1990_2004.csv",
  row.names = FALSE
)

write.csv(
  data.frame(ANI_ID = ref_2005_2014_ids),
  "reference_IDs_2005_2014.csv",
  row.names = FALSE
)

write.csv(
  data.frame(ANI_ID = ref_2015_2024_ids),
  "reference_IDs_2015_2024.csv",
  row.names = FALSE
)

# 13. PREPARE THE COMPLETE PEDIGREE FOR visPedigree ---------------

ped_input <- pedigree[, c(
  "ANI_ID",
  "SIRE_ID",
  "DAM_ID"
)]

names(ped_input) <- c(
  "Ind",
  "Sire",
  "Dam"
)

ped_tidy <- tidyped(
  ped_input,
  addgen = TRUE,
  addnum = TRUE
)

cat("\nClass of prepared pedigree:\n")
print(class(ped_tidy))

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

# 14. CONFIRM THAT ALL REFERENCE ANIMALS ARE PRESENT --------------

tidy_ids <- as.character(
  ped_tidy$Ind
)

missing_reference_ids <- lapply(
  reference_populations,
  function(ids) {
    setdiff(ids, tidy_ids)
  }
)

missing_reference_counts <- vapply(
  missing_reference_ids,
  length,
  integer(1)
)

if (any(missing_reference_counts > 0)) {
  print(missing_reference_counts)
  
  stop(
    "Some reference animals are missing from ped_tidy."
  )
}

# 15. RUN pedcontrib FOR ALL FOUR REFERENCE POPULATIONS -----------

contribution_objects <- lapply(
  reference_populations,
  function(reference_ids) {
    
    pedcontrib(
      ped = ped_tidy,
      reference = reference_ids,
      mode = "ancestor",
      top = TOP_N
    )
  }
)

names(contribution_objects) <- names(
  reference_populations
)

# 16. DISPLAY THE FOUR ANALYSIS SUMMARIES -------------------------

for (reference_name in names(contribution_objects)) {
  
  cat(
    "\n",
    reference_name,
    " contribution summary:\n",
    sep = ""
  )
  
  print(
    contribution_objects[[reference_name]]$summary
  )
}

# 17. EXTRACT AND STANDARDISE THE ANCESTOR TABLES -----------------

prepare_ancestor_table <- function(contribution_object) {
  
  ancestor_table <- as.data.frame(
    contribution_object$ancestors
  )
  
  ancestor_table$Ind <- as.character(
    ancestor_table$Ind
  )
  
  ancestor_table$Contrib <- as.numeric(
    ancestor_table$Contrib
  )
  
  return(ancestor_table)
}

ancestor_tables <- lapply(
  contribution_objects,
  prepare_ancestor_table
)

# 18. CHECK THE CONTRIBUTION REPRESENTED BY THE TOP 20 ------------

returned_top20_contribution_check <- data.frame(
  Reference_population = names(ancestor_tables),
  
  Ancestors_returned = vapply(
    ancestor_tables,
    nrow,
    integer(1)
  ),
  
  Contribution_represented = vapply(
    ancestor_tables,
    function(x) {
      sum(x$Contrib, na.rm = TRUE)
    },
    numeric(1)
  )
)

returned_top20_contribution_check$Contribution_percent_represented <-
  100 *
  returned_top20_contribution_check$Contribution_represented

cat(
  "\nContribution represented by the returned top 20 ancestors:\n"
)

print(
  returned_top20_contribution_check
)

# 19. FUNCTION TO SUM ONE HISTORICAL POOL -------------------------

summarise_historical_pool <- function(
    ancestor_table,
    historical_ids,
    scenario_name,
    reference_name
) {
  
  matched <- ancestor_table[
    ancestor_table$Ind %in% historical_ids,
    ,
    drop = FALSE
  ]
  
  contribution <- sum(
    matched$Contrib,
    na.rm = TRUE
  )
  
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
  
  matched$Reference_population <- rep(
    reference_name,
    nrow(matched)
  )
  
  matched$Scenario <- rep(
    scenario_name,
    nrow(matched)
  )
  
  return(
    list(
      summary = summary,
      matched = matched
    )
  )
}

# 20. CALCULATE ALL TWELVE RESULTS --------------------------------

summary_list <- list()
detail_list <- list()
result_number <- 1L

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
    
    result_number <- result_number + 1L
  }
}

results <- do.call(
  rbind,
  summary_list
)

rownames(results) <- NULL

matched_ancestor_details <- do.call(
  rbind,
  detail_list
)

rownames(matched_ancestor_details) <- NULL

# 21. ORDER AND DISPLAY THE RESULTS -------------------------------

reference_order <- c(
  "1973-1989",
  "1990-2004",
  "2005-2014",
  "2015-2024"
)

results$Reference_population <- factor(
  results$Reference_population,
  levels = reference_order
)

results$Scenario <- factor(
  results$Scenario,
  levels = c("S1", "S2", "S3")
)

results <- results[
  order(
    results$Reference_population,
    results$Scenario
  ),
]

rownames(results) <- NULL

cat("\nComplete final results:\n")
print(results)

cat("\nMain values for interpretation:\n")

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

# 22. CREATE A WIDE PERCENTAGE TABLE ------------------------------

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

cat("\nPercentage table:\n")
print(percentage_table)

# 23. CREATE THE EFFECTIVE-ANCESTOR SUMMARY -----------------------

effective_ancestor_summary <- do.call(
  rbind,
  lapply(
    names(contribution_objects),
    function(reference_name) {
      
      current_summary <-
        contribution_objects[[reference_name]]$summary
      
      data.frame(
        Reference_population = reference_name,
        
        Number_of_reference_animals =
          current_summary$n_ref,
        
        Number_of_contributing_ancestors =
          current_summary$n_ancestor,
        
        Number_of_ancestors_returned =
          current_summary$n_ancestor_show,
        
        Effective_number_of_ancestors_fa =
          current_summary$f_a,
        
        Shannon_effective_ancestors_faH =
          current_summary$f_a_H
      )
    }
  )
)

rownames(effective_ancestor_summary) <- NULL

cat("\nEffective-ancestor summary:\n")
print(effective_ancestor_summary)

# 24. CREATE A SEPARATE OUTPUT FOLDER -----------------------------

output_folder <-
  "marginal_contribution_results_supervisor_periods_TOP20_2015_2024"

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# 25. SAVE THE OUTPUT TABLES --------------------------------------

write.csv(
  historical_pool_check,
  file.path(
    output_folder,
    "historical_pool_check.csv"
  ),
  row.names = FALSE
)

write.csv(
  reference_population_check,
  file.path(
    output_folder,
    "reference_population_check.csv"
  ),
  row.names = FALSE
)

write.csv(
  returned_top20_contribution_check,
  file.path(
    output_folder,
    "returned_top20_contribution_check.csv"
  ),
  row.names = FALSE
)

write.csv(
  results,
  file.path(
    output_folder,
    "marginal_contribution_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  percentage_table,
  file.path(
    output_folder,
    "marginal_contribution_percentage_table.csv"
  ),
  row.names = FALSE
)

write.csv(
  matched_ancestor_details,
  file.path(
    output_folder,
    "matched_historical_ancestor_details.csv"
  ),
  row.names = FALSE
)

write.csv(
  effective_ancestor_summary,
  file.path(
    output_folder,
    "effective_ancestor_summary.csv"
  ),
  row.names = FALSE
)

for (reference_name in names(ancestor_tables)) {
  
  safe_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    reference_name
  )
  
  write.csv(
    ancestor_tables[[reference_name]],
    file.path(
      output_folder,
      paste0(
        "top20_ancestors_",
        safe_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

saveRDS(
  contribution_objects,
  file.path(
    output_folder,
    "pedcontrib_objects_TOP20.rds"
  )
)

# 26. CREATE THE GROUPED BAR GRAPH --------------------------------

contribution_plot <- ggplot(
  results,
  aes(
    x = Reference_population,
    y = Contribution_percent,
    fill = Scenario
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.72,
    colour = "black"
  ) +
  
  geom_text(
    aes(
      label = paste0(
        round(Contribution_percent, 1),
        "%"
      )
    ),
    position = position_dodge(width = 0.8),
    vjust = -0.35,
    size = 3
  ) +
  
  scale_fill_manual(
    values = c(
      "S1" = "#67BFA5",
      "S2" = "#FC8D62",
      "S3" = "#8DA0CB"
    ),
    
    labels = c(
      "S1: Known pre-1973 ancestors",
      "S2: Ancestors including unknown years",
      "S3: All pre-1973 horses"
    )
  ) +
  
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
      paste0(
        "Comparison of three definitions of the pre-1973 genetic pool ",
        "(top ",
        TOP_N,
        " ancestors)"
      ),
    
    x = "Reference population",
    
    y = "Marginal genetic contribution (%)",
    
    fill = "Historical pool definition"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 12
    ),
    
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    
    legend.position = "right"
  )

print(contribution_plot)

# 27. SAVE THE GRAPH ----------------------------------------------

ggsave(
  filename = file.path(
    output_folder,
    "historical_marginal_contribution_TOP20_2015_2024.png"
  ),
  plot = contribution_plot,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_folder,
    "historical_marginal_contribution_TOP20_2015_2024.pdf"
  ),
  plot = contribution_plot,
  width = 11,
  height = 7
)

# 28. SAVE R AND PACKAGE VERSION INFORMATION ----------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_folder,
    "sessionInfo.txt"
  )
)

# 29. FINAL MESSAGE -----------------------------------------------

cat(
  "\nAnalysis complete.\n",
  "The CFC workflow was skipped.\n",
  "The final reference population includes only animals born ",
  "from 2015 through 2024.\n",
  "No animals born in 2025 or later were included.\n",
  "Results were saved in:\n",
  output_folder,
  "\n",
  sep = ""
)

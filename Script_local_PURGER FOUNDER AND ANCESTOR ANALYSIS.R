# ================================================================
# PURGER FOUNDER AND ANCESTOR ANALYSIS
# SOUTH AFRICAN BOERPERD
#
# Input pedigree:
# SAH_ped_PurgeR_input_loops_res_2025_excl.csv
#
# Reference population:
# All animals with STATUS beginning with "REG"
#
# Outputs:
# 1. PurgeR founder/ancestor statistics
# 2. Exact list of contributing founders
# 3. Birth dates/years of contributing founders
# 4. Exact list of ancestors
# 5. Birth dates/years of ancestors
#
# IMPORTANT:
# Nfe and Nae are effective numbers, not literal subsets of horses.
# Therefore birth dates can be extracted for Nf and Na,
# but not directly for Nfe and Nae.
# ================================================================


# ================================================================
# 1. PACKAGES
# ================================================================

# Install once if required:
# install.packages(c("purgeR", "dplyr"))

library(purgeR)
library(dplyr)


# ================================================================
# 2. WORKING DIRECTORY
# ================================================================

setwd(
  "C:/2025/Paper drafts/Pedigree paper/1973 cohort analysis check"
)


# ================================================================
# 3. READ LOOP-RESOLVED PURGER PEDIGREE
# ================================================================

ped_raw <- read.csv(
  "SAH_ped_PurgeR_input_loops_res_2025_excl.csv",
  colClasses = "character",
  stringsAsFactors = FALSE
)


# ================================================================
# 4. CHECK REQUIRED COLUMNS
# ================================================================

required_columns <- c(
  "ANI_ID",
  "ANI_ID_SIRE",
  "ANI_ID_DAM",
  "BIRTH_DTM",
  "SEX",
  "STATUS"
)


if (!all(required_columns %in% names(ped_raw))) {
  
  stop(
    "Missing required columns: ",
    paste(
      setdiff(
        required_columns,
        names(ped_raw)
      ),
      collapse = ", "
    )
  )
}


# Check that individual IDs are unique
if (anyDuplicated(ped_raw$ANI_ID) > 0) {
  
  stop(
    "Duplicate ANI_ID values are present in the pedigree."
  )
}


# ================================================================
# 5. CHECK STATUS CATEGORIES
# ================================================================

cat(
  "\nSTATUS categories in the pedigree:\n"
)

print(
  table(
    ped_raw$STATUS,
    useNA = "ifany"
  )
)


# ================================================================
# 6. DEFINE THE REFERENCE POPULATION
# ================================================================

# Include all animals whose STATUS starts with REG.
#
# This includes coding such as:
# REG
# REG (Registered)

ped_raw$reference <- (
  !is.na(ped_raw$STATUS) &
    grepl(
      "^REG\\b",
      trimws(ped_raw$STATUS),
      ignore.case = TRUE
    )
)


cat(
  "\nNumber of REG animals in reference population:",
  sum(ped_raw$reference),
  "\n"
)


# ================================================================
# 7. PREPARE PEDIGREE FOR PURGER
# ================================================================

# ped_sort() orders ancestors before descendants and renames
# animals internally from 1 to N.
#
# keep_names = TRUE keeps the original ANI_ID in a column
# called "names".

ped <- purgeR::ped_sort(
  ped_raw,
  id = "ANI_ID",
  dam = "ANI_ID_DAM",
  sire = "ANI_ID_SIRE",
  keep_names = TRUE
)


cat(
  "\nNumber of animals in PurgeR pedigree:",
  nrow(ped),
  "\n"
)


cat(
  "Number of REG animals after ped_sort:",
  sum(
    ped$reference,
    na.rm = TRUE
  ),
  "\n"
)


# ================================================================
# 8. RUN PURGER FOUNDER/ANCESTOR ANALYSIS
# ================================================================

set.seed(
  1973
)


ancestor_results <- purgeR::pop_Nancestors(
  ped = ped,
  reference = "reference",
  nboot = 10000,
  seed = 1973,
  skip_Ng = FALSE
)


# ================================================================
# 9. PRINT PURGER RESULTS
# ================================================================

cat(
  "\n============================================\n",
  "PURGER FOUNDER AND ANCESTOR RESULTS\n",
  "============================================\n"
)


print(
  ancestor_results
)


# Save results
write.csv(
  ancestor_results,
  "PurgeR_REG_founder_ancestor_results.csv",
  row.names = FALSE
)


# ================================================================
# 10. IDENTIFY REFERENCE-POPULATION IDS
# ================================================================

# These are the PurgeR-renamed IDs for all REG animals.

rp_ids <- ped$id[
  !is.na(ped$reference) &
    ped$reference == TRUE
]


cat(
  "\nReference population size from pedigree:",
  length(rp_ids),
  "\n"
)


cat(
  "PurgeR Nr:",
  ancestor_results$Nr,
  "\n"
)


# Cross-check Nr
if (length(rp_ids) != ancestor_results$Nr) {
  
  stop(
    "Reference population count does not match PurgeR Nr."
  )
}


# ================================================================
# 11. TRACE ALL PARENTS, GRANDPARENTS, ETC.
# ================================================================

# Start at the REG reference animals.
#
# At each iteration:
# reference animals
# -> parents
# -> grandparents
# -> great-grandparents
# -> etc.
#
# PurgeR uses 0 for unknown parents after ped_sort().

ancestor_ped_ids <- integer(0)

current_ids <- rp_ids


repeat {
  
  # Find rows corresponding to current animals
  current_rows <- match(
    current_ids,
    ped$id
  )
  
  
  current_rows <- current_rows[
    !is.na(current_rows)
  ]
  
  
  # Retrieve dams and sires
  parents <- unique(
    c(
      ped$dam[current_rows],
      ped$sire[current_rows]
    )
  )
  
  
  # Remove unknown parents
  parents <- parents[
    !is.na(parents) &
      parents != 0
  ]
  
  
  # Retain only newly discovered ancestors
  new_ancestors <- setdiff(
    parents,
    ancestor_ped_ids
  )
  
  
  # Stop when no additional ancestors are found
  if (length(new_ancestors) == 0) {
    break
  }
  
  
  # Add newly discovered ancestors
  ancestor_ped_ids <- unique(
    c(
      ancestor_ped_ids,
      new_ancestors
    )
  )
  
  
  # Move another generation backwards
  current_ids <- new_ancestors
}


cat(
  "\nAncestors reached through parent tracing:",
  length(ancestor_ped_ids),
  "\n"
)


# ================================================================
# 12. IDENTIFY FOUNDERS AMONG TRACED ANCESTORS
# ================================================================

# PurgeR defines founders as animals with:
#
# dam = 0
# sire = 0


ancestor_rows <- match(
  ancestor_ped_ids,
  ped$id
)


ancestor_founders <- ancestor_ped_ids[
  ped$dam[ancestor_rows] == 0 &
    ped$sire[ancestor_rows] == 0
]


cat(
  "Founders reached through parent tracing:",
  length(ancestor_founders),
  "\n"
)


# ================================================================
# 13. IDENTIFY REG ANIMALS THAT ARE THEMSELVES FOUNDERS
# ================================================================

# THIS IS THE IMPORTANT DIFFERENCE WE DISCOVERED.
#
# Some animals belong to the REG reference population AND
# have both parents unknown.
#
# These reference-population founders must also be included
# in the PurgeR-consistent founder/ancestor reconstruction.

rp_rows <- match(
  rp_ids,
  ped$id
)


rp_founders <- rp_ids[
  ped$dam[rp_rows] == 0 &
    ped$sire[rp_rows] == 0
]


cat(
  "REG animals that are themselves founders:",
  length(rp_founders),
  "\n"
)


# ================================================================
# 14. FIND REG FOUNDERS NOT REACHED THROUGH PARENT TRACING
# ================================================================

rp_founders_not_traced <- setdiff(
  rp_founders,
  ancestor_ped_ids
)


cat(
  "REG founders not reached through parent tracing:",
  length(rp_founders_not_traced),
  "\n"
)


# ================================================================
# 15. CREATE THE FINAL PURGER-CONSISTENT FOUNDER LIST
# ================================================================

# Combine:
#
# 1. founders reached through ancestry tracing
# 2. REG animals that are themselves founders


founder_ped_ids <- sort(
  unique(
    c(
      ancestor_founders,
      rp_founders
    )
  )
)


cat(
  "\nPurgeR Nf:",
  ancestor_results$Nf,
  "\n"
)


cat(
  "Extracted Nf:",
  length(founder_ped_ids),
  "\n"
)


# ================================================================
# 16. CREATE FINAL PURGER-CONSISTENT ANCESTOR LIST
# ================================================================

# Parent tracing does not necessarily include a REG animal
# that is itself a founder.
#
# Therefore add the RP founders to the ancestor set.

ancestor_ped_ids_purger <- sort(
  unique(
    c(
      ancestor_ped_ids,
      rp_founders
    )
  )
)


cat(
  "\nPurgeR Na:",
  ancestor_results$Na,
  "\n"
)


cat(
  "Extracted Na:",
  length(ancestor_ped_ids_purger),
  "\n"
)


# ================================================================
# 17. CRITICAL CROSS-CHECK
# ================================================================

# Do NOT continue to birth-date extraction unless the manually
# reconstructed lists reproduce PurgeR exactly.


if (
  length(founder_ped_ids) !=
  ancestor_results$Nf
) {
  
  stop(
    paste0(
      "STOP: Extracted founder count (",
      length(founder_ped_ids),
      ") does not match PurgeR Nf (",
      ancestor_results$Nf,
      ")."
    )
  )
}


if (
  length(ancestor_ped_ids_purger) !=
  ancestor_results$Na
) {
  
  stop(
    paste0(
      "STOP: Extracted ancestor count (",
      length(ancestor_ped_ids_purger),
      ") does not match PurgeR Na (",
      ancestor_results$Na,
      ")."
    )
  )
}


cat(
  "\n============================================\n",
  "CROSS-CHECK PASSED\n",
  "Founders and ancestors match PurgeR exactly.\n",
  "============================================\n"
)


# ================================================================
# 18. CREATE PURGER ID -> ORIGINAL ANI_ID LOOKUP
# ================================================================

# ped$id = PurgeR's renamed numerical ID
# ped$names = original ANI_ID because keep_names = TRUE


id_lookup <- data.frame(
  
  PurgeR_ID = ped$id,
  
  ANI_ID = as.character(
    ped[["names"]]
  ),
  
  stringsAsFactors = FALSE
)


# ================================================================
# 19. CREATE ORIGINAL METADATA LOOKUP
# ================================================================

animal_metadata <- ped_raw %>%
  transmute(
    
    ANI_ID = as.character(
      ANI_ID
    ),
    
    BIRTH_DTM = as.character(
      BIRTH_DTM
    ),
    
    SEX = as.character(
      SEX
    ),
    
    STATUS = as.character(
      STATUS
    )
  )


# ================================================================
# 20. EXTRACT EXACT FOUNDER LIST AND BIRTH DATES
# ================================================================

founder_list <- id_lookup %>%
  
  filter(
    PurgeR_ID %in% founder_ped_ids
  ) %>%
  
  left_join(
    animal_metadata,
    by = "ANI_ID"
  ) %>%
  
  mutate(
    
    BirthYear = suppressWarnings(
      as.integer(
        substr(
          BIRTH_DTM,
          1,
          4
        )
      )
    )
  ) %>%
  
  arrange(
    is.na(BirthYear),
    BirthYear,
    ANI_ID
  )


# Final founder count check
cat(
  "\nFounder rows in final date file:",
  nrow(founder_list),
  "\n"
)


if (
  nrow(founder_list) !=
  ancestor_results$Nf
) {
  
  stop(
    "Founder date table does not contain exactly Nf animals."
  )
}


# ================================================================
# 21. SAVE FOUNDER BIRTH-DATE FILE
# ================================================================

write.csv(
  founder_list,
  "PurgeR_REG_founders_birthdates.csv",
  row.names = FALSE
)


# ================================================================
# 22. EXTRACT EXACT ANCESTOR LIST AND BIRTH DATES
# ================================================================

ancestor_list <- id_lookup %>%
  
  filter(
    PurgeR_ID %in% ancestor_ped_ids_purger
  ) %>%
  
  left_join(
    animal_metadata,
    by = "ANI_ID"
  ) %>%
  
  mutate(
    
    BirthYear = suppressWarnings(
      as.integer(
        substr(
          BIRTH_DTM,
          1,
          4
        )
      )
    )
  ) %>%
  
  arrange(
    is.na(BirthYear),
    BirthYear,
    ANI_ID
  )


# Final ancestor count check
cat(
  "\nAncestor rows in final date file:",
  nrow(ancestor_list),
  "\n"
)


if (
  nrow(ancestor_list) !=
  ancestor_results$Na
) {
  
  stop(
    "Ancestor date table does not contain exactly Na animals."
  )
}


# ================================================================
# 23. SAVE ANCESTOR BIRTH-DATE FILE
# ================================================================

write.csv(
  ancestor_list,
  "PurgeR_REG_ancestors_birthdates.csv",
  row.names = FALSE
)


# ================================================================
# 24. SUMMARISE FOUNDER BIRTH YEARS
# ================================================================

cat(
  "\n============================================\n",
  "FOUNDER BIRTH-YEAR SUMMARY\n",
  "============================================\n"
)


print(
  summary(
    founder_list$BirthYear
  )
)


cat(
  "\nTotal founders:",
  nrow(founder_list),
  "\n"
)


cat(
  "Founders with known birth year:",
  sum(
    !is.na(founder_list$BirthYear)
  ),
  "\n"
)


cat(
  "Founders with unknown birth year:",
  sum(
    is.na(founder_list$BirthYear)
  ),
  "\n"
)


if (
  any(
    !is.na(founder_list$BirthYear)
  )
) {
  
  cat(
    "Founder birth-year range:",
    min(
      founder_list$BirthYear,
      na.rm = TRUE
    ),
    "-",
    max(
      founder_list$BirthYear,
      na.rm = TRUE
    ),
    "\n"
  )
}


# ================================================================
# 25. SUMMARISE ANCESTOR BIRTH YEARS
# ================================================================

cat(
  "\n============================================\n",
  "ANCESTOR BIRTH-YEAR SUMMARY\n",
  "============================================\n"
)


print(
  summary(
    ancestor_list$BirthYear
  )
)


cat(
  "\nTotal ancestors:",
  nrow(ancestor_list),
  "\n"
)


cat(
  "Ancestors with known birth year:",
  sum(
    !is.na(ancestor_list$BirthYear)
  ),
  "\n"
)


cat(
  "Ancestors with unknown birth year:",
  sum(
    is.na(ancestor_list$BirthYear)
  ),
  "\n"
)


if (
  any(
    !is.na(ancestor_list$BirthYear)
  )
) {
  
  cat(
    "Ancestor birth-year range:",
    min(
      ancestor_list$BirthYear,
      na.rm = TRUE
    ),
    "-",
    max(
      ancestor_list$BirthYear,
      na.rm = TRUE
    ),
    "\n"
  )
}


# ================================================================
# 26. SAVE A SIMPLE DATE SUMMARY TABLE
# ================================================================

birthdate_summary <- data.frame(
  
  Group = c(
    "Founders",
    "Ancestors"
  ),
  
  Total = c(
    nrow(founder_list),
    nrow(ancestor_list)
  ),
  
  Known_birth_year = c(
    sum(!is.na(founder_list$BirthYear)),
    sum(!is.na(ancestor_list$BirthYear))
  ),
  
  Unknown_birth_year = c(
    sum(is.na(founder_list$BirthYear)),
    sum(is.na(ancestor_list$BirthYear))
  )
)


write.csv(
  birthdate_summary,
  "PurgeR_REG_founder_ancestor_birthdate_summary.csv",
  row.names = FALSE
)


# ================================================================
# 27. PRINT EFFECTIVE NUMBERS SEPARATELY
# ================================================================

cat(
  "\n============================================\n",
  "EFFECTIVE FOUNDER/ANCESTOR NUMBERS\n",
  "============================================\n"
)


cat(
  "Actual founders (Nf):",
  ancestor_results$Nf,
  "\n"
)


cat(
  "Effective number of founders (Nfe):",
  ancestor_results$Nfe,
  "\n"
)


cat(
  "Actual ancestors (Na):",
  ancestor_results$Na,
  "\n"
)


cat(
  "Effective number of ancestors (Nae):",
  ancestor_results$Nae,
  "\n"
)


cat(
  "Founder genome equivalents (Ng):",
  ancestor_results$Ng,
  "\n"
)


# ================================================================
# 28. FINISHED
# ================================================================

cat(
  "\n============================================\n",
  "ANALYSIS COMPLETE\n",
  "============================================\n"
)

cat(
  "\nFiles created:\n",
  "PurgeR_REG_founder_ancestor_results.csv\n",
  "PurgeR_REG_founders_birthdates.csv\n",
  "PurgeR_REG_ancestors_birthdates.csv\n",
  "PurgeR_REG_founder_ancestor_birthdate_summary.csv\n"
)

sum(is.na(founder_list$BirthYear))
sum(is.na(ancestor_list$BirthYear))

range(founder_list$BirthYear, na.rm = TRUE)
range(ancestor_list$BirthYear, na.rm = TRUE)

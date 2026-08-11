# Pedigree analysis: 1973 cohort cross-check
# Create three CFC files using numeric group codes only.

setwd("C:/2025/Paper drafts/Pedigree paper/1973 cohort analysis check")

# Read the pedigree in which loops have already been resolved.
pedigree <- read.csv(
  "SAH_basic_ped_input_ungrouped_loopsres.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Confirm that the columns needed for ancestry tracing are present.
required <- c("ANI_ID", "SIRE_ID", "DAM_ID", "BIRTH_DTM")

if (!all(required %in% names(pedigree))) {
  stop(
    "Missing columns: ",
    paste(setdiff(required, names(pedigree)), collapse = ", ")
  )
}

# Convert IDs to character so animal and parent IDs match correctly.
pedigree$ANI_ID  <- trimws(as.character(pedigree$ANI_ID))
pedigree$SIRE_ID <- trimws(as.character(pedigree$SIRE_ID))
pedigree$DAM_ID  <- trimws(as.character(pedigree$DAM_ID))

# Convert common unknown-parent codes to missing values.
unknown <- c("", "0", "-9", "*", "NA", "N/A", "UNKNOWN")

pedigree$SIRE_ID[pedigree$SIRE_ID %in% unknown] <- NA
pedigree$DAM_ID[pedigree$DAM_ID %in% unknown]   <- NA

# Each animal should occur only once.
if (anyDuplicated(pedigree$ANI_ID)) {
  stop("Duplicated ANI_ID values were found.")
}

# Convert the birth date and extract the birth year.
birth_date <- as.Date(
  pedigree$BIRTH_DTM,
  format = "%Y/%m/%d"
)

birth_year <- as.integer(
  format(birth_date, "%Y")
)

# Define the 1973 cohort as all animals born in 1973.
cohort_1973 <- pedigree$ANI_ID[
  !is.na(birth_year) &
    birth_year == 1973
]

# Trace all recorded parents, grandparents and earlier ancestors.
get_ancestors <- function(start_ids, ped) {
  
  # Begin with the animals whose ancestry must be traced.
  current <- unique(start_ids)
  
  # Create an empty vector in which all discovered ancestors will be stored.
  ancestors <- character()
  
  repeat {
    
    # Find the pedigree rows for the current animals.
    rows <- ped[
      ped$ANI_ID %in% current,
      c("SIRE_ID", "DAM_ID"),
      drop = FALSE
    ]
    
    # Combine their sire and dam IDs into one parent list.
    parents <- unique(c(
      rows$SIRE_ID,
      rows$DAM_ID
    ))
    
    # Remove missing parents.
    parents <- parents[!is.na(parents)]
    
    # Keep only parents that have not already been traced.
    new_parents <- setdiff(
      parents,
      ancestors
    )
    
    # Stop once no additional ancestors can be found.
    if (length(new_parents) == 0) {
      break
    }
    
    # Add newly discovered parents to the complete ancestor list.
    ancestors <- union(
      ancestors,
      new_parents
    )
    
    # Trace the parents of these animals in the next iteration.
    current <- new_parents
  }
  
  # CFC groups can only be assigned to animals with their own pedigree row.
  return(
    intersect(ancestors, ped$ANI_ID)
  )
}

# Find all recorded ancestors of the 1973 cohort.
ancestors_1973 <- get_ancestors(
  cohort_1973,
  pedigree
)

# Scenario 1:
# confirmed ancestors of the 1973 cohort born before 1973.
S1 <- pedigree$ANI_ID[
  pedigree$ANI_ID %in% ancestors_1973 &
    !is.na(birth_year) &
    birth_year < 1973
]

# Scenario 2:
# S1 plus ancestors of the 1973 cohort with unknown birth years.
S2 <- pedigree$ANI_ID[
  pedigree$ANI_ID %in% ancestors_1973 &
    (
      is.na(birth_year) |
        birth_year < 1973
    )
]

# Scenario 3:
# all animals in the full pedigree born before 1973.
S3 <- pedigree$ANI_ID[
  !is.na(birth_year) &
    birth_year < 1973
]

# Store the three historical pools together.
pools <- list(
  S1 = unique(S1),
  S2 = unique(S2),
  S3 = unique(S3)
)

# Compare the observed pool sizes with the expected values.
pool_check <- data.frame(
  Scenario = names(pools),
  Expected = c(29, 39, 343),
  Observed = lengths(pools)
)

print(pool_check)

# Create one CFC input file for each scenario.
for (scenario in names(pools)) {
  
  # Work on a copy of the pedigree.
  out <- pedigree
  
  # Numeric CFC groups:
  # 1 = 1973 cohort
  # 2 = historical pool
  # 3 = 1974–1999
  # 4 = 2000–2024
  # 5 = all remaining pedigree animals
  out$group <- 5L
  
  # Assign the 1973 cohort.
  out$group[
    !is.na(birth_year) &
      birth_year == 1973
  ] <- 1L
  
  # Assign the 1974–1999 cohort.
  out$group[
    !is.na(birth_year) &
      birth_year >= 1974 &
      birth_year <= 1999
  ] <- 3L
  
  # Assign the 2000–2024 reference population.
  out$group[
    !is.na(birth_year) &
      birth_year >= 2000 &
      birth_year <= 2024
  ] <- 4L
  
  # Assign the scenario-specific historical pool.
  # This is assigned last so unknown-year ancestors in S2 receive Group 2.
  out$group[
    out$ANI_ID %in% pools[[scenario]]
  ] <- 2L
  
  # Export a separate file for each scenario.
  write.csv(
    out,
    paste0(
      "SAH_CFC_",
      scenario,
      "_grouped.csv"
    ),
    row.names = FALSE,
    na = ""
  )
  
  # Display the number of animals in each numeric group.
  cat(
    "\n",
    scenario,
    " group counts:\n",
    sep = ""
  )
  
  print(
    table(out$group)
  )
}

# Save the historical-pool size check.
write.csv(
  pool_check,
  "CFC_historical_pool_size_check.csv",
  row.names = FALSE
)

#output specifically for cfc:
# Set the folder containing the three grouped CSV files.

# Names of the three files created previously.
files <- c(
  "SAH_CFC_S1_grouped.csv",
  "SAH_CFC_S2_grouped.csv",
  "SAH_CFC_S3_grouped.csv"
)

# Process each file.
for (file in files) {
  
  # Read the CSV file.
  pedigree <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Rename the columns to the format required for CFC.
  names(pedigree)[names(pedigree) == "ANI_ID"]  <- "Progeny"
  names(pedigree)[names(pedigree) == "SIRE_ID"] <- "Sire"
  names(pedigree)[names(pedigree) == "DAM_ID"]  <- "Dam"
  
  # Keep only the four columns required for CFC.
  # This removes BIRTH_DTM and all other unnecessary columns.
  pedigree <- pedigree[, c(
    "Progeny",
    "Sire",
    "Dam",
    "group"
  )]
  
  # Create a new output filename.
  output_file <- sub(
    "_grouped.csv",
    "_CFC_format.csv",
    file
  )
  
  # Save the reformatted file without overwriting the grouped file.
  write.csv(
    pedigree,
    output_file,
    row.names = FALSE,
    na = ""
  )
}

#Now we go to CFC for analysis
#Scenario 1
#create a special cut Scenario 1 pedigree in which
#For every Group 2 animal:
#Sire = 0
#Dam = 0

# Read the Scenario 1 CFC input file.
S1 <- read.table(
  "CFC_S1_crosscheck_input.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Confirm that the group column is named Group.
names(S1)[tolower(names(S1)) == "group"] <- "Group"

# Check that Scenario 1 contains exactly 29 historical animals.
sum(S1$Group == 2)

# Cut the pedigree at the Scenario 1 historical pool.
S1$Sire[S1$Group == 2] <- 0
S1$Dam[S1$Group == 2]  <- 0

# Confirm that all Group 2 animals now have unknown parents.
table(
  Sire_zero = S1$Sire[S1$Group == 2] == 0,
  Dam_zero  = S1$Dam[S1$Group == 2] == 0
)

# Export the cut Scenario 1 file.
write.table(
  S1,
  "CFC_S1_crosscheck_CUT_input.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "0"
)

# Read the full Scenario 2 CFC pedigree.
S2_ped <- read.table(
  "CFC_S2_crosscheck_input.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Standardise the group-column name in case it is written as "group".
names(S2_ped)[tolower(names(S2_ped)) == "group"] <- "Group"

# Check that the required columns were read correctly.
names(S2_ped)

# Confirm that Group 2 contains the expected 39 historical animals.
sum(S2_ped$Group == 2)

# Stop if the Scenario 2 pool is not exactly 39 animals.
stopifnot(sum(S2_ped$Group == 2) == 39)

# Cut the pedigree at the 39 Scenario 2 historical animals.
S2_ped$Sire[S2_ped$Group == 2] <- 0
S2_ped$Dam[S2_ped$Group == 2]  <- 0

# Confirm that all 39 historical animals now have Sire = 0 and Dam = 0.
table(
  Sire_zero = S2_ped$Sire[S2_ped$Group == 2] == 0,
  Dam_zero  = S2_ped$Dam[S2_ped$Group == 2] == 0
)

# Save a new cut file without overwriting the original.
write.table(
  S2_ped,
  "CFC_S2_crosscheck_CUT_input.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "0"
)

# Read the full Scenario 3 CFC pedigree.
S3_ped <- read.table(
  "CFC_S3_crosscheck_input.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Standardise the group-column name.
names(S3_ped)[tolower(names(S3_ped)) == "group"] <- "Group"

# Check that the file was read correctly.
names(S3_ped)

# Confirm that Group 2 contains exactly 343 historical animals.
sum(S3_ped$Group == 2)

# Stop if the historical pool is not the expected size.
stopifnot(sum(S3_ped$Group == 2) == 343)

# Cut the pedigree at the 343 Scenario 3 historical animals.
# This makes them analytical founders in the CFC analysis.
S3_ped$Sire[S3_ped$Group == 2] <- 0
S3_ped$Dam[S3_ped$Group == 2]  <- 0

# Confirm that every Group 2 animal now has unknown parents.
table(
  Sire_zero = S3_ped$Sire[S3_ped$Group == 2] == 0,
  Dam_zero  = S3_ped$Dam[S3_ped$Group == 2] == 0
)

# Save the cut Scenario 3 pedigree without overwriting the original.
write.table(
  S3_ped,
  "CFC_S3_crosscheck_CUT_input.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "0"
)

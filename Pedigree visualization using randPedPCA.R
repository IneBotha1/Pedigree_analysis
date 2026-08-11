#Pedigree visualization using randPedPCA

# ================================================================
# 1. CALCULATE ENOUGH PEDIGREE PRINCIPAL COMPONENTS
# ================================================================

library(randPedPCA)
library(adegenet)
library(ggplot2)

set.seed(1973)

ped_pca <- randPedPCA::rppca(
  X = ped_object,
  
  # We need more than six PCs because clustering should not be
  # based only on the first two plotted dimensions.
  rank = 20,
  
  # rspec is generally faster when eight or more PCs are requested.
  method = "rspec",
  
  center = TRUE
)

# ================================================================
# 2. EXTRACT THE PEDIGREE PC SCORES
# ================================================================

pca_scores <- as.data.frame(
  ped_pca$x
)

pca_scores$ANI_ID <- as.character(
  ped_object@label
)

# ================================================================
# 3. TEST K ACROSS SEVERAL NUMBERS OF RETAINED PCs
# ================================================================

pc_numbers <- c(
  5,
  10,
  15,
  20
)

max_k <- 10

cluster_tests <- list()

k_sensitivity <- data.frame(
  Number_of_PCs = pc_numbers,
  Best_K = NA_integer_
)


for (i in seq_along(pc_numbers)) {
  
  number_of_pcs <- pc_numbers[i]
  
  # Select the required pedigree PCs.
  pc_matrix <- pca_scores[
    ,
    paste0(
      "PC",
      1:number_of_pcs
    ),
    drop = FALSE
  ]
  
  set.seed(1973)
  
  # find.clusters tests successive K-means solutions.
  cluster_tests[[i]] <- adegenet::find.clusters(
    x = pc_matrix,
    
    # Retain all supplied pedigree PCs during the internal
    # dimension-reduction step.
    n.pca = number_of_pcs,
    
    # Test K = 1 through K = 10.
    max.n.clust = max_k,
    
    # Do not ask for an interactive manual choice.
    choose.n.clust = FALSE,
    
    # Use BIC to compare clustering solutions.
    stat = "BIC",
    
    # Select the K associated with the lowest BIC.
    criterion = "min",
    
    # Repeat K-means from multiple starting positions.
    n.start = 100,
    
    # Do not give low-variance PCs equal weight by scaling them.
    scale = FALSE
  )
  
  # Kstat contains BIC values for K = 1, 2, ..., max_k.
  k_sensitivity$Best_K[i] <- which.min(
    cluster_tests[[i]]$Kstat
  )
}


print(
  k_sensitivity
)


# ================================================================
# EXTRACT THE CLUSTERING RESULT FOR 10 PCs
# ================================================================

result_10pc <- cluster_tests[[
  which(pc_numbers == 10)
]]


# Check what is inside the result
names(result_10pc)

# ================================================================
# CREATE BIC TABLE
# ================================================================

bic_table <- data.frame(
  
  K = seq_along(
    result_10pc$Kstat
  ),
  
  BIC = as.numeric(
    result_10pc$Kstat
  )
)


print(bic_table)

# ================================================================
# GRAPH BIC AGAINST K
# ================================================================

library(ggplot2)

bic_plot <- ggplot(
  bic_table,
  aes(
    x = K,
    y = BIC
  )
) +
  
  geom_line(
    colour = "blue",
    linewidth = 1
  ) +
  
  geom_point(
    colour = "blue",
    size = 2.8
  ) +
  
  scale_x_continuous(
    breaks = bic_table$K
  ) +
  
  labs(
    title = "BIC for pedigree-based clustering",
    x = "Number of clusters (K)",
    y = "Bayesian information criterion"
  ) +
  
  theme_classic(
    base_size = 13
  )


print(bic_plot)

# ================================================================
# CALCULATE THE REDUCTION IN BIC WHEN EACH CLUSTER IS ADDED
# ================================================================

bic_change <- data.frame(
  
  K = 2:nrow(bic_table),
  
  BIC_reduction = -diff(
    bic_table$BIC
  )
)


print(
  bic_change
)

# ================================================================
# GRAPH THE REDUCTION IN BIC
# ================================================================

library(ggplot2)


bic_change_plot <- ggplot(
  bic_change,
  aes(
    x = K,
    y = BIC_reduction
  )
) +
  
  geom_line(
    colour = "blue",
    linewidth = 1
  ) +
  
  geom_point(
    colour = "blue",
    size = 3
  ) +
  
  scale_x_continuous(
    breaks = bic_change$K
  ) +
  
  labs(
    title = "Improvement in clustering fit with increasing K",
    
    subtitle =
      "Reduction in BIC produced by adding each additional pedigree cluster",
    
    x = "Number of clusters (K)",
    
    y = "Reduction in BIC"
  ) +
  
  theme_classic(
    base_size = 13
  )


print(
  bic_change_plot
)

# ================================================================
# PEDIGREE PCA AND CLUSTERING
# SOUTH AFRICAN BOERPERD
# randPedPCA + adegenet
# ================================================================


# ================================================================
# 1. PACKAGES
# ================================================================

# Install once if required:
# install.packages(c(
#   "optiSel",
#   "pedigreeTools",
#   "randPedPCA",
#   "adegenet",
#   "ggplot2"
# ))

library(optiSel)
library(pedigreeTools)
library(randPedPCA)
library(adegenet)
library(ggplot2)


# ================================================================
# 2. WORKING DIRECTORY
# ================================================================

setwd(
  "C:/2025/Paper drafts/Pedigree paper/1973 cohort analysis check"
)


# ================================================================
# 3. READ COMPLETE LOOP-RESOLVED PEDIGREE
# ================================================================

ped_raw <- read.csv(
  "SAH_basic_ped_input_ungrouped_loopsres.csv",
  colClasses = "character",
  stringsAsFactors = FALSE
)


# Check important columns are present
required <- c(
  "ANI_ID",
  "SIRE_ID",
  "DAM_ID",
  "BIRTH_DTM"
)

if (!all(required %in% names(ped_raw))) {
  
  stop(
    "Missing columns: ",
    paste(
      setdiff(required, names(ped_raw)),
      collapse = ", "
    )
  )
}


# ================================================================
# 4. PREPARE PEDIGREE
# ================================================================

ped <- data.frame(
  
  Indiv = trimws(
    ped_raw$ANI_ID
  ),
  
  Sire = trimws(
    ped_raw$SIRE_ID
  ),
  
  Dam = trimws(
    ped_raw$DAM_ID
  ),
  
  stringsAsFactors = FALSE
)


# Remove rows without an individual ID
ped <- ped[
  !is.na(ped$Indiv) &
    ped$Indiv != "",
]


# Codes representing unknown parents
unknown <- c(
  "",
  "0",
  "-9",
  "*",
  "NA",
  "N/A",
  "UNKNOWN",
  "UNK"
)


ped$Sire[
  toupper(ped$Sire) %in% unknown
] <- NA


ped$Dam[
  toupper(ped$Dam) %in% unknown
] <- NA


# ================================================================
# 5. PREPARE AND ORDER PEDIGREE
# ================================================================

# prePed() ensures the pedigree is correctly structured
# for pedigree analysis and places parents before descendants.

Pedig <- optiSel::prePed(
  ped
)


cat(
  "Number of animals in prepared pedigree:",
  nrow(Pedig),
  "\n"
)


# ================================================================
# 6. CREATE pedigreeTools PEDIGREE OBJECT
# ================================================================

ped_object <- pedigreeTools::pedigree(
  
  sire = Pedig$Sire,
  
  dam = Pedig$Dam,
  
  label = Pedig$Indiv
)


# ================================================================
# 7. RUN CENTRED PEDIGREE PCA
# ================================================================

# IMPORTANT:
# center = TRUE follows the recommendation of Lee et al. (2025).
#
# Without centring, PC1 can strongly reflect time/generation.
# With centring, population/pedigree structure is more likely
# to dominate the leading PCs.
#
# We calculate 20 PCs so that clustering is not restricted
# to only PC1 and PC2.

set.seed(
  1973
)


ped_pca <- randPedPCA::rppca(
  
  X = ped_object,
  
  method = "rspec",
  
  rank = 20,
  
  center = TRUE
)


# ================================================================
# 8. EXAMINE PCA
# ================================================================

summary(
  ped_pca
)


# Confirm that centring was applied
ped_pca$center


# ================================================================
# 9. EXTRACT PCA SCORES
# ================================================================

pca_scores <- as.data.frame(
  ped_pca$x
)


# Make sure PC columns have simple names
colnames(pca_scores) <- paste0(
  "PC",
  seq_len(
    ncol(pca_scores)
  )
)


# Add animal IDs
pca_scores$ANI_ID <- as.character(
  ped_object@label
)


head(
  pca_scores
)


# Save all PCA scores
write.csv(
  pca_scores,
  "randPedPCA_scores_20PC_centered.csv",
  row.names = FALSE
)


# ================================================================
# 10. VARIANCE EXPLAINED
# ================================================================

# Because the input is a pedigree object, randPedPCA should
# calculate variance proportions for the centred analysis.

if (!is.null(ped_pca$varProps)) {
  
  variance_table <- data.frame(
    
    PC = paste0(
      "PC",
      seq_along(
        ped_pca$varProps
      )
    ),
    
    Variance_percent =
      100 *
      as.numeric(
        ped_pca$varProps
      )
  )
  
  
  print(
    variance_table
  )
  
  
  write.csv(
    variance_table,
    "randPedPCA_variance_centered.csv",
    row.names = FALSE
  )
  
} else {
  
  message(
    "Variance proportions were not stored in $varProps. ",
    "Use summary(ped_pca) to inspect the PCA."
  )
}


# ================================================================
# 11. PLOT PC1 VS PC2 BEFORE CLUSTERING
# ================================================================

# This is important.
# First inspect the pedigree structure WITHOUT imposing clusters.

if (!is.null(ped_pca$varProps)) {
  
  pc1_label <- paste0(
    "Pedigree PC1 (",
    round(
      100 * ped_pca$varProps[1],
      2
    ),
    "%)"
  )
  
  
  pc2_label <- paste0(
    "Pedigree PC2 (",
    round(
      100 * ped_pca$varProps[2],
      2
    ),
    "%)"
  )
  
} else {
  
  pc1_label <- "Pedigree PC1"
  pc2_label <- "Pedigree PC2"
}


pca_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2
  )
) +
  
  geom_point(
    alpha = 0.30,
    size = 0.8
  ) +
  
  labs(
    
    title =
      "Pedigree PCA of the South African Boerperd",
    
    subtitle =
      "Centred pedigree PCA using the complete loop-resolved pedigree",
    
    x = pc1_label,
    
    y = pc2_label
  ) +
  
  theme_classic(
    base_size = 13
  )


print(
  pca_plot
)


ggsave(
  "randPedPCA_PC1_PC2_unclustered.png",
  pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)


# ================================================================
# 12. SELECT PCs FOR CLUSTERING
# ================================================================

# We calculated 20 PCs.
#
# For the primary clustering analysis we will initially use
# the first 10 pedigree PCs.
#
# We do NOT cluster using only PC1 and PC2.

n_pc_cluster <- 10


cluster_matrix <- pca_scores[
  ,
  paste0(
    "PC",
    1:n_pc_cluster
  )
]


# ================================================================
# 13. TEST K = 1 TO 20
# ================================================================

# find.clusters performs successive K-means clustering.
#
# IMPORTANT:
# We are NOT going to automatically accept the K with
# the smallest BIC.
#
# We use diffNgroup only as an additional indication.
#
# The actual decision will be based on:
#   1. shape of the BIC curve
#   2. change in BIC
#   3. biological interpretation
#   4. stability across numbers of retained PCs

set.seed(
  1973
)


cluster_test <- adegenet::find.clusters(
  
  x = cluster_matrix,
  
  n.pca = n_pc_cluster,
  
  max.n.clust = 20,
  
  choose.n.clust = FALSE,
  
  method = "kmeans",
  
  stat = "BIC",
  
  criterion = "diffNgroup",
  
  n.start = 100,
  
  center = FALSE,
  
  scale = FALSE
)


# ================================================================
# 14. AUTOMATIC diffNgroup CANDIDATE
# ================================================================

# This is a candidate only.
# It is NOT automatically accepted as the final K.

candidate_k <- nlevels(
  cluster_test$grp
)


cat(
  "\nCandidate K suggested by diffNgroup:",
  candidate_k,
  "\n\n"
)


# ================================================================
# 15. CREATE BIC TABLE
# ================================================================

bic_table <- data.frame(
  
  K = seq_along(
    cluster_test$Kstat
  ),
  
  BIC = as.numeric(
    cluster_test$Kstat
  )
)


# Calculate the improvement in BIC when each additional
# cluster is introduced.

bic_table$BIC_reduction <- c(
  
  NA,
  
  -diff(
    bic_table$BIC
  )
)


print(
  bic_table
)


write.csv(
  bic_table,
  "randPedPCA_BIC_K1_to_K20.csv",
  row.names = FALSE
)


# ================================================================
# 16. GRAPH BIC
# ================================================================

bic_plot <- ggplot(
  bic_table,
  aes(
    x = K,
    y = BIC
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 2.5
  ) +
  
  geom_vline(
    xintercept = candidate_k,
    linetype = "dashed"
  ) +
  
  scale_x_continuous(
    breaks = bic_table$K
  ) +
  
  labs(
    
    title =
      "BIC for pedigree-based clustering",
    
    subtitle = paste0(
      "Dashed line = diffNgroup candidate K = ",
      candidate_k
    ),
    
    x =
      "Number of clusters (K)",
    
    y =
      "Bayesian information criterion"
  ) +
  
  theme_classic(
    base_size = 13
  )


print(
  bic_plot
)


ggsave(
  "randPedPCA_BIC_curve.png",
  bic_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ================================================================
# 17. GRAPH CHANGE IN BIC
# ================================================================

# This graph is especially important for our dataset.
#
# It shows how much improvement is gained every time another
# cluster is added.
#
# We are looking for the point after which the improvement
# becomes substantially smaller.

bic_change_plot <- ggplot(
  
  bic_table[
    !is.na(
      bic_table$BIC_reduction
    ),
  ],
  
  aes(
    x = K,
    y = BIC_reduction
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 2.8
  ) +
  
  scale_x_continuous(
    breaks = 2:20
  ) +
  
  labs(
    
    title =
      "Change in BIC with increasing number of clusters",
    
    subtitle =
      "Reduction in BIC produced by each additional pedigree cluster",
    
    x =
      "Number of clusters (K)",
    
    y =
      "Reduction in BIC"
  ) +
  
  theme_classic(
    base_size = 13
  )


print(
  bic_change_plot
)


ggsave(
  "randPedPCA_BIC_change.png",
  bic_change_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ================================================================
# 18. SENSITIVITY TO NUMBER OF RETAINED PCs
# ================================================================

# We should make sure our inferred K is not simply a consequence
# of choosing exactly 10 PCs.
#
# Test clustering using:
# 5 PCs
# 10 PCs
# 15 PCs
# 20 PCs

pc_numbers <- c(
  5,
  10,
  15,
  20
)


k_sensitivity <- data.frame(
  
  Number_of_PCs = pc_numbers,
  
  diffNgroup_candidate_K = NA_integer_
)


for (i in seq_along(pc_numbers)) {
  
  npc <- pc_numbers[i]
  
  
  test_matrix <- pca_scores[
    ,
    paste0(
      "PC",
      1:npc
    )
  ]
  
  
  set.seed(
    1973
  )
  
  
  test_result <- adegenet::find.clusters(
    
    x = test_matrix,
    
    n.pca = npc,
    
    max.n.clust = 20,
    
    choose.n.clust = FALSE,
    
    method = "kmeans",
    
    stat = "BIC",
    
    criterion = "diffNgroup",
    
    n.start = 50,
    
    center = FALSE,
    
    scale = FALSE
  )
  
  
  k_sensitivity$diffNgroup_candidate_K[i] <-
    nlevels(
      test_result$grp
    )
}


print(
  k_sensitivity
)


write.csv(
  k_sensitivity,
  "randPedPCA_K_sensitivity.csv",
  row.names = FALSE
)


# ================================================================
# 19. CHOOSE FINAL K
# ================================================================

# DO NOT choose K simply because it has the lowest BIC.
#
# Look at:
#
#   bic_table
#   bic_plot
#   bic_change_plot
#   k_sensitivity
#
# Then enter the K that is scientifically supported.

final_k <- as.integer(
  readline(
    prompt =
      "Enter the final K after examining BIC and BIC change: "
  )
)

#choose 2 and then later 6

cat(
  "\nFinal K selected:",
  final_k,
  "\n"
)


# ================================================================
# 20. RUN FINAL CLUSTERING
# ================================================================

set.seed(
  1973
)


final_clusters <- adegenet::find.clusters(
  
  x = cluster_matrix,
  
  n.pca = n_pc_cluster,
  
  n.clust = final_k,
  
  method = "kmeans",
  
  n.start = 100,
  
  center = FALSE,
  
  scale = FALSE
)


# Add cluster assignments to PCA scores

pca_scores$Cluster <- factor(
  final_clusters$grp
)


# Number of horses in each cluster

cat(
  "\nNumber of individuals in each cluster:\n"
)


print(
  table(
    pca_scores$Cluster
  )
)



# ================================================================
# 21. PLOT FINAL CLUSTERS
# ================================================================

cluster_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    colour = Cluster
  )
) +
  
  geom_point(
    alpha = 0.55,
    size = 0.9
  ) +
  
  labs(
    
    title =
      "Pedigree-based structure of the South African Boerperd",
    
    subtitle = paste0(
      "Centred randPedPCA; K = ",
      final_k
    ),
    
    x = pc1_label,
    
    y = pc2_label,
    
    colour =
      "Pedigree cluster"
  ) +
  
  theme_classic(
    base_size = 13
  )



print(
  cluster_plot
)


ggsave(
  paste0(
    "randPedPCA_clusters_K",
    final_k,
    ".png"
  ),
  cluster_plot,
  width = 9,
  height = 7,
  dpi = 300
)


# ================================================================
# 22. ADD BIRTH YEAR
# ================================================================

# This allows us to test the important Lee et al. (2025)
# question:
#
# Are our pedigree PCs describing pedigree structure,
# or are they still strongly tracking time/generation?


# Birth year from BIRTH_DTM
birth_year <- suppressWarnings(
  as.integer(
    substr(
      ped_raw$BIRTH_DTM,
      1,
      4
    )
  )
)


birth_lookup <- data.frame(
  
  ANI_ID = trimws(
    ped_raw$ANI_ID
  ),
  
  BirthYear = birth_year,
  
  stringsAsFactors = FALSE
)


pca_scores$BirthYear <- birth_lookup$BirthYear[
  
  match(
    pca_scores$ANI_ID,
    birth_lookup$ANI_ID
  )
]


# ================================================================
# 23. TEST CORRELATION BETWEEN PCs AND BIRTH YEAR
# ================================================================

year_correlations <- data.frame(
  PC = paste0("PC", 1:10),
  Spearman_rho = NA_real_
)

for (i in 1:10) {
  
  year_correlations$Spearman_rho[i] <- cor(
    
    pca_scores[[paste0("PC", i)]],
    
    pca_scores$BirthYear,
    
    method = "spearman",
    
    use = "complete.obs"
  )
}

print(year_correlations)


write.csv(
  year_correlations,
  "randPedPCA_PC_birthyear_correlations.csv",
  row.names = FALSE
)


# ================================================================
# 24. VISUALISE BIRTH YEAR ON THE PCA
# ================================================================

year_plot <- ggplot(
  
  pca_scores[
    !is.na(
      pca_scores$BirthYear
    ),
  ],
  
  aes(
    x = PC1,
    y = PC2,
    colour = BirthYear
  )
) +
  
  geom_point(
    alpha = 0.5,
    size = 0.9
  ) +
  
  labs(
    
    title =
      "Pedigree PCA according to year of birth",
    
    subtitle =
      "Used to assess temporal structure in the centred pedigree PCA",
    
    x = pc1_label,
    
    y = pc2_label,
    
    colour =
      "Birth year"
  ) +
  
  theme_classic(
    base_size = 13
  )


print(
  year_plot
)


ggsave(
  "randPedPCA_PC1_PC2_birthyear.png",
  year_plot,
  width = 9,
  height = 7,
  dpi = 300
)


# ================================================================
# 25. BIRTH YEAR DISTRIBUTION WITHIN CLUSTERS
# ================================================================

cluster_year_summary <- aggregate(
  
  BirthYear ~ Cluster,
  
  data = pca_scores,
  
  FUN = function(x) {
    
    c(
      N = length(x),
      Mean = mean(
        x,
        na.rm = TRUE
      ),
      Median = median(
        x,
        na.rm = TRUE
      ),
      Min = min(
        x,
        na.rm = TRUE
      ),
      Max = max(
        x,
        na.rm = TRUE
      )
    )
  }
)


print(
  cluster_year_summary
)


# ================================================================
# 26. SAVE FINAL DATASET
# ================================================================

write.csv(
  pca_scores,
  paste0(
    "randPedPCA_scores_clusters_K",
    final_k,
    ".csv"
  ),
  row.names = FALSE
)


# Save complete PCA object
saveRDS(
  ped_pca,
  "randPedPCA_centered_result.rds"
)


# Save clustering object
saveRDS(
  final_clusters,
  paste0(
    "randPedPCA_clustering_K",
    final_k,
    ".rds"
  )
)


# ================================================================
# DONE
# ================================================================

cat(
  "\nPedigree PCA and clustering analysis complete.\n"
)

# ============================================================================
# STAT 467 - CLUSTERING ANALYSIS
# ============================================================================
# Purpose: Identify natural groupings of countries based on health/socioeconomic
#          indicators using unsupervised learning methods
#
# Input:   data_country_level.csv (country-level aggregated data)
#
# Output:  - Elbow plot (figures/clustering_elbow.png)
#          - Silhouette method (figures/clustering_silhouette.png)
#          - Gap statistic (figures/clustering_gap.png)
#          - Dendrogram (figures/clustering_dendrogram.png)
#          - K-means visualization (figures/clustering_kmeans.png)
#          - Silhouette plot (figures/clustering_silhouette_plot.png)
#          - Cluster profiles (figures/cluster_profiles.csv, figures/cluster_profiles_plot.png)
#          - K-means centers (figures/kmeans_centers.csv)
#          - Cluster membership (figures/cluster_membership.csv)
#
# Methods:
#   1. Optimal k determination: Elbow, Silhouette, Gap statistic
#   2. Hierarchical clustering (Ward's method)
#   3. K-means partitioning
#   4. Cluster validation (silhouette, ARI)
#
# Key Concepts:
#   Clustering Goal:
#     - Maximize within-cluster similarity (minimize WSS)
#     - Maximize between-cluster dissimilarity
#     - Find "natural" groupings without using class labels
#
#   Optimal k Selection:
#     1. Elbow method: Look for "bend" in WSS plot
#        - Subjective: the "elbow" is not always clear
#     2. Silhouette: Measures how similar objects are to their own cluster
#        - Ranges from -1 to 1; higher is better
#        - Interpretation: > 0.7 strong, > 0.5 reasonable, > 0.25 weak
#     3. Gap statistic (Tibshirani et al., 2001): Compares log(WSS) to expectation
#        under null (uniform) distribution
#        - Most rigorous method; provides confidence interval
#        - firstSEmax: first k where gap(k) >= gap(k+1) - SE(k+1)
#
#   Hierarchical vs K-Means:
#     - Hierarchical: Builds tree structure, no need to pre-specify k
#       * Ward's method minimizes within-cluster variance (like K-means)
#     - K-means: Fast, spherical clusters, requires k upfront
#       * Sensitive to initialization → use nstart > 1
#
# Dependencies: tidyverse, cluster, factoextra, dendextend, mclust
# ============================================================================

library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(mclust)

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n\n")

# Select variables for clustering
# These capture health outcomes and socioeconomic development
cluster_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
                  "Income_composition", "BMI", "HIV_AIDS", "Diphtheria")

df_cluster <- df %>% select(Country, Status, all_of(cluster_vars)) %>% drop_na()
cat("Complete cases:", nrow(df_cluster), "\n")

# ============================================================================
# SECTION 1: DATA STANDARDIZATION
# ============================================================================
# Standardization is CRITICAL for clustering:
#   - Variables with larger scales would dominate distance calculations
#   - Z-score: (x - mean) / sd → mean = 0, sd = 1
#   - Ensures all variables contribute equally to distance
#
# Note: We standardize across ALL observations (global standardization)
#       This preserves relative differences between observations

X_scaled <- scale(df_cluster[, cluster_vars])
rownames(X_scaled) <- df_cluster$Country

# ============================================================================
# SECTION 2: OPTIMAL NUMBER OF CLUSTERS (k)
# ============================================================================
# Three methods to determine optimal k:
#   1. Elbow method (WSS)
#   2. Silhouette method
#   3. Gap statistic (most rigorous)

cat("\n=== OPTIMAL K ===\n")

# ----------------------------------------------------------------------------
# 2.1 Within-Sum-of-Squares (WSS) - Elbow Method
# ----------------------------------------------------------------------------
# WSS = Σᵢ Σⱼ ||xⱼ - μᵢ||² (sum of squared distances to cluster centers)
# WSS decreases as k increases (more clusters = less variance within)
# Look for "elbow" where adding clusters provides diminishing returns

wss <- sapply(1:10, function(k) kmeans(X_scaled, k, nstart = 25)$tot.withinss)

# ----------------------------------------------------------------------------
# 2.2 Silhouette Method
# ----------------------------------------------------------------------------
# For each point i:
#   a(i) = average distance to other points in same cluster
#   b(i) = average distance to points in nearest neighboring cluster
#   s(i) = (b(i) - a(i)) / max(a(i), b(i))
#
# Silhouette ranges from -1 to +1:
#   +1: well matched to own cluster, poorly matched to neighbors
#    0: on the border between clusters
#   -1: probably assigned to wrong cluster

sil_width <- sapply(2:10, function(k) {
  km <- kmeans(X_scaled, k, nstart = 25)
  mean(silhouette(km$cluster, dist(X_scaled))[, 3])
})

cat("Max silhouette at k =", which.max(sil_width) + 1, "\n")

# Plot elbow
p_elbow <- ggplot(data.frame(k = 1:10, WSS = wss), aes(x = k, y = WSS)) +
  geom_line(color = "steelblue", linewidth = 1.2) + geom_point(size = 3) +
  scale_x_continuous(breaks = 1:10) + labs(title = "Elbow Method") + theme_minimal()
ggsave("figures/clustering_elbow.png", p_elbow, width = 10, height = 6, dpi = 150)

# Plot silhouette
p_sil <- ggplot(data.frame(k = 2:10, Silhouette = sil_width), aes(x = k, y = Silhouette)) +
  geom_line(color = "steelblue", linewidth = 1.2) + geom_point(size = 3) +
  scale_x_continuous(breaks = 2:10) + labs(title = "Silhouette Method") + theme_minimal()
ggsave("figures/clustering_silhouette.png", p_sil, width = 10, height = 6, dpi = 150)

# ----------------------------------------------------------------------------
# 2.3 Gap Statistic (Tibshirani et al., 2001)
# ----------------------------------------------------------------------------
# Gap(k) = E*[log(WSS_k)] - log(WSS_k)
#   where E* is expectation under null reference (uniform distribution)
#
# Compares actual WSS to expected WSS from random uniform data
# Optimal k: first k where Gap(k) ≥ Gap(k+1) - SE(k+1) (firstSEmax rule)
#
# This is the most rigorous method as it accounts for sampling variability

set.seed(123)
gap <- clusGap(X_scaled, FUN = kmeans, nstart = 25, K.max = 10, B = 50)
ggsave("figures/clustering_gap.png", fviz_gap_stat(gap) + theme_minimal(), width = 10, height = 6, dpi = 150)

# Determine optimal k using Gap statistic (Tibshirani et al., 2001)
# firstSEmax method: first k where gap(k) >= gap(k+1) - SE(k+1)
optimal_k_gap <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")
optimal_k_sil <- which.max(sil_width) + 1

cat("\nOptimal k selection:\n")
cat("  Gap statistic (firstSEmax):", optimal_k_gap, "\n")
cat("  Silhouette method:", optimal_k_sil, "\n")

# Use Gap statistic as primary criterion (more rigorous)
optimal_k <- optimal_k_gap
cat("Using k =", optimal_k, "(based on Gap statistic)\n")

# ============================================================================
# SECTION 3: HIERARCHICAL CLUSTERING
# ============================================================================
# Agglomerative (bottom-up) clustering:
#   1. Start with n clusters (each observation is a cluster)
#   2. Merge closest pair of clusters
#   3. Repeat until one cluster remains
#
# Ward's Method (ward.D2):
#   - Minimizes increase in total within-cluster variance at each merge
#   - Produces compact, spherical clusters (similar to K-means)
#   - Uses squared Euclidean distance
#
# Other linkage methods:
#   - single: nearest neighbor (can produce "chaining")
#   - complete: farthest neighbor (compact but unequal sizes)
#   - average: UPGMA (balanced approach)

cat("\n=== HIERARCHICAL CLUSTERING ===\n")
dist_mat <- dist(X_scaled)
hc <- hclust(dist_mat, method = "ward.D2")

# Dendrogram with cluster rectangles
png("figures/clustering_dendrogram.png", width = 1800, height = 800, res = 100)
plot(hc, main = "Dendrogram (Ward's Method)", xlab = "Countries", cex = 0.4)
rect.hclust(hc, k = optimal_k, border = c("#2E86AB", "#E94F37", "#2EAB5B", "#AB8F2E"))
dev.off()

# Cut tree to get cluster assignments
df_cluster$HC_Cluster <- factor(cutree(hc, k = optimal_k))

# ============================================================================
# SECTION 4: K-MEANS CLUSTERING
# ============================================================================
# Algorithm (Lloyd's):
#   1. Initialize k random centroids
#   2. Assign each point to nearest centroid
#   3. Update centroids as cluster means
#   4. Repeat steps 2-3 until convergence
#
# Important parameters:
#   - nstart = 25: Run algorithm 25 times with different initializations
#                  Helps avoid local minima (K-means is sensitive to initialization)
#   - centers = k: Number of clusters
#
# Assumptions:
#   - Clusters are spherical (equal variance in all directions)
#   - Clusters are similar in size (can fail with very unequal groups)

cat("\n=== K-MEANS ===\n")
set.seed(123)
km <- kmeans(X_scaled, centers = optimal_k, nstart = 25)

cat("Cluster sizes:\n")
print(table(km$cluster))

# BSS/TSS: Between-cluster sum of squares / Total sum of squares
# Higher ratio indicates better separation (more variance explained by clustering)
cat("BSS/TSS:", round(km$betweenss / km$totss * 100, 1), "%\n")

df_cluster$KM_Cluster <- factor(km$cluster)

# Print cluster centers (in standardized units)
# Values > 0 indicate above average; < 0 indicates below average
cat("\nCluster centers:\n")
print(round(km$centers, 3))
write.csv(km$centers, "figures/kmeans_centers.csv")

# Visualize K-means with PCA projection
p_km <- fviz_cluster(km, data = X_scaled,
                     palette = c("#2E86AB", "#E94F37", "#2EAB5B", "#AB8F2E"),
                     ellipse.type = "convex", repel = TRUE, ggtheme = theme_minimal()) +
  labs(title = paste("K-Means (k =", optimal_k, ")"))
ggsave("figures/clustering_kmeans.png", p_km, width = 14, height = 10, dpi = 150)

# ============================================================================
# SECTION 5: CLUSTER VALIDATION
# ============================================================================
# Internal validation: How well-defined are the clusters?
#   - Silhouette width: Measures cohesion vs separation
#   - Cluster size balance: Highly unequal sizes may indicate poor clustering
#
# External validation: How well do clusters match known labels (Status)?
#   - Adjusted Rand Index (ARI): Measures agreement with external classification

cat("\n=== VALIDATION ===\n")

# ----------------------------------------------------------------------------
# 5.1 Silhouette Analysis
# ----------------------------------------------------------------------------
sil_km <- silhouette(km$cluster, dist_mat)
avg_sil <- mean(sil_km[, 3])
cat("Avg silhouette:", round(avg_sil, 3), "\n")

# Check silhouette quality by cluster
# Negative average silhouette for a cluster indicates poor separation
cat("\nSilhouette by cluster:\n")
sil_by_cluster <- aggregate(sil_km[, 3], by = list(Cluster = km$cluster), FUN = mean)
colnames(sil_by_cluster)[2] <- "Avg_Silhouette"
print(round(sil_by_cluster, 3))

# Warn if any cluster has negative average silhouette
neg_sil_clusters <- sil_by_cluster$Cluster[sil_by_cluster$Avg_Silhouette < 0]
if (length(neg_sil_clusters) > 0) {
  cat("WARNING: Cluster(s)", paste(neg_sil_clusters, collapse = ", "),
      "have negative avg silhouette (poor separation)\n")
}

# ----------------------------------------------------------------------------
# 5.2 Cluster Size Balance
# ----------------------------------------------------------------------------
cluster_sizes <- table(km$cluster)
size_ratio <- max(cluster_sizes) / min(cluster_sizes)
cat("\nCluster size ratio (max/min):", round(size_ratio, 2), "\n")
if (size_ratio > 5) {
  cat("WARNING: High size imbalance - consider different k or clustering method\n")
}

# Silhouette plot showing individual observation contributions
ggsave("figures/clustering_silhouette_plot.png",
       fviz_silhouette(sil_km, palette = c("#2E86AB", "#E94F37", "#2EAB5B", "#AB8F2E")) + theme_minimal(),
       width = 12, height = 6, dpi = 150)

# ----------------------------------------------------------------------------
# 5.3 External Validation: Comparison with Status
# ----------------------------------------------------------------------------
# Cross-tabulation shows how clusters align with Developed/Developing status
cat("\nCluster vs Status:\n")
cross_tab <- table(df_cluster$KM_Cluster, df_cluster$Status)
print(cross_tab)

# Adjusted Rand Index (ARI)
# - Measures similarity between cluster assignments and true labels
# - Ranges from -1 to 1: 1 = perfect match, 0 = random, <0 = worse than random
# - "Adjusted" for chance: accounts for expected agreement by random assignment
ari <- adjustedRandIndex(km$cluster, as.numeric(df_cluster$Status))
cat("Adjusted Rand Index:", round(ari, 3), "\n")

# ============================================================================
# SECTION 6: CLUSTER PROFILES
# ============================================================================
# Cluster profiles characterize what makes each cluster unique
# Using GLOBALLY standardized data for proper cross-cluster comparison
#
# CRITICAL: We must use the globally scaled data (X_scaled) for profiles
# Scaling within each cluster separately would make cross-cluster comparison invalid

cat("\n=== CLUSTER PROFILES ===\n")

# Raw (unstandardized) profiles for interpretability
profiles <- df_cluster %>%
  group_by(KM_Cluster) %>%
  summarise(n = n(), across(all_of(cluster_vars), ~round(mean(.), 2)))
print(profiles)
write.csv(profiles, "figures/cluster_profiles.csv", row.names = FALSE)

# CRITICAL FIX: Use globally-scaled data for cluster profiles
# Previous code scaled WITHIN each cluster, making cross-cluster comparison invalid
# Now we use X_scaled (already scaled across ALL data) for proper comparison
X_scaled_df <- as.data.frame(X_scaled)
X_scaled_df$KM_Cluster <- df_cluster$KM_Cluster

# Standardized profiles for comparison
# Positive values = above global average; Negative = below global average
profiles_std <- X_scaled_df %>%
  group_by(KM_Cluster) %>%
  summarise(across(all_of(cluster_vars), mean)) %>%
  pivot_longer(-KM_Cluster, names_to = "Variable", values_to = "Mean")

p_prof <- ggplot(profiles_std, aes(x = Variable, y = Mean, fill = KM_Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#2E86AB", "#E94F37", "#2EAB5B", "#AB8F2E")) +
  labs(title = "Cluster Profiles") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/cluster_profiles_plot.png", p_prof, width = 14, height = 7, dpi = 150)

# ============================================================================
# SECTION 7: CLUSTER MEMBERSHIP
# ============================================================================
# List countries in each cluster for domain interpretation

for (cl in 1:optimal_k) {
  cat("\n--- Cluster", cl, "---\n")
  countries <- df_cluster %>% filter(KM_Cluster == cl) %>% pull(Country)
  status <- df_cluster %>% filter(KM_Cluster == cl) %>% dplyr::count(Status)
  cat("n =", length(countries), ", Status:", paste(paste(status$Status, status$n, sep = "="), collapse = ", "), "\n")
  cat("Examples:", paste(head(countries, 5), collapse = ", "), "\n")
}

write.csv(df_cluster %>% select(Country, Status, KM_Cluster) %>% arrange(KM_Cluster),
          "figures/cluster_membership.csv", row.names = FALSE)

cat("\n=== Clustering Complete ===\n")
cat("Run 08_cca.R next.\n")

# STAT 467 - Clustering
# Input: data_country_level.csv

library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(mclust)

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n\n")

cluster_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
                  "Income_composition", "BMI", "HIV_AIDS", "Diphtheria")

df_cluster <- df %>% select(Country, Status, all_of(cluster_vars)) %>% drop_na()
cat("Complete cases:", nrow(df_cluster), "\n")

X_scaled <- scale(df_cluster[, cluster_vars])
rownames(X_scaled) <- df_cluster$Country

# Optimal k
cat("\n=== OPTIMAL K ===\n")

wss <- sapply(1:10, function(k) kmeans(X_scaled, k, nstart = 25)$tot.withinss)
sil_width <- sapply(2:10, function(k) {
  km <- kmeans(X_scaled, k, nstart = 25)
  mean(silhouette(km$cluster, dist(X_scaled))[, 3])
})

cat("Max silhouette at k =", which.max(sil_width) + 1, "\n")

p_elbow <- ggplot(data.frame(k = 1:10, WSS = wss), aes(x = k, y = WSS)) +
  geom_line(color = "steelblue", linewidth = 1.2) + geom_point(size = 3) +
  scale_x_continuous(breaks = 1:10) + labs(title = "Elbow Method") + theme_minimal()
ggsave("figures/clustering_elbow.png", p_elbow, width = 10, height = 6, dpi = 150)

p_sil <- ggplot(data.frame(k = 2:10, Silhouette = sil_width), aes(x = k, y = Silhouette)) +
  geom_line(color = "steelblue", linewidth = 1.2) + geom_point(size = 3) +
  scale_x_continuous(breaks = 2:10) + labs(title = "Silhouette Method") + theme_minimal()
ggsave("figures/clustering_silhouette.png", p_sil, width = 10, height = 6, dpi = 150)

set.seed(123)
gap <- clusGap(X_scaled, FUN = kmeans, nstart = 25, K.max = 10, B = 50)
ggsave("figures/clustering_gap.png", fviz_gap_stat(gap) + theme_minimal(), width = 10, height = 6, dpi = 150)

optimal_k <- 4
cat("Using k =", optimal_k, "\n")

# Hierarchical clustering
cat("\n=== HIERARCHICAL CLUSTERING ===\n")
dist_mat <- dist(X_scaled)
hc <- hclust(dist_mat, method = "ward.D2")

png("figures/clustering_dendrogram.png", width = 1800, height = 800, res = 100)
plot(hc, main = "Dendrogram (Ward's Method)", xlab = "Countries", cex = 0.4)
rect.hclust(hc, k = optimal_k, border = c("#2E86AB", "#E94F37", "#2EAB5B", "#AB8F2E"))
dev.off()

df_cluster$HC_Cluster <- factor(cutree(hc, k = optimal_k))

# K-means
cat("\n=== K-MEANS ===\n")
set.seed(123)
km <- kmeans(X_scaled, centers = optimal_k, nstart = 25)

cat("Cluster sizes:\n")
print(table(km$cluster))
cat("BSS/TSS:", round(km$betweenss / km$totss * 100, 1), "%\n")

df_cluster$KM_Cluster <- factor(km$cluster)

cat("\nCluster centers:\n")
print(round(km$centers, 3))
write.csv(km$centers, "figures/kmeans_centers.csv")

p_km <- fviz_cluster(km, data = X_scaled,
                     palette = c("#2E86AB", "#E94F37", "#2EAB5B", "#AB8F2E"),
                     ellipse.type = "convex", repel = TRUE, ggtheme = theme_minimal()) +
  labs(title = paste("K-Means (k =", optimal_k, ")"))
ggsave("figures/clustering_kmeans.png", p_km, width = 14, height = 10, dpi = 150)

# Validation
cat("\n=== VALIDATION ===\n")
sil_km <- silhouette(km$cluster, dist_mat)
cat("Avg silhouette:", round(mean(sil_km[, 3]), 3), "\n")

ggsave("figures/clustering_silhouette_plot.png",
       fviz_silhouette(sil_km, palette = c("#2E86AB", "#E94F37", "#2EAB5B", "#AB8F2E")) + theme_minimal(),
       width = 12, height = 6, dpi = 150)

cat("\nCluster vs Status:\n")
cross_tab <- table(df_cluster$KM_Cluster, df_cluster$Status)
print(cross_tab)

ari <- adjustedRandIndex(km$cluster, as.numeric(df_cluster$Status))
cat("Adjusted Rand Index:", round(ari, 3), "\n")

# Cluster profiles
cat("\n=== CLUSTER PROFILES ===\n")
profiles <- df_cluster %>%
  group_by(KM_Cluster) %>%
  summarise(n = n(), across(all_of(cluster_vars), ~round(mean(.), 2)))
print(profiles)
write.csv(profiles, "figures/cluster_profiles.csv", row.names = FALSE)

profiles_std <- df_cluster %>%
  group_by(KM_Cluster) %>%
  summarise(across(all_of(cluster_vars), ~mean(scale(.)[,1]))) %>%
  pivot_longer(-KM_Cluster, names_to = "Variable", values_to = "Mean")

p_prof <- ggplot(profiles_std, aes(x = Variable, y = Mean, fill = KM_Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#2E86AB", "#E94F37", "#2EAB5B", "#AB8F2E")) +
  labs(title = "Cluster Profiles") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/cluster_profiles_plot.png", p_prof, width = 14, height = 7, dpi = 150)

# Cluster membership
for (cl in 1:optimal_k) {
  cat("\n--- Cluster", cl, "---\n")
  countries <- df_cluster %>% filter(KM_Cluster == cl) %>% pull(Country)
  status <- df_cluster %>% filter(KM_Cluster == cl) %>% count(Status)
  cat("n =", length(countries), ", Status:", paste(paste(status$Status, status$n, sep = "="), collapse = ", "), "\n")
  cat("Examples:", paste(head(countries, 5), collapse = ", "), "\n")
}

write.csv(df_cluster %>% select(Country, Status, KM_Cluster) %>% arrange(KM_Cluster),
          "figures/cluster_membership.csv", row.names = FALSE)

cat("\n=== Clustering Complete ===\n")
cat("Run 08_cca.R next.\n")

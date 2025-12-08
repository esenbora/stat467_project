# STAT 467 - Utility Functions

# Mahalanobis distance
calc_mahalanobis <- function(X, center = NULL, cov_mat = NULL) {
  if (is.null(center)) center <- colMeans(X)
  if (is.null(cov_mat)) cov_mat <- cov(X)
  mahalanobis(X, center, cov_mat)
}

# Identify outliers
identify_outliers <- function(X, alpha = 0.001) {
  p <- ncol(X)
  mahal <- calc_mahalanobis(X)
  critical <- qchisq(1 - alpha, df = p)
  list(distances = mahal, critical = critical, outliers = mahal > critical)
}

# One-sample Hotelling's T²
hotelling_one_sample <- function(X, mu0) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  x_bar <- colMeans(X); S <- cov(X)
  diff <- x_bar - mu0
  T2 <- n * t(diff) %*% solve(S) %*% diff
  F_stat <- ((n - p) / (p * (n - 1))) * T2
  p_value <- 1 - pf(F_stat, p, n - p)
  list(T2 = as.numeric(T2), F = as.numeric(F_stat), df1 = p, df2 = n - p, p_value = as.numeric(p_value))
}

# Two-sample Hotelling's T²
hotelling_two_sample <- function(X1, X2) {
  X1 <- as.matrix(X1); X2 <- as.matrix(X2)
  n1 <- nrow(X1); n2 <- nrow(X2); p <- ncol(X1)
  x_bar1 <- colMeans(X1); x_bar2 <- colMeans(X2)
  Sp <- ((n1 - 1) * cov(X1) + (n2 - 1) * cov(X2)) / (n1 + n2 - 2)
  diff <- x_bar1 - x_bar2
  T2 <- (n1 * n2 / (n1 + n2)) * t(diff) %*% solve(Sp) %*% diff
  F_stat <- ((n1 + n2 - p - 1) / (p * (n1 + n2 - 2))) * T2
  p_value <- 1 - pf(F_stat, p, n1 + n2 - p - 1)
  list(T2 = as.numeric(T2), F = as.numeric(F_stat), df1 = p, df2 = n1 + n2 - p - 1, p_value = as.numeric(p_value))
}

# Bonferroni CI
bonferroni_ci <- function(X, alpha = 0.05) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  x_bar <- colMeans(X); se <- apply(X, 2, sd) / sqrt(n)
  t_crit <- qt(1 - alpha / (2 * p), df = n - 1)
  data.frame(Variable = colnames(X), Mean = x_bar, SE = se,
             Lower = x_bar - t_crit * se, Upper = x_bar + t_crit * se)
}

# Correlation heatmap
create_cor_heatmap <- function(cor_mat, title = "Correlation Matrix", filename = NULL) {
  require(corrplot)
  if (!is.null(filename)) png(filename, width = 1200, height = 1000, res = 120)
  corrplot(cor_mat, method = "color", type = "lower", order = "hclust",
           tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.6,
           title = title, mar = c(0, 0, 2, 0))
  if (!is.null(filename)) { dev.off(); cat("Saved:", filename, "\n") }
}

# MVN Q-Q plot
create_mvn_qq <- function(X, title = "MVN Q-Q Plot", filename = NULL) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  mahal <- calc_mahalanobis(X)
  if (!is.null(filename)) png(filename, width = 800, height = 600, res = 120)
  qqplot(qchisq(ppoints(n), df = p), mahal, main = title,
         xlab = "Chi-squared Quantiles", ylab = "Mahalanobis Distances", pch = 19, col = "steelblue")
  abline(0, 1, col = "red", lwd = 2)
  if (!is.null(filename)) { dev.off(); cat("Saved:", filename, "\n") }
}

# Eta-squared
calc_eta_squared <- function(aov_result) {
  ss <- summary(aov_result)[[1]][, "Sum Sq"]
  eta_sq <- ss / sum(ss)
  names(eta_sq) <- rownames(summary(aov_result)[[1]])
  eta_sq
}

interpret_eta_sq <- function(eta_sq) {
  if (eta_sq < 0.01) return("Negligible")
  if (eta_sq < 0.06) return("Small")
  if (eta_sq < 0.14) return("Medium")
  "Large"
}

# Optimal k for clustering
find_optimal_k <- function(X, max_k = 10) {
  require(cluster)
  dist_mat <- dist(X)
  data.frame(
    k = 2:max_k,
    WSS = sapply(2:max_k, function(k) kmeans(X, k, nstart = 25)$tot.withinss),
    Silhouette = sapply(2:max_k, function(k) mean(silhouette(kmeans(X, k, nstart = 25)$cluster, dist_mat)[, 3]))
  )
}

# Cluster profile
cluster_profile <- function(data, cluster_var, vars) {
  data %>%
    group_by(!!sym(cluster_var)) %>%
    summarise(n = n(), across(all_of(vars), list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE))))
}

# Compare classifiers
compare_classifiers <- function(actual, predictions) {
  require(caret)
  results <- data.frame(Model = names(predictions))
  for (i in seq_along(predictions)) {
    cm <- confusionMatrix(predictions[[i]], actual)
    results$Accuracy[i] <- cm$overall["Accuracy"]
    results$Kappa[i] <- cm$overall["Kappa"]
  }
  results
}

# Standardize
standardize <- function(X) {
  X <- as.matrix(X)
  X_std <- scale(X)
  attr(X_std, "center") <- colMeans(X, na.rm = TRUE)
  attr(X_std, "scale") <- apply(X, 2, sd, na.rm = TRUE)
  X_std
}

unstandardize <- function(X_std) {
  sweep(sweep(X_std, 2, attr(X_std, "scale"), "*"), 2, attr(X_std, "center"), "+")
}

# Print helpers
print_section <- function(title) {
  cat("\n", rep("=", 70), "\n", format(title, width = 70, justify = "centre"), "\n", rep("=", 70), "\n\n", sep = "")
}

print_subsection <- function(title) cat("\n=== ", title, " ===\n", sep = "")

# Load packages
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  cat("All packages loaded.\n")
}

cat("=== STAT 467 Utils Loaded ===\n")

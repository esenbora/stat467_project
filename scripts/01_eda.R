# ============================================================================
# STAT 467 - EXPLORATORY DATA ANALYSIS (EDA)
# ============================================================================
# Purpose: Explore data distributions, correlations, and multivariate normality
#          before applying parametric multivariate methods
#
# Input:   data_country_level.csv (country-level aggregated data with Continent)
#
# Output:  - Correlation heatmaps (figures/correlation_*.png)
#          - Boxplots by Status and Continent (figures/boxplots_*.png)
#          - Scatter plots (figures/scatter_life_expectancy.png)
#          - Scatterplot matrix (figures/scatterplot_matrix.png)
#          - MVN Q-Q plot (figures/mvn_qq_plot.png)
#          - Outlier detection plot (figures/mvn_outlier_detection.png)
#
# Focus:   Understanding Developed vs Developing and Continent differences
#          Assessing MVN assumption critical for Hotelling's T², MANOVA, LDA
#
# Key Concepts:
#   - Multivariate normality (MVN) is required for many parametric methods
#   - MVN tests: Mardia (skewness/kurtosis), Royston, Henze-Zirkler
#   - Each test has different sensitivity and sample size requirements
#
# Dependencies: tidyverse, corrplot, ggcorrplot, GGally, gridExtra, MVN, rstatix, car
# ============================================================================

library(tidyverse)
library(corrplot)
library(ggcorrplot)
library(GGally)
library(gridExtra)
library(MVN)
library(rstatix)
library(car)

# Load country-level data (aggregated means across years 2000-2015)
df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)

# Set grouping variables as factors
df$Status <- as.factor(df$Status)
df$Continent <- as.factor(df$Continent)

cat("Data:", nrow(df), "countries,", ncol(df), "variables\n\n")
cat("=== STATUS DISTRIBUTION ===\n")
print(table(df$Status))
cat("\n=== CONTINENT DISTRIBUTION ===\n")
print(table(df$Continent))
cat("\n")

# Define key variables for analysis
# These represent health outcomes, mortality, vaccination, and socioeconomic factors
key_vars <- c("Life_expectancy", "Adult_Mortality", "infant_deaths", "Alcohol",
              "Hepatitis_B", "BMI", "under_five_deaths", "Polio", "Total_expenditure",
              "Diphtheria", "HIV_AIDS", "thinness_10_19_years", "thinness_5_9_years",
              "Income_composition", "Schooling", "GDP_log", "Population_log")

# ============================================================================
# SECTION 1: Descriptive Statistics
# ============================================================================
cat("=== Descriptive Statistics ===\n")
print(summary(df[, key_vars]))

# Summary statistics by development status
summary_by_status <- df %>%
  group_by(Status) %>%
  summarise(n = n(),
            Life_exp = mean(Life_expectancy),
            Mortality = mean(Adult_Mortality),
            Schooling = mean(Schooling),
            GDP = mean(GDP_log))
cat("\n--- Summary by Status ---\n")
print(summary_by_status)

# Summary statistics by Continent
summary_by_continent <- df %>%
  group_by(Continent) %>%
  summarise(n = n(),
            Life_exp = round(mean(Life_expectancy), 1),
            Mortality = round(mean(Adult_Mortality), 1),
            Schooling = round(mean(Schooling), 1),
            GDP = round(mean(GDP_log), 2))
cat("\n--- Summary by Continent ---\n")
print(summary_by_continent, n = 6)
write.csv(summary_by_continent, "figures/summary_by_continent.csv", row.names = FALSE)

# ============================================================================
# SECTION 2: Correlation Analysis
# ============================================================================
# Calculate pairwise correlations using complete observations
# High correlations (|r| > 0.7) suggest potential multicollinearity issues
# for regression and may indicate latent factors for factor analysis
cor_matrix <- cor(df[, key_vars], use = "pairwise.complete.obs")
write.csv(cor_matrix, "figures/correlation_matrix.csv")

# 2.1 Correlation heatmap with hierarchical clustering
png("figures/correlation_heatmap.png", width = 1400, height = 1200, res = 120)
corrplot(cor_matrix, method = "color", type = "lower", order = "hclust",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.5)
dev.off()
cat("Saved: figures/correlation_heatmap.png\n")

# 2.2 Correlation heatmap with numbers only (recitation style)
png("figures/correlation_heatmap_numbers.png", width = 1400, height = 1200, res = 120)
corrplot(cor_matrix, method = "number", type = "full",
         tl.col = "black", tl.srt = 45, tl.cex = 0.7,
         number.cex = 0.6, cl.cex = 0.8)
dev.off()
cat("Saved: figures/correlation_heatmap_numbers.png\n")

# 2.3 Correlation lower triangle only (recitation style)
png("figures/correlation_heatmap_lower.png", width = 1400, height = 1200, res = 120)
corrplot(cor_matrix, method = "circle", type = "lower",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         order = "hclust", diag = FALSE)
dev.off()
cat("Saved: figures/correlation_heatmap_lower.png\n")

# 2.4 Mixed visualization (upper = circle, lower = number)
png("figures/correlation_heatmap_mixed.png", width = 1400, height = 1200, res = 120)
corrplot.mixed(cor_matrix, lower = "number", upper = "circle",
               tl.col = "black", tl.cex = 0.7, number.cex = 0.6,
               order = "hclust")
dev.off()
cat("Saved: figures/correlation_heatmap_mixed.png\n")

# Identify highly correlated variable pairs
high_cor <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("\nHigh correlations (|r| > 0.7):\n")
  for (i in 1:nrow(high_cor)) {
    if (high_cor[i,1] < high_cor[i,2]) {
      cat(rownames(cor_matrix)[high_cor[i,1]], "-",
          colnames(cor_matrix)[high_cor[i,2]], ":",
          round(cor_matrix[high_cor[i,1], high_cor[i,2]], 3), "\n")
    }
  }
}

# ============================================================================
# SECTION 3: Boxplots by Development Status and Continent
# ============================================================================
# Boxplots reveal:
#   - Central tendency differences between groups
#   - Spread/variability differences
#   - Potential outliers
#   - Distribution shape (symmetry)
selected_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
                   "Income_composition", "HIV_AIDS", "BMI", "Alcohol")

# 3.1 Boxplots by Status
p_list_status <- lapply(selected_vars, function(var) {
  ggplot(df, aes_string(x = "Status", y = var, fill = "Status")) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
    labs(title = var) +
    theme_minimal() +
    theme(legend.position = "none")
})

p_boxplots_status <- do.call(grid.arrange, c(p_list_status, ncol = 4))
ggsave("figures/boxplots_by_status.png", p_boxplots_status, width = 14, height = 8, dpi = 150)
cat("Saved: figures/boxplots_by_status.png\n")

# 3.2 Boxplots by Continent (6 groups)
continent_colors <- c("Africa" = "#E41A1C", "Asia" = "#377EB8",
                      "Europe" = "#4DAF4A", "North America" = "#984EA3",
                      "South America" = "#FF7F00", "Oceania" = "#FFFF33")

p_list_continent <- lapply(selected_vars, function(var) {
  ggplot(df, aes_string(x = "Continent", y = var, fill = "Continent")) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    scale_fill_manual(values = continent_colors) +
    labs(title = var, x = "") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

p_boxplots_continent <- do.call(grid.arrange, c(p_list_continent, ncol = 4))
ggsave("figures/boxplots_by_continent.png", p_boxplots_continent, width = 16, height = 10, dpi = 150)
cat("Saved: figures/boxplots_by_continent.png\n")

# ============================================================================
# SECTION 4: Scatter Plots - Life Expectancy Relationships
# ============================================================================
# Scatter plots show bivariate relationships with Life_expectancy
# Separate regression lines by Status reveal if relationships differ by group
p1 <- ggplot(df, aes(x = Schooling, y = Life_expectancy, color = Status)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()

p2 <- ggplot(df, aes(x = GDP_log, y = Life_expectancy, color = Status)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()

p3 <- ggplot(df, aes(x = Adult_Mortality, y = Life_expectancy, color = Status)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()

p4 <- ggplot(df, aes(x = Income_composition, y = Life_expectancy, color = Status)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()

ggsave("figures/scatter_life_expectancy.png", grid.arrange(p1, p2, p3, p4, ncol = 2),
       width = 12, height = 10, dpi = 150)

# ============================================================================
# SECTION 5: Pairs Plot (Scatterplot Matrix)
# ============================================================================
# Comprehensive view of all pairwise relationships
# Upper triangle: correlation coefficients
# Diagonal: density plots by group
# Lower triangle: scatter plots with points
pairs_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log", "Income_composition", "HIV_AIDS")

# 5.1 GGally pairs plot by Status
p_pairs <- ggpairs(df, columns = pairs_vars, mapping = aes(color = Status, alpha = 0.5),
                   upper = list(continuous = wrap("cor", size = 3)),
                   lower = list(continuous = wrap("points", alpha = 0.5, size = 1))) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()
ggsave("figures/pairs_plot.png", p_pairs, width = 14, height = 12, dpi = 150)
cat("Saved: figures/pairs_plot.png\n")

# 5.2 GGally pairs plot by Continent
p_pairs_continent <- ggpairs(df, columns = pairs_vars,
                              mapping = aes(color = Continent, alpha = 0.5),
                              upper = list(continuous = wrap("cor", size = 2.5)),
                              lower = list(continuous = wrap("points", alpha = 0.4, size = 0.8))) +
  scale_color_manual(values = continent_colors) +
  theme_minimal()
ggsave("figures/pairs_plot_continent.png", p_pairs_continent, width = 16, height = 14, dpi = 150)
cat("Saved: figures/pairs_plot_continent.png\n")

# 5.3 car::scatterplotMatrix (from recitation guidelines)
# This provides histograms on diagonal and regression lines
scatterplot_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log")

png("figures/scatterplot_matrix.png", width = 1200, height = 1200, res = 120)
car::scatterplotMatrix(df[, scatterplot_vars],
                       diagonal = "histogram",
                       smooth = FALSE,
                       regLine = list(col = "red", lwd = 2),
                       pch = 19, cex = 0.5,
                       col = "steelblue",
                       main = "Scatterplot Matrix with Histograms")
dev.off()
cat("Saved: figures/scatterplot_matrix.png\n")

# 5.4 car::scatterplotMatrix by Status
png("figures/scatterplot_matrix_status.png", width = 1200, height = 1200, res = 120)
car::scatterplotMatrix(~ Life_expectancy + Adult_Mortality + Schooling + GDP_log | Status,
                       data = df,
                       diagonal = "histogram",
                       smooth = FALSE,
                       regLine = list(lwd = 2),
                       pch = c(16, 17),
                       col = c("#2E86AB", "#E94F37"),
                       main = "Scatterplot Matrix by Status")
dev.off()
cat("Saved: figures/scatterplot_matrix_status.png\n")

# ============================================================================
# SECTION 6: MULTIVARIATE NORMALITY ASSESSMENT
# ============================================================================
# WHY MVN MATTERS:
#   - Hotelling's T², MANOVA, LDA assume multivariate normality
#   - Violations can inflate Type I error rates
#   - However, large samples provide some robustness (CLT)
#
# MVN TEST SELECTION (Güler, 2025 guidelines):
#   - Mardia: Good for detecting skewness/kurtosis, n > 20
#   - Royston: Extension of Shapiro-Wilk, n < 5000, sensitive to outliers
#   - Henze-Zirkler: Consistent power across sample sizes, based on moment generating function
#
# INTERPRETATION:
#   - If p < 0.05: Reject H0 of MVN (evidence against normality)
#   - If p >= 0.05: Fail to reject MVN (no evidence against normality)
#   - Note: Large samples may reject MVN even for minor deviations

norm_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
               "Income_composition", "BMI", "HIV_AIDS")
X_norm <- df[, norm_vars]
n <- nrow(X_norm)

cat("\n=== MULTIVARIATE NORMALITY ASSESSMENT ===\n")
cat("Variables:", paste(norm_vars, collapse = ", "), "\n")
cat("Sample size: n =", n, "\n\n")

# ----------------------------------------------------------------------------
# 6.1 Univariate Normality Tests (Shapiro-Wilk)
# ----------------------------------------------------------------------------
# Check each variable individually first
# MVN requires ALL marginal distributions to be normal (necessary but not sufficient)
cat("--- 1. Univariate Normality (Shapiro-Wilk) ---\n")
univar_results <- X_norm %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  shapiro_test(Value)
univar_results$Normal <- ifelse(univar_results$p > 0.05, "Yes", "No")
print(univar_results)

# Count violations - many violations suggest MVN is unlikely
n_violations <- sum(univar_results$p < 0.05)
cat("\nUnivariate normality violations:", n_violations, "of", length(norm_vars), "variables\n")

# Univariate histograms for visual inspection
# Note: univariatePlot argument may not be available in all MVN versions
tryCatch({
  png("figures/eda_univariate_histograms.png", width = 1200, height = 800, res = 120)
  mvn_univ <- mvn(X_norm, univariatePlot = "histogram")
  dev.off()
}, error = function(e) {
  dev.off()
  cat("Note: Univariate histogram plot skipped (MVN version compatibility)\n")
})

# ----------------------------------------------------------------------------
# 6.2 Multivariate Normality Tests
# ----------------------------------------------------------------------------
# Three complementary tests provide robust assessment
cat("\n--- 2. Multivariate Normality Tests ---\n")

# Mardia's test: Tests multivariate skewness and kurtosis
# - Skewness: Measures asymmetry of the joint distribution
# - Kurtosis: Measures tail heaviness (outlier proneness)
mvn_mardia <- mvn(X_norm, mvn_test = "mardia")
cat("\nMardia's Test:\n")
print(mvn_mardia$multivariate_normality)

# Royston's test: Extension of Shapiro-Wilk to multivariate case
# - Based on Shapiro-Wilk statistics for each variable
# - Combines using a normalized statistic
mvn_royston <- mvn(X_norm, mvn_test = "royston")
cat("\nRoyston's Test:\n")
print(mvn_royston$multivariate_normality)

# Henze-Zirkler test: Based on the moment generating function
# - Compares empirical and theoretical moment generating functions
# - No sample size limitations, consistent power
mvn_hz <- mvn(X_norm, mvn_test = "hz")
cat("\nHenze-Zirkler Test:\n")
print(mvn_hz$multivariate_normality)

# Summary of all MVN tests
cat("\n--- MVN Test Summary ---\n")
mvn_summary <- data.frame(
  Test = c("Mardia Skewness", "Mardia Kurtosis", "Royston", "Henze-Zirkler"),
  p_value = c(
    mvn_mardia$multivariate_normality$p.value[1],
    mvn_mardia$multivariate_normality$p.value[2],
    mvn_royston$multivariate_normality$p.value,
    mvn_hz$multivariate_normality$p.value
  )
)
mvn_summary$Conclusion <- ifelse(mvn_summary$p_value < 0.05,
                                  "Reject MVN", "Fail to reject MVN")
print(mvn_summary)

# ----------------------------------------------------------------------------
# 6.3 Chi-Square Q-Q Plot
# ----------------------------------------------------------------------------
# Under MVN, Mahalanobis distances follow chi-squared distribution with p df
# Q-Q plot compares observed distances to theoretical quantiles
# Points should fall along the diagonal line if MVN holds
cat("\n--- 3. Chi-Square Q-Q Plot ---\n")
png("figures/mvn_qq_plot.png", width = 800, height = 600, res = 120)
result_qq <- mvn(X_norm, mvn_test = "mardia")
plot(result_qq, diagnostic = "multivariate", type = "qq")
dev.off()
cat("Saved: figures/mvn_qq_plot.png\n")

# ----------------------------------------------------------------------------
# 6.4 Multivariate Outlier Detection
# ----------------------------------------------------------------------------
# Outliers can severely distort MVN tests and parameter estimates
# Two methods:
#   - Classical Mahalanobis: Uses sample mean and covariance
#   - Adjusted Mahalanobis: Uses robust (outlier-resistant) estimates
cat("\n--- 4. Multivariate Outlier Detection ---\n")

# Try outlier detection with tryCatch since covariance may be singular
# (e.g., when HIV_AIDS has many identical values)
tryCatch({
  # Mahalanobis distance based outliers (classical)
  mvn_outlier_mah <- mvn(X_norm, mvn_test = "mardia",
                         multivariate_outlier_method = "quan")
  n_outliers_mah <- sum(mvn_outlier_mah$multivariate_outliers$Outlier == "TRUE",
                        na.rm = TRUE)
  cat("Outliers (Mahalanobis distance, alpha = 0.025):", n_outliers_mah, "\n")

  # Adjusted Mahalanobis distance (robust to masking effects)
  mvn_outlier_adj <- mvn(X_norm, mvn_test = "mardia",
                         multivariate_outlier_method = "adj")
  n_outliers_adj <- sum(mvn_outlier_adj$multivariate_outliers$Outlier == "TRUE",
                        na.rm = TRUE)
  cat("Outliers (Adjusted Mahalanobis):", n_outliers_adj, "\n")

  # Save outlier detection plot
  png("figures/mvn_outlier_detection.png", width = 900, height = 600, res = 120)
  plot(mvn_outlier_mah, diagnostic = "outlier")
  dev.off()
  cat("Saved: figures/mvn_outlier_detection.png\n")

  # List outlier countries for domain interpretation
  if (n_outliers_mah > 0) {
    outlier_indices <- which(
      mvn_outlier_mah$multivariate_outliers$Outlier == "TRUE")
    cat("\nOutlier countries (Mahalanobis):\n")
    print(df$Country[outlier_indices])
  }
}, error = function(e) {
  cat("Note: Outlier detection skipped - covariance matrix is singular.\n")
  cat("This can happen when variables have limited variability.\n")
})

# ----------------------------------------------------------------------------
# 6.5 Implications for Downstream Analyses
# ----------------------------------------------------------------------------
cat("\n--- 5. Implications for Downstream Analyses ---\n")
# Check if any MVN test rejected the null hypothesis
# (p_value may be a string like "<0.001", so check Conclusion column instead)
if (any(grepl("Reject", mvn_summary$Conclusion))) {
  cat("WARNING: MVN assumption violated.\n")
  cat("Recommendations:\n")
  cat("  - Use robust methods where available\n")
  cat("  - Consider transformations for skewed variables\n")
  cat("  - Interpret parametric tests with caution\n")
  cat("  - Large sample size (n=", n, ") provides some robustness\n")
} else {
  cat("MVN assumption satisfied - parametric methods appropriate\n")
}

cat("\n=== EDA Complete ===\n")
cat("Figures saved. Run 02_mean_vector_inference.R next.\n")

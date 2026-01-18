# STAT 467 - Exploratory Data Analysis
# Input: data_country_level.csv

library(tidyverse)
library(corrplot)
library(ggcorrplot)
library(GGally)
library(gridExtra)
library(MVN)
library(rstatix)  # For univariate normality tests

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries,", ncol(df), "variables\n")
cat("Status:", table(df$Status), "\n\n")

key_vars <- c("Life_expectancy", "Adult_Mortality", "infant_deaths", "Alcohol",
              "Hepatitis_B", "BMI", "under_five_deaths", "Polio", "Total_expenditure",
              "Diphtheria", "HIV_AIDS", "thinness_10_19_years", "thinness_5_9_years",
              "Income_composition", "Schooling", "GDP_log", "Population_log")

# Descriptive statistics
cat("=== Descriptive Statistics ===\n")
print(summary(df[, key_vars]))

# Summary by Status
summary_by_status <- df %>%
  group_by(Status) %>%
  summarise(n = n(),
            Life_exp = mean(Life_expectancy),
            Mortality = mean(Adult_Mortality),
            Schooling = mean(Schooling),
            GDP = mean(GDP_log))
print(summary_by_status)

# Correlation matrix
cor_matrix <- cor(df[, key_vars], use = "pairwise.complete.obs")
write.csv(cor_matrix, "figures/correlation_matrix.csv")

png("figures/correlation_heatmap.png", width = 1400, height = 1200, res = 120)
corrplot(cor_matrix, method = "color", type = "lower", order = "hclust",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.5)
dev.off()

# High correlations
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

# Boxplots by Status
selected_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
                   "Income_composition", "HIV_AIDS", "BMI", "Alcohol")

p_list <- lapply(selected_vars, function(var) {
  ggplot(df, aes_string(x = "Status", y = var, fill = "Status")) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2 , alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
    labs(title = var) +
    theme_minimal() +
    theme(legend.position = "none")
})

p_boxplots <- do.call(grid.arrange, c(p_list, ncol = 4))
ggsave("figures/boxplots_by_status.png", p_boxplots, width = 14, height = 8, dpi = 150)

# Scatter plots - Life expectancy relationships
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

# Pairs plot
pairs_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log", "Income_composition", "HIV_AIDS")
p_pairs <- ggpairs(df, columns = pairs_vars, mapping = aes(color = Status, alpha = 0.5),
                   upper = list(continuous = wrap("cor", size = 3)),
                   lower = list(continuous = wrap("points", alpha = 0.5, size = 1))) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()
ggsave("figures/pairs_plot.png", p_pairs, width = 14, height = 12, dpi = 150)

# ============================================================================
# MULTIVARIATE NORMALITY ASSESSMENT (Following MVN Package Guidelines)
# ============================================================================
norm_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
               "Income_composition", "BMI", "HIV_AIDS")
X_norm <- df[, norm_vars]
n <- nrow(X_norm)

cat("\n=== MULTIVARIATE NORMALITY ASSESSMENT ===\n")
cat("Variables:", paste(norm_vars, collapse = ", "), "\n")
cat("Sample size: n =", n, "\n\n")

# 1. UNIVARIATE NORMALITY TESTS (Shapiro-Wilk)
cat("--- 1. Univariate Normality (Shapiro-Wilk) ---\n")
univar_results <- X_norm %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  shapiro_test(Value)
univar_results$Normal <- ifelse(univar_results$p > 0.05, "Yes", "No")
print(univar_results)

# Count violations
n_violations <- sum(univar_results$p < 0.05)
cat("\nUnivariate normality violations:", n_violations, "of", length(norm_vars), "variables\n")

# Univariate histograms
png("figures/eda_univariate_histograms.png", width = 1200, height = 800, res = 120)
mvn_univ <- mvn(X_norm, univariatePlot = "histogram")
dev.off()

# 2. MULTIVARIATE NORMALITY TESTS
# Test selection based on sample size (Güler, 2026 recommendations):
# - Mardia: n < 20
# - Royston: n < 5000 (based on Shapiro-Wilk)
# - Henze-Zirkler: Any sample size

cat("\n--- 2. Multivariate Normality Tests ---\n")

# Mardia's test (skewness and kurtosis)
mvn_mardia <- mvn(X_norm, mvnTest = "mardia")
cat("\nMardia's Test:\n")
print(mvn_mardia$multivariateNormality)

# Royston's test (extension of Shapiro-Wilk)
mvn_royston <- mvn(X_norm, mvnTest = "royston")
cat("\nRoyston's Test:\n")
print(mvn_royston$multivariateNormality)

# Henze-Zirkler test (no sample size limitation)
mvn_hz <- mvn(X_norm, mvnTest = "hz")
cat("\nHenze-Zirkler Test:\n")
print(mvn_hz$multivariateNormality)

# Summary
cat("\n--- MVN Test Summary ---\n")
mvn_summary <- data.frame(
  Test = c("Mardia Skewness", "Mardia Kurtosis", "Royston", "Henze-Zirkler"),
  p_value = c(
    mvn_mardia$multivariateNormality$`p value`[1],
    mvn_mardia$multivariateNormality$`p value`[2],
    mvn_royston$multivariateNormality$`p value`,
    mvn_hz$multivariateNormality$`p value`
  )
)
mvn_summary$Conclusion <- ifelse(as.numeric(mvn_summary$p_value) < 0.05,
                                  "Reject MVN", "Fail to reject MVN")
print(mvn_summary)

# 3. CHI-SQUARE Q-Q PLOT
cat("\n--- 3. Chi-Square Q-Q Plot ---\n")
png("figures/mvn_qq_plot.png", width = 800, height = 600, res = 120)
result_qq <- mvn(X_norm, mvnTest = "mardia")
plot(result_qq, diagnostic = "multivariate", type = "qq")
dev.off()
cat("Saved: figures/mvn_qq_plot.png\n")

# 4. MULTIVARIATE OUTLIER DETECTION
cat("\n--- 4. Multivariate Outlier Detection ---\n")

# Mahalanobis distance based outliers
mvn_outlier_mah <- mvn(X_norm, mvnTest = "mardia",
                       multivariateOutlierMethod = "quan")
n_outliers_mah <- sum(mvn_outlier_mah$multivariateOutliers$Outlier == "TRUE",
                      na.rm = TRUE)
cat("Outliers (Mahalanobis distance, alpha = 0.025):", n_outliers_mah, "\n")

# Adjusted Mahalanobis distance
mvn_outlier_adj <- mvn(X_norm, mvnTest = "mardia",
                       multivariateOutlierMethod = "adj")
n_outliers_adj <- sum(mvn_outlier_adj$multivariateOutliers$Outlier == "TRUE",
                      na.rm = TRUE)
cat("Outliers (Adjusted Mahalanobis):", n_outliers_adj, "\n")

# Save outlier plot
png("figures/mvn_outlier_detection.png", width = 900, height = 600, res = 120)
plot(mvn_outlier_mah, diagnostic = "outlier")
dev.off()
cat("Saved: figures/mvn_outlier_detection.png\n")

# List outlier countries
if (n_outliers_mah > 0) {
  outlier_indices <- which(mvn_outlier_mah$multivariateOutliers$Outlier == "TRUE")
  cat("\nOutlier countries (Mahalanobis):\n")
  print(df$Country[outlier_indices])
}

# 5. IMPLICATIONS FOR DOWNSTREAM ANALYSES
cat("\n--- 5. Implications for Downstream Analyses ---\n")
if (any(as.numeric(mvn_summary$p_value) < 0.05)) {
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

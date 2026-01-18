# ============================================================================
# STAT 467 - PRINCIPAL COMPONENT ANALYSIS AND PC REGRESSION
# ============================================================================
# Purpose: Reduce dimensionality while preserving variance; use PCs as
#          predictors in regression to address multicollinearity
#
# Input:   data_country_level.csv (country-level aggregated data)
#
# Output:  - Scree plot (figures/pca_screeplot.png)
#          - Cumulative variance (figures/pca_cumulative_variance.png)
#          - Loadings heatmap (figures/pca_loadings_heatmap.png)
#          - Biplot (figures/pca_biplot.png)
#          - PC individuals plot (figures/pca_individuals.png)
#          - PC contributions (figures/pca_contributions.png)
#          - PCR diagnostics (figures/pcr_diagnostics.png)
#          - Cross-validation results (figures/pcr_cv_results.csv)
#
# Methods:
#   - PCA via singular value decomposition (prcomp)
#   - KMO and Bartlett's test for sampling adequacy
#   - Principal Component Regression (PCR)
#   - 5-fold cross-validation for model selection
#
# Key Concepts:
#   - PCA finds orthogonal linear combinations that maximize variance
#   - First PC captures most variance, each subsequent PC captures remaining
#   - Loadings: correlations between original variables and PCs
#   - PCR uses PCs as predictors, eliminating multicollinearity
#
# Component Retention Criteria:
#   1. Kaiser criterion: eigenvalue > 1 (often retains too many)
#   2. Variance explained: cumulative variance > 70-80%
#   3. Scree plot: "elbow" where eigenvalues level off
#   4. Parallel analysis: compare to random data (for FA)
#
# Dependencies: tidyverse, FactoMineR, factoextra, corrplot, gridExtra,
#               caret, car, psych
# ============================================================================

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(gridExtra)
library(caret)
library(car)
library(psych)

# Load data
df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n\n")

# Select variables for PCA
# Exclude response variable (Life_expectancy) - it will be predicted later
pca_vars <- c("Adult_Mortality", "infant_deaths", "Alcohol", "Hepatitis_B", "BMI",
              "under_five_deaths", "Polio", "Total_expenditure", "Diphtheria", "HIV_AIDS",
              "thinness_10_19_years", "thinness_5_9_years", "Income_composition",
              "Schooling", "GDP_log", "Population_log")

# Prepare data
X_pca <- df %>% select(all_of(pca_vars)) %>% drop_na()
df_complete <- df %>% select(Country, Status, Life_expectancy, all_of(pca_vars)) %>% drop_na()

cat("Variables:", length(pca_vars), "\n")
cat("Complete obs:", nrow(X_pca), "\n\n")

# ============================================================================
# SECTION 1: SAMPLING ADEQUACY TESTS
# ============================================================================
# Before PCA, check if correlation structure supports dimension reduction

# Kaiser-Meyer-Olkin (KMO) Measure of Sampling Adequacy
# - Measures proportion of variance that might be common variance
# - KMO > 0.8: meritorious; > 0.6: mediocre; < 0.5: unacceptable
cor_matrix <- cor(X_pca)
kmo <- KMO(cor_matrix)
cat("KMO:", round(kmo$MSA, 3), "\n")

# Bartlett's Test of Sphericity
# - Tests H₀: correlation matrix = identity (no correlations)
# - Significant p-value indicates correlations exist (PCA appropriate)
bartlett <- cortest.bartlett(cor_matrix, n = nrow(X_pca))
cat("Bartlett p-value:", format(bartlett$p.value, scientific = TRUE), "\n\n")

# ============================================================================
# SECTION 2: PRINCIPAL COMPONENT ANALYSIS
# ============================================================================
# PCA Algorithm:
#   1. Standardize data (center and scale)
#   2. Compute covariance/correlation matrix
#   3. Find eigenvalues and eigenvectors
#   4. Eigenvalues = variance explained by each PC
#   5. Eigenvectors = loadings (weights for linear combinations)
#
# prcomp uses SVD (numerically stable):
#   X = UDV', where columns of V are eigenvectors (loadings)

pca_result <- prcomp(X_pca, center = TRUE, scale. = TRUE)

# Extract eigenvalues and variance explained
eigenvalues <- pca_result$sdev^2          # Eigenvalue = variance of each PC
var_exp <- eigenvalues / sum(eigenvalues) # Proportion of variance
cum_var <- cumsum(var_exp)                # Cumulative variance

# Summary table
pca_summary <- data.frame(
  PC = paste0("PC", 1:length(eigenvalues)),
  Eigenvalue = round(eigenvalues, 4),
  Var_Pct = round(var_exp * 100, 2),
  Cum_Pct = round(cum_var * 100, 2)
)
cat("=== EIGENVALUES ===\n")
print(pca_summary)
write.csv(pca_summary, "figures/pca_eigenvalues.csv", row.names = FALSE)

# ============================================================================
# SECTION 3: COMPONENT RETENTION
# ============================================================================
# Multiple criteria to decide how many PCs to retain

# Kaiser criterion: eigenvalue > 1
# Rationale: Each PC should explain more variance than a single original variable
n_kaiser <- sum(eigenvalues > 1)

# Variance explained criterion
n_70var <- which(cum_var >= 0.70)[1]  # First PC reaching 70%
n_80var <- which(cum_var >= 0.80)[1]  # First PC reaching 80%

cat("\nKaiser criterion:", n_kaiser, "components\n")
cat("70% variance:", n_70var, "components\n")
cat("80% variance:", n_80var, "components\n")

# Use 70% variance as primary criterion
# Kaiser often retains too many components
n_components <- n_70var
cat("\nUsing", n_components, "components (70% variance criterion)\n")
cat("Note: Kaiser suggests", n_kaiser, "(reference only)\n")

# ============================================================================
# SECTION 4: VISUALIZATION
# ============================================================================

# ----------------------------------------------------------------------------
# 4.1 Scree Plot
# ----------------------------------------------------------------------------
# Shows eigenvalues in decreasing order
# Look for "elbow" where eigenvalues level off
# Horizontal line = Kaiser criterion (eigenvalue = 1)
p_scree <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50),
                    barfill = "steelblue", linecolor = "red") +
  geom_hline(yintercept = 100/length(eigenvalues), linetype = "dashed") +
  labs(title = "Scree Plot") + theme_minimal()
ggsave("figures/pca_screeplot.png", p_scree, width = 10, height = 6, dpi = 150)

# ----------------------------------------------------------------------------
# 4.2 Cumulative Variance Plot
# ----------------------------------------------------------------------------
# Shows how variance accumulates with each added PC
# Horizontal lines at 70% and 80% thresholds
p_cumvar <- ggplot(pca_summary, aes(x = 1:nrow(pca_summary), y = Cum_Pct)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 3) +
  geom_hline(yintercept = c(70, 80), linetype = "dashed", color = c("orange", "red")) +
  labs(title = "Cumulative Variance", x = "Components", y = "Cumulative %") +
  theme_minimal()
ggsave("figures/pca_cumulative_variance.png", p_cumvar, width = 10, height = 6, dpi = 150)

# ============================================================================
# SECTION 5: LOADINGS ANALYSIS
# ============================================================================
# Loadings: correlations between original variables and PCs
# - High positive loading: variable increases with PC
# - High negative loading: variable decreases with PC
# - Near-zero loading: variable unrelated to PC
#
# Interpretation: loadings help name/interpret each PC

loadings <- pca_result$rotation[, 1:n_components]
cat("\n=== LOADINGS (first", n_components, "PCs) ===\n")
print(round(loadings, 3))
write.csv(round(loadings, 4), "figures/pca_loadings.csv")

# Loadings heatmap
loadings_df <- as.data.frame(loadings) %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = starts_with("PC"), names_to = "Component", values_to = "Loading")

p_loadings <- ggplot(loadings_df, aes(x = Component, y = Variable, fill = Loading)) +
  geom_tile() + geom_text(aes(label = round(Loading, 2)), size = 3) +
  scale_fill_gradient2(low = "#E94F37", mid = "white", high = "#2E86AB", limits = c(-1, 1)) +
  labs(title = "PCA Loadings Heatmap") + theme_minimal()
ggsave("figures/pca_loadings_heatmap.png", p_loadings, width = 10, height = 10, dpi = 150)

# ----------------------------------------------------------------------------
# 5.1 Biplot
# ----------------------------------------------------------------------------
# Shows both observations (points) and variables (arrows) in PC space
# Arrow direction: variable's correlation with PCs
# Arrow length: strength of relationship
# Angle between arrows: correlation between variables (small = high correlation)
p_biplot <- fviz_pca_biplot(pca_result, geom.ind = "point", col.ind = df_complete$Status,
                             palette = c("#2E86AB", "#E94F37"), addEllipses = TRUE,
                             ellipse.level = 0.95, repel = TRUE, legend.title = "Status") +
  labs(title = "PCA Biplot") + theme_minimal()
ggsave("figures/pca_biplot.png", p_biplot, width = 14, height = 10, dpi = 150)

# ----------------------------------------------------------------------------
# 5.2 Individual Countries in PC Space
# ----------------------------------------------------------------------------
pc_scores <- as.data.frame(pca_result$x[, 1:n_components])
pc_scores$Status <- df_complete$Status

p_ind <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(alpha = 0.7, size = 2.5) + stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = "Countries in PC Space",
       x = paste0("PC1 (", round(var_exp[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2]*100, 1), "%)")) +
  theme_minimal()
ggsave("figures/pca_individuals.png", p_ind, width = 10, height = 8, dpi = 150)

# ----------------------------------------------------------------------------
# 5.3 Variable Contributions
# ----------------------------------------------------------------------------
# Shows which variables contribute most to each PC
p_contrib1 <- fviz_contrib(pca_result, choice = "var", axes = 1, top = 10) + theme_minimal()
p_contrib2 <- fviz_contrib(pca_result, choice = "var", axes = 2, top = 10) + theme_minimal()
ggsave("figures/pca_contributions.png",
       grid.arrange(p_contrib1, p_contrib2, ncol = 2), width = 14, height = 6, dpi = 150)

# ============================================================================
# SECTION 6: PRINCIPAL COMPONENT REGRESSION (PCR)
# ============================================================================
# PCR uses PC scores as predictors instead of original variables
# Advantages:
#   - Eliminates multicollinearity (PCs are orthogonal)
#   - Reduces number of predictors (dimension reduction)
#   - Can improve prediction when original variables are highly correlated
#
# Model: Y = β₀ + β₁PC₁ + β₂PC₂ + ... + βₖPCₖ + ε

cat("\n=== PC REGRESSION ===\n")
pc_data <- as.data.frame(pca_result$x)
pc_data$Life_expectancy <- df_complete$Life_expectancy
pc_data$Status <- df_complete$Status

# Fit models with different numbers of PCs
pcr_m3 <- lm(Life_expectancy ~ PC1 + PC2 + PC3, data = pc_data)
pcr_m5 <- lm(Life_expectancy ~ PC1 + PC2 + PC3 + PC4 + PC5, data = pc_data)

cat("Model (3 PCs): R² =", round(summary(pcr_m3)$r.squared, 4), "\n")
cat("Model (5 PCs): R² =", round(summary(pcr_m5)$r.squared, 4), "\n")

cat("\n=== 5-PC MODEL SUMMARY ===\n")
print(summary(pcr_m5))

# Save coefficients
write.csv(data.frame(
  Term = names(coef(pcr_m5)),
  Coefficient = round(coef(pcr_m5), 4),
  p_value = summary(pcr_m5)$coefficients[, 4]
), "figures/pcr_coefficients.csv", row.names = FALSE)

# ============================================================================
# SECTION 7: CROSS-VALIDATION
# ============================================================================
# Use k-fold CV to select optimal number of PCs
# Prevents overfitting by evaluating on held-out data

cat("\n=== CROSS-VALIDATION ===\n")
set.seed(123)
train_control <- trainControl(method = "cv", number = 5)

cv_results <- data.frame(n_PCs = 1:10, RMSE = NA, R2 = NA)
for (k in 1:10) {
  formula_k <- as.formula(paste("Life_expectancy ~",
                                paste(paste0("PC", 1:k), collapse = " + ")))
  model_cv <- train(formula_k, data = pc_data, method = "lm",
                    trControl = train_control)
  # Extract CV results
  cv_results$RMSE[k] <- min(model_cv$results$RMSE)
  cv_results$R2[k] <- max(model_cv$results$Rsquared)
}
print(cv_results)
write.csv(cv_results, "figures/pcr_cv_results.csv", row.names = FALSE)

# CV R² plot
p_cv <- ggplot(cv_results, aes(x = n_PCs, y = R2)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 3) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "CV R² by Number of PCs", x = "PCs", y = "R²") + theme_minimal()
ggsave("figures/pcr_cv_r2.png", p_cv, width = 10, height = 6, dpi = 150)

cat("Optimal PCs:", cv_results$n_PCs[which.max(cv_results$R2)], "\n")

# ============================================================================
# SECTION 8: COMPARISON WITH STANDARD REGRESSION
# ============================================================================
# Compare PCR to standard OLS with all original predictors
# Check for multicollinearity using Variance Inflation Factor (VIF)
# VIF > 5 or 10 indicates problematic multicollinearity

std_model <- lm(Life_expectancy ~ .,
                data = df_complete %>% select(Life_expectancy, all_of(pca_vars)))
cat("\nStandard Regression R²:", round(summary(std_model)$r.squared, 4), "\n")

vif_vals <- car::vif(std_model)
cat("Variables with VIF > 5:", sum(vif_vals > 5), "\n")

# ============================================================================
# SECTION 9: REGRESSION DIAGNOSTICS
# ============================================================================
cat("\n=== REGRESSION DIAGNOSTICS ===\n")

# Shapiro-Wilk test for residual normality
resid_5pc <- residuals(pcr_m5)
shapiro_result <- shapiro.test(resid_5pc)
cat("Shapiro-Wilk test for residual normality:\n")
cat("  W =", round(shapiro_result$statistic, 4),
    ", p =", format(shapiro_result$p.value, digits = 3), "\n")
if (shapiro_result$p.value < 0.05) {
  cat("  WARNING: Residuals may not be normally distributed (p < 0.05)\n")
} else {
  cat("  Residual normality assumption is met\n")
}

# Breusch-Pagan test for heteroscedasticity
# H₀: Homoscedasticity (constant variance)
bp_result <- car::ncvTest(pcr_m5)
cat("\nBreusch-Pagan test for heteroscedasticity:\n")
cat("  Chi-sq =", round(bp_result$ChiSquare, 2),
    ", df =", bp_result$Df,
    ", p =", format(bp_result$p, digits = 3), "\n")
if (bp_result$p < 0.05) {
  cat("  WARNING: Heteroscedasticity detected (p < 0.05)\n")
} else {
  cat("  Homoscedasticity assumption is met\n")
}

# Diagnostic plots
p_fitted <- ggplot(pc_data, aes(x = fitted(pcr_m5), y = Life_expectancy)) +
  geom_point(aes(color = Status), alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = "Fitted vs Actual (5-PC Model)") + theme_minimal()

p_resid <- ggplot(data.frame(Fitted = fitted(pcr_m5), Residuals = resid_5pc),
                  aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "orange", linetype = "dashed") +
  labs(title = "Residual Plot") + theme_minimal()

# Q-Q plot for residuals
p_qq <- ggplot(data.frame(Residuals = resid_5pc), aes(sample = Residuals)) +
  stat_qq(color = "steelblue") +
  stat_qq_line(color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Residuals") + theme_minimal()

ggsave("figures/pcr_diagnostics.png",
       grid.arrange(p_fitted, p_resid, p_qq, ncol = 3),
       width = 18, height = 6, dpi = 150)

# ============================================================================
# SECTION 10: ROBUST STANDARD ERRORS
# ============================================================================
# When heteroscedasticity is detected (Breusch-Pagan p < 0.05),
# robust standard errors provide valid inference
# HC3 is recommended for small samples (better finite-sample performance)

cat("\n=== ROBUST STANDARD ERRORS ===\n")

library(sandwich)
library(lmtest)

# Calculate robust standard errors (HC3 - best for small samples)
robust_vcov <- vcovHC(pcr_m5, type = "HC3")
robust_test <- coeftest(pcr_m5, vcov = robust_vcov)

cat("Comparing OLS vs Robust (HC3) Standard Errors:\n\n")

# Create comparison table
ols_se <- summary(pcr_m5)$coefficients[, "Std. Error"]
robust_se <- sqrt(diag(robust_vcov))
ols_p <- summary(pcr_m5)$coefficients[, "Pr(>|t|)"]
robust_p <- robust_test[, "Pr(>|t|)"]

se_comparison <- data.frame(
  Term = names(coef(pcr_m5)),
  Estimate = round(coef(pcr_m5), 4),
  OLS_SE = round(ols_se, 4),
  Robust_SE = round(robust_se, 4),
  SE_Ratio = round(robust_se / ols_se, 3),
  OLS_p = round(ols_p, 4),
  Robust_p = round(robust_p, 4)
)
print(se_comparison)

cat("\nInterpretation:\n")
cat("  SE_Ratio > 1: Robust SE larger (OLS underestimates uncertainty)\n")
cat("  SE_Ratio < 1: Robust SE smaller (OLS overestimates uncertainty)\n")

# Check if conclusions change
sig_changes <- sum((ols_p < 0.05) != (robust_p < 0.05))
if (sig_changes > 0) {
  cat("\n  WARNING:", sig_changes, "coefficient(s) changed significance status\n")
  cat("  with robust SE - results should be interpreted cautiously\n")
} else {
  cat("\n  All coefficients retain significance status with robust SE\n")
  cat("  Results are ROBUST to heteroscedasticity\n")
}

# Save robust results
write.csv(se_comparison, "figures/pcr_robust_se_comparison.csv", row.names = FALSE)

# Full robust coefficient test output
cat("\n=== ROBUST COEFFICIENT TESTS (HC3) ===\n")
print(robust_test)

# ============================================================================
# SECTION 11: STATISTICAL VALIDITY SUMMARY
# ============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("REGRESSION DIAGNOSTICS SUMMARY\n")
cat(strrep("=", 60), "\n")

resid_normality <- ifelse(shapiro_result$p.value >= 0.05, "Met", "Violated")
homoscedasticity <- ifelse(bp_result$p >= 0.05, "Met", "Violated")
robust_needed <- ifelse(bp_result$p < 0.05, "Yes (HC3 SE provided)", "No")

cat("Residual Normality:        ", resid_normality,
    " (Shapiro-Wilk p =", round(shapiro_result$p.value, 3), ")\n")
cat("Homoscedasticity:          ", homoscedasticity,
    " (Breusch-Pagan p =", round(bp_result$p, 3), ")\n")
cat("Robust SE Needed:          ", robust_needed, "\n")
cat("Model R²:                  ", round(summary(pcr_m5)$r.squared, 4), "\n")
cat(strrep("=", 60), "\n")

cat("\n=== PCA and PC Regression Complete ===\n")
cat("Run 05_factor_analysis.R next.\n")

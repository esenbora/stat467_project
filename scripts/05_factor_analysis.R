# ============================================================================
# STAT 467 - EXPLORATORY FACTOR ANALYSIS (EFA)
# ============================================================================
# Purpose: Identify latent constructs (factors) underlying observed variables
#          Understand the structure of health/socioeconomic indicators
#
# Input:   data_country_level.csv (country-level aggregated data)
#
# Output:  - Parallel analysis (figures/fa_parallel_analysis.png)
#          - Loadings heatmap (figures/fa_loadings_heatmap.png)
#          - Factor correlations (figures/fa_factor_correlations.png)
#          - Factor scores plot (figures/fa_factor_scores_plot.png)
#          - Factor loadings (figures/fa_loadings.csv)
#          - Factor scores (figures/fa_factor_scores.csv)
#
# Methods:
#   - Parallel analysis for factor number determination
#   - Principal Axis (PA) or Maximum Likelihood (ML) extraction
#   - Varimax (orthogonal) and Promax (oblique) rotation
#   - Factor score estimation
#
# Key Concepts:
#   - PCA vs Factor Analysis:
#     * PCA: data reduction, maximizes variance, components are linear combos
#     * FA: theory-driven, identifies latent factors causing observed correlations
#   - Common vs. unique variance:
#     * FA models only common variance (shared among variables)
#     * Unique variance = specific + error variance
#
# Factor Number Determination:
#   1. Parallel analysis (gold standard): compare to random data eigenvalues
#   2. Kaiser criterion (eigenvalue > 1): often overestimates
#   3. Scree test: subjective "elbow" interpretation
#
# Rotation Methods:
#   - Varimax: Orthogonal (factors uncorrelated), simpler interpretation
#   - Promax: Oblique (factors can correlate), more realistic
#
# Dependencies: tidyverse, psych, GPArotation, corrplot, gridExtra, nFactors
# ============================================================================

library(tidyverse)
library(psych)
library(GPArotation)
library(corrplot)
library(gridExtra)
library(nFactors)

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# Load data
df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n\n")

# Select variables for factor analysis
# These are hypothesized to reflect underlying latent constructs
fa_vars <- c("Adult_Mortality", "infant_deaths", "Alcohol", "Hepatitis_B", "BMI",
             "under_five_deaths", "Polio", "Total_expenditure", "Diphtheria", "HIV_AIDS",
             "thinness_10_19_years", "thinness_5_9_years", "Income_composition",
             "Schooling", "GDP_log")

X_fa <- df %>% select(all_of(fa_vars)) %>% drop_na()
df_complete <- df %>% select(Country, Status, all_of(fa_vars)) %>% drop_na()

cat("Variables:", length(fa_vars), ", Obs:", nrow(X_fa), "\n\n")

# Correlation matrix
R <- cor(X_fa)

# ============================================================================
# SECTION 1: SAMPLING ADEQUACY TESTS
# ============================================================================

# Kaiser-Meyer-Olkin (KMO) Measure
# - Proportion of variance that might be common variance
# - Interpretation: > 0.9 marvelous, > 0.8 meritorious, > 0.7 middling,
#                   > 0.6 mediocre, > 0.5 miserable, < 0.5 unacceptable
kmo <- KMO(R)
cat("=== KMO ===\n")
cat("Overall MSA:", round(kmo$MSA, 3), "\n")

# Bartlett's Test of Sphericity
# - Tests H₀: correlation matrix = identity (no correlations)
# - If significant (p < 0.05): correlations exist, FA is appropriate
bartlett <- cortest.bartlett(R, n = nrow(X_fa))
cat("\n=== Bartlett's Test ===\n")
cat("Chi-sq:", round(bartlett$chisq, 2), ", df:", bartlett$df,
    ", p:", format(bartlett$p.value, scientific = TRUE), "\n")

# ============================================================================
# SECTION 2: FACTOR NUMBER DETERMINATION
# ============================================================================

# Kaiser criterion (reference only)
eigenvalues <- eigen(R)$values
n_kaiser <- sum(eigenvalues > 1)
cat("\nKaiser criterion:", n_kaiser, "factors\n")

# Parallel Analysis (recommended method)
# Compares actual eigenvalues to eigenvalues from random data
# Retain factors where actual > random (simulated percentile)
cat("\n=== Parallel Analysis ===\n")
png("figures/fa_parallel_analysis.png", width = 900, height = 600, res = 120)
pa <- fa.parallel(X_fa, fa = "fa", n.iter = 100, main = "Parallel Analysis")
dev.off()
cat("Parallel analysis suggests:", pa$nfact, "factors\n")

# Use parallel analysis result
n_factors <- pa$nfact
cat("Using", n_factors, "factors (from parallel analysis)\n")

# ============================================================================
# SECTION 3: EXTRACTION METHOD SELECTION
# ============================================================================
# Two common extraction methods:
#   - Maximum Likelihood (ML): Assumes MVN, provides fit indices
#   - Principal Axis (PA): Robust, no distributional assumptions
#
# If MVN is violated, use PA (robust)

cat("\n=== MVN CHECK FOR EXTRACTION METHOD ===\n")
mvn_result <- MVN::mvn(X_fa, mvn_test = "mardia")
# Check MVN column - new MVN package uses "✓ Normal" instead of "YES"
mvn_pass <- grepl("Normal", mvn_result$multivariate_normality$MVN[1])

if (mvn_pass) {
  cat("MVN assumption met - using Maximum Likelihood (ML) extraction\n")
  fm_method <- "ml"
} else {
  cat("MVN assumption violated - using Principal Axis (PA) extraction (robust)\n")
  fm_method <- "pa"
}

# ============================================================================
# SECTION 4: FACTOR ANALYSIS WITH VARIMAX ROTATION
# ============================================================================
# Varimax: Orthogonal rotation that maximizes variance of squared loadings
# within each factor. Makes factors easier to interpret by pushing loadings
# toward 0 or ±1.

cat("\n=== FA WITH VARIMAX ===\n")
fa_varimax <- fa(X_fa, nfactors = n_factors, rotate = "varimax",
                 fm = fm_method, scores = "regression")
print(fa_varimax, cut = 0.3, sort = TRUE)  # Show loadings > 0.3

# ============================================================================
# SECTION 5: FACTOR ANALYSIS WITH PROMAX ROTATION
# ============================================================================
# Promax: Oblique rotation that allows factors to correlate
# More realistic when factors represent related constructs
# Factor correlations shown in fa_promax$Phi

cat("\n=== FA WITH PROMAX ===\n")
fa_promax <- fa(X_fa, nfactors = n_factors, rotate = "promax",
                fm = fm_method, scores = "regression")
print(fa_promax, cut = 0.3, sort = TRUE)

# ============================================================================
# SECTION 6: LOADINGS INTERPRETATION
# ============================================================================
# Loadings represent the correlation between each variable and each factor
# Interpretation guidelines:
#   |loading| >= 0.7: Excellent
#   |loading| >= 0.5: Good
#   |loading| >= 0.4: Acceptable minimum
#   |loading| < 0.3: Usually ignored

loadings_mat <- fa_varimax$loadings
class(loadings_mat) <- "matrix"
loadings_df <- as.data.frame(loadings_mat)
colnames(loadings_df) <- paste0("F", 1:n_factors)
loadings_df$Variable <- rownames(loadings_df)

# Communality: proportion of variance in each variable explained by all factors
# h² = sum of squared loadings for each variable
loadings_df$Communality <- fa_varimax$communality

cat("\n=== LOADINGS ===\n")
loadings_print <- loadings_df[, c("Variable", paste0("F", 1:n_factors), "Communality")]
loadings_print[, -1] <- round(loadings_print[, -1], 3)
print(loadings_print)
write.csv(loadings_df, "figures/fa_loadings.csv", row.names = FALSE)

# Check for low communalities
# Variables with h² < 0.4 may not belong to the factor structure
cat("\n=== COMMUNALITY CHECK ===\n")
low_comm_threshold <- 0.4
low_comm_vars <- loadings_df$Variable[loadings_df$Communality < low_comm_threshold]
if (length(low_comm_vars) > 0) {
  cat("WARNING: Variables with low communality (< 0.4):\n")
  for (v in low_comm_vars) {
    comm_val <- loadings_df$Communality[loadings_df$Variable == v]
    cat("  ", v, ": ", round(comm_val, 3), "\n", sep = "")
  }
  cat("Consider removing these variables or using more factors\n")
} else {
  cat("All variables have adequate communality (>= 0.4)\n")
}

# Loadings heatmap
loadings_long <- loadings_df %>%
  select(Variable, paste0("F", 1:n_factors)) %>%
  pivot_longer(cols = starts_with("F"), names_to = "Factor", values_to = "Loading")

p_loadings <- ggplot(loadings_long, aes(x = Factor, y = Variable, fill = Loading)) +
  geom_tile() +
  geom_text(aes(label = ifelse(abs(Loading) >= 0.3, round(Loading, 2), "")), size = 3) +
  scale_fill_gradient2(low = "#E94F37", mid = "white", high = "#2E86AB", limits = c(-1, 1)) +
  labs(title = "Factor Loadings (Varimax)") + theme_minimal()
ggsave("figures/fa_loadings_heatmap.png", p_loadings, width = 10, height = 10, dpi = 150)

# ============================================================================
# SECTION 7: FACTOR INTERPRETATION
# ============================================================================
# Name factors based on variables with high loadings
# This is subjective but grounded in domain knowledge

cat("\n=== FACTOR INTERPRETATION ===\n")
cat("F1: Economic Development (GDP, Schooling, Income)\n")
cat("F2: Healthcare Access (Diphtheria, Polio, Hep B)\n")
cat("F3: Mortality Burden (Adult/infant/under-5 deaths)\n")
cat("F4: Nutritional Status (BMI, thinness, Alcohol)\n")

# ============================================================================
# SECTION 8: FACTOR CORRELATIONS (FROM PROMAX)
# ============================================================================
# Promax allows factors to correlate
# Phi matrix shows these correlations
# High correlations suggest factors may not be distinct

if (!is.null(fa_promax$Phi)) {
  cat("\n=== FACTOR CORRELATIONS (Promax) ===\n")
  print(round(fa_promax$Phi, 3))

  png("figures/fa_factor_correlations.png", width = 600, height = 500, res = 120)
  corrplot(fa_promax$Phi, method = "color", type = "lower", addCoef.col = "black",
           tl.col = "black", title = "Factor Correlations", mar = c(0, 0, 2, 0))
  dev.off()
}

# ============================================================================
# SECTION 9: FACTOR SCORES
# ============================================================================
# Factor scores: estimated values for each observation on each factor
# Can be used as variables in subsequent analyses

factor_scores <- as.data.frame(fa_varimax$scores)
colnames(factor_scores) <- paste0("F", 1:n_factors)
factor_scores$Country <- df_complete$Country
factor_scores$Status <- df_complete$Status

cat("\n=== FACTOR SCORE MEANS BY STATUS ===\n")
factor_means <- factor_scores %>%
  group_by(Status) %>%
  summarise(across(starts_with("F"), ~round(mean(.), 3)))
print(factor_means)

write.csv(factor_scores, "figures/fa_factor_scores.csv", row.names = FALSE)

# Factor score plots
p_f1f2 <- ggplot(factor_scores, aes(x = F1, y = F2, color = Status)) +
  geom_point(alpha = 0.7, size = 2.5) + stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = "F1 vs F2", x = "F1 (Economic Development)", y = "F2 (Healthcare)") +
  theme_minimal()

p_f1f3 <- ggplot(factor_scores, aes(x = F1, y = F3, color = Status)) +
  geom_point(alpha = 0.7, size = 2.5) + stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = "F1 vs F3", x = "F1 (Economic Development)", y = "F3 (Mortality)") +
  theme_minimal()

ggsave("figures/fa_factor_scores_plot.png",
       grid.arrange(p_f1f2, p_f1f3, ncol = 2), width = 14, height = 6, dpi = 150)

# ============================================================================
# SECTION 10: MODEL FIT
# ============================================================================
# Model fit indices assess how well the factor model reproduces correlations
#
# TLI (Tucker-Lewis Index): > 0.95 excellent, > 0.90 acceptable
# RMSEA (Root Mean Square Error of Approximation):
#   < 0.05 good, < 0.08 acceptable, > 0.10 poor
# BIC (Bayesian Information Criterion): lower is better (for model comparison)

cat("\n=== MODEL FIT ===\n")
tli_val <- fa_varimax$TLI
rmsea_val <- fa_varimax$RMSEA[1]
bic_val <- fa_varimax$BIC

cat("TLI:", round(tli_val, 3))
if (tli_val >= 0.95) {
  cat(" (Excellent, >= 0.95)\n")
} else if (tli_val >= 0.90) {
  cat(" (Acceptable, >= 0.90)\n")
} else {
  cat(" (Poor, < 0.90 - consider more/fewer factors)\n")
}

cat("RMSEA:", round(rmsea_val, 3))
if (rmsea_val <= 0.05) {
  cat(" (Good, <= 0.05)\n")
} else if (rmsea_val <= 0.08) {
  cat(" (Acceptable, <= 0.08)\n")
} else {
  cat(" (Mediocre/Poor, > 0.08)\n")
}

cat("BIC:", round(bic_val, 2), "(lower is better)\n")

# ============================================================================
# SECTION 11: MODEL COMPARISON
# ============================================================================
# If model fit is poor, compare alternative specifications:
#   Model A: Remove Total_expenditure (low communality, h² < 0.4)
#   Model B: Increase number of factors (k + 1)
#   Model C: Original specification (baseline)

cat("\n=== FACTOR MODEL COMPARISON ===\n")
cat("Comparing alternative model specifications to improve fit\n\n")

# Model A: Remove low-communality variable (Total_expenditure if present)
# Check which variables have low communality
low_comm_vars_check <- names(fa_varimax$communality[fa_varimax$communality < 0.4])

if (length(low_comm_vars_check) > 0) {
  fa_vars_A <- setdiff(fa_vars, low_comm_vars_check[1])  # Remove worst
  cat("Model A: Removed", low_comm_vars_check[1], "(lowest communality)\n")
} else {
  fa_vars_A <- fa_vars[-length(fa_vars)]  # Remove last variable as fallback
  cat("Model A: Removed last variable as no low h² detected\n")
}

tryCatch({
  fa_A <- fa(X_fa[, fa_vars_A], nfactors = n_factors, rotate = "varimax", fm = fm_method)
  model_A_fit <- data.frame(
    TLI = fa_A$TLI,
    RMSEA = fa_A$RMSEA[1],
    BIC = fa_A$BIC
  )
}, error = function(e) {
  cat("Model A failed:", conditionMessage(e), "\n")
  model_A_fit <<- data.frame(TLI = NA, RMSEA = NA, BIC = NA)
})

# Model B: More factors (n_factors + 1)
n_factors_B <- min(n_factors + 1, floor(length(fa_vars) / 3))  # At most 1/3 of variables
cat("Model B: Using", n_factors_B, "factors (increased from", n_factors, ")\n")

tryCatch({
  fa_B <- fa(X_fa, nfactors = n_factors_B, rotate = "varimax", fm = fm_method)
  model_B_fit <- data.frame(
    TLI = fa_B$TLI,
    RMSEA = fa_B$RMSEA[1],
    BIC = fa_B$BIC
  )
}, error = function(e) {
  cat("Model B failed:", conditionMessage(e), "\n")
  model_B_fit <<- data.frame(TLI = NA, RMSEA = NA, BIC = NA)
})

# Model C: Original (baseline)
model_C_fit <- data.frame(
  TLI = fa_varimax$TLI,
  RMSEA = fa_varimax$RMSEA[1],
  BIC = fa_varimax$BIC
)

# Combine results
model_comparison <- data.frame(
  Model = c(paste0("A: Remove ", if(length(low_comm_vars_check) > 0) low_comm_vars_check[1] else "last var"),
            paste0("B: ", n_factors_B, " factors"),
            paste0("C: Original (", n_factors, " factors)")),
  Variables = c(length(fa_vars_A), length(fa_vars), length(fa_vars)),
  Factors = c(n_factors, n_factors_B, n_factors),
  TLI = round(c(model_A_fit$TLI, model_B_fit$TLI, model_C_fit$TLI), 3),
  RMSEA = round(c(model_A_fit$RMSEA, model_B_fit$RMSEA, model_C_fit$RMSEA), 3),
  BIC = round(c(model_A_fit$BIC, model_B_fit$BIC, model_C_fit$BIC), 1)
)

cat("\n=== MODEL COMPARISON TABLE ===\n")
print(model_comparison)
write.csv(model_comparison, "figures/fa_model_comparison.csv", row.names = FALSE)

# Determine best model
best_idx <- which.max(model_comparison$TLI)
cat("\nBest TLI:", model_comparison$Model[best_idx], "=", model_comparison$TLI[best_idx], "\n")

best_rmsea_idx <- which.min(model_comparison$RMSEA)
cat("Best RMSEA:", model_comparison$Model[best_rmsea_idx], "=", model_comparison$RMSEA[best_rmsea_idx], "\n")

best_bic_idx <- which.min(model_comparison$BIC)
cat("Best BIC:", model_comparison$Model[best_bic_idx], "=", model_comparison$BIC[best_bic_idx], "\n")

# Recommendation
cat("\n--- RECOMMENDATION ---\n")
if (model_comparison$TLI[best_idx] >= 0.90) {
  cat("Model", model_comparison$Model[best_idx], "achieves acceptable fit (TLI >= 0.90)\n")
} else if (model_comparison$TLI[best_idx] > model_C_fit$TLI) {
  cat("Model", model_comparison$Model[best_idx], "improves fit but still below threshold\n")
  cat("Consider: (1) removing more low-h² variables, (2) using more factors,\n")
  cat("or (3) accepting the FA as exploratory rather than confirmatory\n")
} else {
  cat("Original model has best or comparable fit\n")
  cat("Factor structure may not adequately represent the data\n")
}

cat("\n=== Factor Analysis Complete ===\n")
cat("Run 06_classification.R next.\n")

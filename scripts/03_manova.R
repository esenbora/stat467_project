# ============================================================================
# STAT 467 - MULTIVARIATE ANALYSIS OF VARIANCE (MANOVA)
# ============================================================================
# Purpose: Test whether Developed and Developing countries differ
#          across multiple health and socioeconomic outcomes simultaneously
#
# Input:   data_country_level.csv (country-level aggregated data)
#
# Output:  - Profile plot (figures/manova_profile_plot.png)
#          - Effect sizes (figures/manova_effect_sizes.png)
#          - Centroids plot (figures/manova_centroids.png)
#          - Test statistics (figures/manova_test_statistics.csv)
#          - Univariate ANOVAs (figures/univariate_anovas.csv)
#          - Discriminant coefficients (figures/discriminant_coefficients.csv)
#
# Methods:
#   - Four MANOVA test statistics (Pillai, Wilks, Hotelling-Lawley, Roy)
#   - Univariate follow-up ANOVAs with Bonferroni correction
#   - Discriminant function coefficients
#
# Key Concepts:
#   - MANOVA tests H₀: μ₁ = μ₂ (equal mean vectors across groups)
#   - Accounts for correlations between dependent variables
#   - Controls Type I error when testing multiple outcomes
#   - Different test statistics have different robustness properties
#
# Dependencies: tidyverse, MASS, car, biotools, MVN, ggpubr, rstatix, heplots
# ============================================================================

library(tidyverse)
library(MASS)
library(car)
# biotools for Box's M test (heplots has rgl dependency issues)
library(biotools)
library(MVN)
library(ggpubr)
library(rstatix)

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# Load data
df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n")
cat("Developed:", sum(df$Status == "Developed"),
    ", Developing:", sum(df$Status == "Developing"), "\n\n")

# Select dependent variables for MANOVA
# These are the health/socioeconomic outcomes we hypothesize differ by Status
# Note: HIV_AIDS excluded as it is constant (0.1) for all Developed countries,
#       causing singular covariance matrix in Box's M test
dv_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "Income_composition",
             "BMI", "GDP_log", "Diphtheria", "Polio")

# Prepare data with complete cases only
df_manova <- df %>% select(Country, Status, all_of(dv_vars)) %>% drop_na()
cat("Complete cases:", nrow(df_manova), "\n\n")

# ============================================================================
# SECTION 1: MANOVA ASSUMPTIONS
# ============================================================================
# MANOVA requires:
#   1. Independent observations (satisfied by country-level aggregation)
#   2. Multivariate normality within each group
#   3. Homogeneity of covariance matrices (Σ₁ = Σ₂)
#   4. Adequate sample size (n per group > p)

cat("=== ASSUMPTIONS ===\n")

# ----------------------------------------------------------------------------
# 1.1 Sample Size Check
# ----------------------------------------------------------------------------
# Rule of thumb: n per group should exceed number of dependent variables
n_per_group <- table(df_manova$Status)
print(n_per_group)
cat("n_min > p:", min(n_per_group), ">", length(dv_vars), "\n\n")

# ----------------------------------------------------------------------------
# 1.2 Multivariate Outlier Detection
# ----------------------------------------------------------------------------
# Outliers can severely distort MANOVA results
# Use Mahalanobis distance with chi-squared cutoff
cat("--- Outliers (Mahalanobis) ---\n")
X_outlier <- as.matrix(df_manova[, dv_vars])
center <- colMeans(X_outlier)
cov_mat <- cov(X_outlier)
mahal_dist <- mahalanobis(X_outlier, center, cov_mat)

# Chi-squared cutoff at alpha = 0.025 (two-tailed)
cutoff <- qchisq(0.975, df = length(dv_vars))
outliers <- which(mahal_dist > cutoff)
cat("Chi-sq cutoff (df=", length(dv_vars), ", alpha=0.025):", round(cutoff, 2), "\n")
cat("Potential outliers:", length(outliers), "countries\n")
if (length(outliers) > 0 && length(outliers) <= 10) {
  cat("Countries:", paste(df_manova$Country[outliers], collapse = ", "), "\n")
} else if (length(outliers) > 10) {
  cat("First 10:", paste(df_manova$Country[outliers[1:10]], collapse = ", "), "...\n")
}
cat("\n")

# ----------------------------------------------------------------------------
# 1.3 Multivariate Normality by Group
# ----------------------------------------------------------------------------
# Using Royston test (extension of Shapiro-Wilk for multivariate case)
# Appropriate for n < 5000
cat("--- Multivariate Normality by Group ---\n")
cat("Test selection: Royston (extension of Shapiro-Wilk, appropriate for n < 5000)\n")

mvn_by_group <- data.frame(Group = character(), Test = character(),
                            Statistic = numeric(), p_value = character(),
                            MVN = character(), stringsAsFactors = FALSE)

for (grp in levels(df_manova$Status)) {
  cat("\n", grp, " (n = ", sum(df_manova$Status == grp), "):\n", sep = "")
  grp_data <- df_manova %>% dplyr::filter(Status == grp) %>%
    dplyr::select(all_of(dv_vars))

  tryCatch({
    # Use Royston test
    mvn_royston <- MVN::mvn(grp_data, mvn_test = "royston")
    print(mvn_royston$multivariate_normality)

    mvn_by_group <- rbind(mvn_by_group, data.frame(
      Group = grp,
      Test = "Royston",
      Statistic = as.numeric(mvn_royston$multivariate_normality$Statistic),
      p_value = mvn_royston$multivariate_normality$p.value,
      MVN = mvn_royston$multivariate_normality$MVN
    ))
  }, error = function(e) {
    cat("Royston test failed:", conditionMessage(e), "\n")
    cat("Trying Mardia test...\n")
    tryCatch({
      mvn_mardia <- MVN::mvn(grp_data, mvn_test = "mardia")
      print(mvn_mardia$multivariate_normality)
    }, error = function(e2) {
      cat("MVN test failed:", conditionMessage(e2), "\n")
    })
  })
}

# ----------------------------------------------------------------------------
# 1.4 Univariate Normality by Group
# ----------------------------------------------------------------------------
# Check each variable within each group using Shapiro-Wilk
# Filter out constant variables (e.g., HIV_AIDS may be constant in Developed)
cat("\n--- Univariate Normality by Group (Shapiro-Wilk) ---\n")
df_for_sw <- df_manova %>%
  dplyr::select(Status, all_of(dv_vars))
univar_by_group <- df_for_sw %>%
  pivot_longer(cols = -Status, names_to = "Variable", values_to = "Value") %>%
  group_by(Status, Variable) %>%
  dplyr::filter(sd(Value) > 0) %>%  # Remove constant variables
  shapiro_test(Value)
univar_by_group$Normal <- ifelse(univar_by_group$p > 0.05, "Yes", "No")

# Count violations
violations <- univar_by_group %>%
  dplyr::filter(p < 0.05) %>%
  dplyr::select(Status, Variable, p)
if (nrow(violations) > 0) {
  cat("Variables violating univariate normality (p < 0.05):\n")
  print(violations)
} else {
  cat("All variables satisfy univariate normality in both groups\n")
}
cat("\n")

# ----------------------------------------------------------------------------
# 1.5 Box's M Test for Homogeneity of Covariance Matrices
# ----------------------------------------------------------------------------
# Tests H₀: Σ₁ = Σ₂ (equal covariance matrices)
# Box's M is VERY sensitive to MVN violations
# If rejected: use Pillai's Trace (most robust)
cat("--- Box's M Test (Covariance Homogeneity) ---\n")
X_matrix <- as.matrix(df_manova[, dv_vars])

# Using biotools::boxM (heplots has rgl dependency issues)
box_m <- biotools::boxM(X_matrix, df_manova$Status)
print(box_m)

cat("\n--- Box's M Interpretation ---\n")
if (is.infinite(box_m$statistic) || is.na(box_m$p.value)) {
  cat("Box's M test returned Inf or NA - possible singularity in covariance matrix\n")
  cat("This can occur with near-collinear variables or small group sizes\n")
  cat("Recommendation: Use Pillai's Trace (most robust to violations)\n")
  use_pillai <- TRUE
} else if (box_m$p.value < 0.05) {
  cat("Box's M significant (p < 0.05)\n")
  cat("Covariance matrices are NOT equal across groups\n")
  cat("Recommendation: Use Pillai's Trace (most robust to violations)\n")
  use_pillai <- TRUE
} else {
  cat("Box's M not significant (p >= 0.05)\n")
  cat("Homogeneity of covariance matrices assumption is met\n")
  cat("All test statistics are appropriate\n")
  use_pillai <- FALSE
}
cat("\n")

# ============================================================================
# SECTION 2: DESCRIPTIVE STATISTICS
# ============================================================================
cat("=== GROUP MEANS ===\n")
means_table <- df_manova %>%
  group_by(Status) %>%
  summarise(across(all_of(dv_vars), ~round(mean(.), 2)))
print(means_table)

# ============================================================================
# SECTION 3: MANOVA TEST
# ============================================================================
# MANOVA tests H₀: μ_Developed = μ_Developing (equal mean vectors)
#
# Four test statistics (each with different properties):
#   1. Pillai's Trace: Most robust to assumption violations
#   2. Wilks' Lambda: Most common, likelihood ratio test
#   3. Hotelling-Lawley: Generalization of T², powerful when k=2
#   4. Roy's Largest Root: Most powerful, but least robust
#
# For two groups, all four should give similar results

cat("\n=== MANOVA RESULTS ===\n")
dv_formula <- as.formula(paste("cbind(", paste(dv_vars, collapse = ", "), ") ~ Status"))
manova_fit <- manova(dv_formula, data = df_manova)

# Pillai's Trace: V = Σ λᵢ/(1+λᵢ) where λᵢ are eigenvalues
# Range: 0 to 1, larger = more separation
cat("\n1. PILLAI'S TRACE:\n")
print(summary(manova_fit, test = "Pillai"))

# Wilks' Lambda: Λ = Π 1/(1+λᵢ) = |W|/|W+B|
# Range: 0 to 1, smaller = more separation
cat("\n2. WILKS' LAMBDA:\n")
print(summary(manova_fit, test = "Wilks"))

# Hotelling-Lawley Trace: H = Σ λᵢ
# Larger = more separation
cat("\n3. HOTELLING-LAWLEY:\n")
print(summary(manova_fit, test = "Hotelling-Lawley"))

# Roy's Largest Root: θ = λ_max/(1+λ_max)
# Based on first eigenvalue only - most powerful but least robust
cat("\n4. ROY'S LARGEST ROOT:\n")
print(summary(manova_fit, test = "Roy"))

# Create summary table
pillai <- summary(manova_fit, test = "Pillai")
wilks <- summary(manova_fit, test = "Wilks")
hotelling <- summary(manova_fit, test = "Hotelling-Lawley")
roy <- summary(manova_fit, test = "Roy")

test_summary <- data.frame(
  Test = c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
  Statistic = round(c(pillai$stats[1,2], wilks$stats[1,2],
                      hotelling$stats[1,2], roy$stats[1,2]), 4),
  F = round(c(pillai$stats[1,3], wilks$stats[1,3],
              hotelling$stats[1,3], roy$stats[1,3]), 2),
  p_value = format(c(pillai$stats[1,6], wilks$stats[1,6],
                     hotelling$stats[1,6], roy$stats[1,6]),
                   scientific = TRUE, digits = 3)
)
cat("\n=== TEST SUMMARY ===\n")
print(test_summary)
write.csv(test_summary, "figures/manova_test_statistics.csv", row.names = FALSE)

# ============================================================================
# SECTION 4: UNIVARIATE FOLLOW-UP ANOVAS
# ============================================================================
# After significant MANOVA, conduct univariate ANOVAs to identify
# which specific dependent variables show group differences
# Apply Bonferroni correction: α* = 0.05/p

cat("\n=== UNIVARIATE ANOVAs ===\n")
anova_summary <- summary.aov(manova_fit)

p <- length(dv_vars)
alpha_bonf <- 0.05 / p  # Bonferroni-corrected alpha

univar_results <- data.frame(Variable = dv_vars, F = NA, p = NA, eta_sq = NA, Sig = NA)
for (i in 1:p) {
  aov_res <- anova_summary[[i]]
  univar_results$F[i] <- round(aov_res["Status", "F value"], 2)
  univar_results$p[i] <- aov_res["Status", "Pr(>F)"]
  # Eta-squared: η² = SS_effect / SS_total (effect size)
  univar_results$eta_sq[i] <- round(aov_res["Status", "Sum Sq"] /
                                     sum(aov_res[, "Sum Sq"]), 3)
  univar_results$Sig[i] <- ifelse(aov_res["Status", "Pr(>F)"] < alpha_bonf, "*", "")
}
print(univar_results)
write.csv(univar_results, "figures/univariate_anovas.csv", row.names = FALSE)

# ============================================================================
# SECTION 5: VISUALIZATIONS
# ============================================================================

# ----------------------------------------------------------------------------
# 5.1 Profile Plot (Standardized Means by Group)
# ----------------------------------------------------------------------------
# Shows pattern of differences across variables on common scale
df_profile <- df_manova %>%
  mutate(across(all_of(dv_vars), ~scale(.)[,1])) %>%
  pivot_longer(cols = all_of(dv_vars), names_to = "Variable", values_to = "Value")

profile_summary <- df_profile %>%
  group_by(Status, Variable) %>%
  summarise(Mean = mean(Value), SE = sd(Value)/sqrt(n()), .groups = "drop")

p_profile <- ggplot(profile_summary,
                    aes(x = Variable, y = Mean, group = Status, color = Status)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  scale_color_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
  labs(title = "Mean Profile Plot by Status", y = "Standardized Mean") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/manova_profile_plot.png", p_profile, width = 12, height = 7, dpi = 150)

# ----------------------------------------------------------------------------
# 5.2 Effect Size Plot
# ----------------------------------------------------------------------------
# Eta-squared interpretation:
#   η² ≈ 0.01: Small effect
#   η² ≈ 0.06: Medium effect
#   η² ≈ 0.14: Large effect
p_effect <- ggplot(univar_results, aes(x = reorder(Variable, eta_sq), y = eta_sq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = c(0.01, 0.06, 0.14),
             linetype = "dashed", color = c("green", "orange", "red")) +
  coord_flip() +
  labs(title = "Effect Sizes (Eta-squared)", x = "", y = expression(eta^2)) +
  theme_minimal()
ggsave("figures/manova_effect_sizes.png", p_effect, width = 10, height = 6, dpi = 150)

# ----------------------------------------------------------------------------
# 5.3 Centroids Plot (2D Projection)
# ----------------------------------------------------------------------------
# Shows group separation in 2D using two key variables
centroids <- df_manova %>%
  group_by(Status) %>%
  summarise(Life_expectancy = mean(Life_expectancy), Schooling = mean(Schooling))

p_centroids <- ggplot(df_manova,
                      aes(x = Life_expectancy, y = Schooling, color = Status)) +
  geom_point(alpha = 0.5) +
  geom_point(data = centroids, size = 5, shape = 18) +  # Diamond for centroids
  stat_ellipse(level = 0.95, linewidth = 1.2) +
  scale_color_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
  labs(title = "Group Centroids with 95% Ellipses") +
  theme_minimal()
ggsave("figures/manova_centroids.png", p_centroids, width = 10, height = 8, dpi = 150)

# ============================================================================
# SECTION 6: DISCRIMINANT COEFFICIENTS
# ============================================================================
# The discriminant function: D = a'x maximizes group separation
# Coefficients a = Sp⁻¹(x̄₁ - x̄₂) indicate variable contributions
#
# Raw coefficients: for original scale variables
# Standardized: a_std = a × √(diag(Sp)) - for comparing across variables

n1 <- sum(df_manova$Status == "Developed")
n2 <- sum(df_manova$Status == "Developing")

# Get group data
grp1_data <- df_manova %>% filter(Status == "Developed") %>% select(all_of(dv_vars))
grp2_data <- df_manova %>% filter(Status == "Developing") %>% select(all_of(dv_vars))

# Calculate covariance matrices
S1 <- cov(grp1_data)
S2 <- cov(grp2_data)

# Pooled covariance matrix: Sp = ((n1-1)S1 + (n2-1)S2) / (n1+n2-2)
Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)

# Mean difference from raw data (not rounded)
x_bar1 <- colMeans(grp1_data)
x_bar2 <- colMeans(grp2_data)
mean_diff <- x_bar1 - x_bar2

# Raw discriminant coefficients: a = Sp⁻¹ × (x̄₁ - x̄₂)
disc_coef <- solve(Sp) %*% mean_diff

# Standardized coefficients: a_std = a × √(diag(Sp))
# This gives coefficients for standardized (z-score) variables
disc_std <- disc_coef * sqrt(diag(Sp))

cat("\n=== DISCRIMINANT COEFFICIENTS ===\n")
cat("Raw coefficients: for original scale variables\n")
cat("Standardized: for comparison across variables with different scales\n\n")
disc_df <- data.frame(
  Variable = dv_vars,
  Raw = round(as.numeric(disc_coef), 6),
  Standardized = round(as.numeric(disc_std), 4)
)
disc_df <- disc_df %>% arrange(desc(abs(Standardized)))
print(disc_df)
write.csv(disc_df, "figures/discriminant_coefficients.csv", row.names = FALSE)

# ============================================================================
# SECTION 7: PERMUTATION MANOVA (ROBUST ALTERNATIVE)
# ============================================================================
# When MVN or covariance homogeneity assumptions are violated,
# permutation-based MANOVA provides robust inference
# adonis2 performs distance-based multivariate analysis of variance

cat("\n=== PERMUTATION MANOVA (ROBUST ALTERNATIVE) ===\n")
cat("Permutation tests do not require MVN or equal covariances\n\n")

library(vegan)

# Standardize variables for distance-based analysis
X_scaled <- scale(X_matrix)

# Permutation MANOVA using Euclidean distance
set.seed(123)
perm_manova <- adonis2(X_scaled ~ Status, data = df_manova,
                       method = "euclidean", permutations = 999)
cat("Permutation MANOVA (adonis2, 999 permutations):\n")
print(perm_manova)

# Compare with parametric results
pillai_p <- pillai$stats[1, 6]
perm_p <- perm_manova$`Pr(>F)`[1]

cat("\nComparison with parametric MANOVA (Pillai's Trace):\n")
cat("  Parametric p-value:   ", format(pillai_p, scientific = TRUE, digits = 3), "\n")
cat("  Permutation p-value:  ", format(perm_p, scientific = TRUE, digits = 3), "\n")

if (perm_p < 0.05 && pillai_p < 0.05) {
  cat("  Both tests significant - results are ROBUST\n")
  robust_status <- "Robust"
} else if (perm_p >= 0.05 && pillai_p >= 0.05) {
  cat("  Both tests non-significant - results are ROBUST\n")
  robust_status <- "Robust"
} else {
  cat("  Tests disagree - interpret with CAUTION\n")
  robust_status <- "Sensitive"
}

# ============================================================================
# SECTION 8: OUTLIER SENSITIVITY ANALYSIS
# ============================================================================
# Examine whether MANOVA conclusions change when outliers are removed

cat("\n=== OUTLIER SENSITIVITY ANALYSIS ===\n")

# Identify multivariate outliers using Mahalanobis distance
p_vars <- length(dv_vars)
cutoff <- qchisq(0.975, df = p_vars)

cat("Mahalanobis cutoff (chi-sq, df=", p_vars, ", alpha=0.025):", round(cutoff, 2), "\n")
cat("Outliers identified:", length(outliers), "countries\n")

if (length(outliers) > 0 && length(outliers) <= 15) {
  cat("Countries:", paste(df_manova$Country[outliers], collapse = ", "), "\n")
} else if (length(outliers) > 15) {
  cat("First 15:", paste(df_manova$Country[outliers[1:15]], collapse = ", "), "...\n")
}

# Re-run MANOVA without outliers
if (length(outliers) > 0 && length(outliers) < nrow(df_manova) * 0.2) {
  df_no_out <- df_manova[-outliers, ]

  cat("\nMANOVA without outliers (n =", nrow(df_no_out), "):\n")
  manova_no_out <- manova(dv_formula, data = df_no_out)
  pillai_no_out <- summary(manova_no_out, test = "Pillai")

  cat("  WITH outliers:    Pillai =", round(pillai$stats[1, 2], 4),
      ", p =", format(pillai_p, scientific = TRUE, digits = 3), "\n")
  cat("  WITHOUT outliers: Pillai =", round(pillai_no_out$stats[1, 2], 4),
      ", p =", format(pillai_no_out$stats[1, 6], scientific = TRUE, digits = 3), "\n")

  if (pillai_p < 0.05 && pillai_no_out$stats[1, 6] < 0.05) {
    cat("  Conclusion: Results are ROBUST to outlier removal\n")
    outlier_status <- "Robust"
  } else if (pillai_p >= 0.05 && pillai_no_out$stats[1, 6] >= 0.05) {
    cat("  Conclusion: Results are ROBUST (both non-significant)\n")
    outlier_status <- "Robust"
  } else {
    cat("  Conclusion: Results are SENSITIVE to outliers\n")
    outlier_status <- "Sensitive"
  }
} else if (length(outliers) == 0) {
  cat("No outliers detected - sensitivity analysis not needed\n")
  outlier_status <- "N/A"
} else {
  cat("Too many outliers (>20%) - sensitivity analysis not meaningful\n")
  outlier_status <- "N/A"
}

# ============================================================================
# SECTION 9: STATISTICAL VALIDITY SUMMARY
# ============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("STATISTICAL VALIDITY SUMMARY\n")
cat(strrep("=", 60), "\n")

# MVN status (check if any group failed)
mvn_status <- if (all(mvn_by_group$MVN == "YES" | grepl("Normal", mvn_by_group$MVN))) {
  "Met (both groups)"
} else {
  "Violated (use Pillai's Trace)"
}
cat("MVN Assumption:            ", mvn_status, "\n")

# Box's M status
if (is.infinite(box_m$statistic) || is.na(box_m$p.value)) {
  boxm_status <- "Indeterminate (singularity)"
} else if (box_m$p.value < 0.05) {
  boxm_status <- "Violated (use Pillai's Trace)"
} else {
  boxm_status <- "Met"
}
cat("Covariance Homogeneity:    ", boxm_status, "\n")
cat("Outlier Sensitivity:       ", outlier_status, "\n")
cat("Permutation Robustness:    ", robust_status, "\n")
cat("Recommended Test:           Pillai's Trace (most robust)\n")
cat(strrep("=", 60), "\n")

cat("\n=== MANOVA Complete ===\n")
cat("Run 04_pca_regression.R next.\n")

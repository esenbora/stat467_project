# ============================================================================
# STAT 467 - MEAN VECTOR INFERENCE
# ============================================================================
# Purpose: Perform hypothesis tests on multivariate mean vectors
#          Compare Developed vs Developing countries using Hotelling's T²
#
# Input:   data_country_level.csv (country-level aggregated data)
#
# Output:  - Q-Q plot (figures/mean_vector_qq_plot.png)
#          - Confidence ellipse (figures/confidence_ellipse.png)
#
# Methods:
#   1. One-sample Hotelling's T² test (compare to hypothesized mean)
#   2. Two-sample Hotelling's T² test (compare two groups)
#   3. Bonferroni simultaneous confidence intervals
#
# Key Concepts:
#   - Hotelling's T² is the multivariate extension of the t-test
#   - T² = n(x̄ - μ₀)'S⁻¹(x̄ - μ₀) follows scaled F-distribution
#   - Assumptions: MVN within groups, equal covariances (for two-sample)
#   - Bonferroni correction controls family-wise error rate
#
# Dependencies: tidyverse, MASS, Hotelling, MVN, ellipse, rstatix, heplots
# ============================================================================

library(tidyverse)
library(MASS)
library(Hotelling)
library(MVN)
library(ellipse)
library(rstatix)
# heplots has rgl dependency which may fail on some systems
# biotools provides alternative boxM implementation
library(biotools)

# Resolve namespace conflicts (MASS::select vs dplyr::select)
select <- dplyr::select
filter <- dplyr::filter

# Load country-level data
df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n")
cat("Developed:", sum(df$Status == "Developed"),
    ", Developing:", sum(df$Status == "Developing"), "\n\n")

# Select variables for mean vector inference
# These represent key health and socioeconomic indicators
# Note: HIV_AIDS excluded as it is constant (0.1) for all Developed countries,
#       causing singular covariance matrix in Box's M test
selected_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling",
                   "GDP_log", "BMI", "Diphtheria")

# Create data matrix
X <- df %>% select(all_of(selected_vars)) %>% as.matrix()
n <- nrow(X)  # Total sample size
p <- ncol(X)  # Number of variables

# Calculate sample statistics
x_bar <- colMeans(X)  # Sample mean vector (p x 1)
S <- cov(X)           # Sample covariance matrix (p x p)

cat("=== Sample Mean Vector ===\n")
print(round(x_bar, 3))

# ============================================================================
# SECTION 1: ASSUMPTION CHECKING FOR HOTELLING'S T²
# ============================================================================
# Hotelling's T² requires:
#   1. Random sample from multivariate normal population
#   2. For two-sample: equal covariance matrices (Σ₁ = Σ₂)
#
# Robustness: T² is somewhat robust to MVN violations when n is large

# ----------------------------------------------------------------------------
# 1.1 Multivariate Normality (Overall Sample)
# ----------------------------------------------------------------------------
cat("\n=== MULTIVARIATE NORMALITY CHECK ===\n")

# Royston test: recommended for n < 5000
mvn_royston <- MVN::mvn(as.data.frame(X), mvn_test = "royston")
cat("Royston's MVN Test (overall):\n")
print(mvn_royston$multivariate_normality)

# Mardia's test: provides separate tests for skewness and kurtosis
mvn_mardia <- MVN::mvn(as.data.frame(X), mvn_test = "mardia")
cat("\nMardia's MVN Test (overall):\n")
print(mvn_mardia$multivariate_normality)

# ----------------------------------------------------------------------------
# 1.2 Univariate Normality
# ----------------------------------------------------------------------------
# If univariate normality fails, MVN cannot hold (necessary condition)
cat("\n--- Univariate Normality (Shapiro-Wilk) ---\n")
univar_df <- as.data.frame(X) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  shapiro_test(Value)
univar_df$Normal <- ifelse(univar_df$p > 0.05, "Yes", "No")
print(univar_df)

# Q-Q plot for visual MVN assessment
png("figures/mean_vector_qq_plot.png", width = 800, height = 600, res = 120)
plot(mvn_mardia, diagnostic = "multivariate", type = "qq")
dev.off()

# ============================================================================
# SECTION 2: ONE-SAMPLE HOTELLING'S T² TEST
# ============================================================================
# Tests H₀: μ = μ₀ (population mean equals hypothesized value)
# vs. H₁: μ ≠ μ₀ (population mean differs)
#
# Test statistic: T² = n(x̄ - μ₀)'S⁻¹(x̄ - μ₀)
#
# Under H₀: ((n-p)/(p(n-1))) * T² ~ F(p, n-p)
#
# This transformation allows use of F-distribution for p-values

cat("\n=== ONE-SAMPLE HOTELLING'S T² TEST ===\n")

# Define hypothesized mean vector (based on domain knowledge)
# These represent "reasonable" benchmark values for global averages
mu0 <- c(Life_expectancy = 70,     # Years
         Adult_Mortality = 150,    # Deaths per 1000 adults
         Schooling = 12,           # Years of schooling
         GDP_log = 8.5,            # log(1 + GDP per capita)
         BMI = 25,                 # Body Mass Index
         Diphtheria = 85)          # % immunization coverage

cat("Hypothesized mean:\n")
print(mu0)

# Calculate T² statistic
diff <- x_bar - mu0                           # Deviation from hypothesis
T2 <- n * t(diff) %*% solve(S) %*% diff       # Hotelling's T²
F_stat <- ((n - p) / (p * (n - 1))) * T2      # Transform to F
p_value <- 1 - pf(F_stat, p, n - p)           # Calculate p-value

cat("\nT² =", round(T2, 4), ", F =", round(F_stat, 4),
    ", p =", format(p_value, scientific = TRUE), "\n")

if (p_value < 0.05) {
  cat("Reject H0: Mean vector differs from hypothesized values\n")
} else {
  cat("Fail to reject H0\n")
}

# ============================================================================
# SECTION 3: BONFERRONI SIMULTANEOUS CONFIDENCE INTERVALS
# ============================================================================
# After rejecting H₀ in Hotelling's T², we want to identify which
# individual means differ from their hypothesized values.
#
# Bonferroni correction: α_individual = α/p
# This ensures family-wise error rate ≤ α (here, α = 0.05)
#
# 100(1-α)% simultaneous CI for μᵢ:
#   x̄ᵢ ± t_{α/(2p), n-1} × √(sᵢᵢ/n)

cat("\n=== BONFERRONI CONFIDENCE INTERVALS ===\n")
alpha <- 0.05
t_crit <- qt(1 - alpha/(2*p), df = n - 1)  # Bonferroni-adjusted critical value

# Calculate CIs for each variable
ci_bonf <- data.frame(
  Variable = selected_vars,
  Mean = round(x_bar, 3),
  Lower = round(x_bar - t_crit * sqrt(diag(S)/n), 3),
  Upper = round(x_bar + t_crit * sqrt(diag(S)/n), 3),
  H0 = mu0
)
print(ci_bonf)

# ============================================================================
# SECTION 4: CONFIDENCE ELLIPSE VISUALIZATION
# ============================================================================
# For 2D visualization, project onto two key variables
# The ellipse represents the 95% confidence region for the mean vector
#
# Ellipse equation: n(x̄ - μ)'S⁻¹(x̄ - μ) ≤ c²
# where c² = (2(n-1)/(n-2)) × F_{0.95, 2, n-2}

png("figures/confidence_ellipse.png", width = 900, height = 700, res = 120)

# Select two variables for 2D plot
X_2d <- df %>% select(Life_expectancy, Schooling) %>% as.matrix()
x_bar_2d <- colMeans(X_2d)
S_2d <- cov(X_2d)

# Calculate critical value for 95% confidence ellipse
# c² = (p(n-1)/(n-p)) × F_{α, p, n-p} where p=2 for 2D
c2 <- (2 * (n - 1) / (n - 2)) * qf(0.95, 2, n - 2)

# Generate ellipse points using eigenvalue decomposition
# Ellipse axes are aligned with eigenvectors, lengths proportional to sqrt(eigenvalues)
theta <- seq(0, 2*pi, length.out = 100)
eigen_S <- eigen(S_2d)
axes <- sqrt(c2/n) * sqrt(eigen_S$values)
rotation <- eigen_S$vectors
ellipse_pts <- t(rotation %*% rbind(axes[1]*cos(theta), axes[2]*sin(theta))) +
               matrix(x_bar_2d, nrow = 100, ncol = 2, byrow = TRUE)

# Create plot
plot(X_2d, pch = 19, col = "gray50",
     xlab = "Life Expectancy", ylab = "Schooling",
     main = "95% Confidence Ellipse for Mean Vector")
lines(ellipse_pts, col = "steelblue", lwd = 3)
points(x_bar_2d[1], x_bar_2d[2], pch = 19, col = "red", cex = 2)
points(mu0["Life_expectancy"], mu0["Schooling"],
       pch = 4, col = "darkgreen", cex = 2, lwd = 3)
legend("bottomright",
       c("Countries", "Sample Mean", "Hypothesized Mean", "95% CR"),
       col = c("gray50", "red", "darkgreen", "steelblue"),
       pch = c(19, 19, 4, NA), lty = c(NA, NA, NA, 1), lwd = c(NA, NA, 3, 3))
dev.off()

# ============================================================================
# SECTION 5: TWO-SAMPLE HOTELLING'S T² TEST
# ============================================================================
# Tests H₀: μ₁ = μ₂ (Developed and Developing have same mean vector)
# vs. H₁: μ₁ ≠ μ₂
#
# This is the multivariate analog of the two-sample t-test
# Uses pooled covariance matrix when Σ₁ = Σ₂ (checked via Box's M)

cat("\n=== TWO-SAMPLE HOTELLING'S T² TEST ===\n")

# Separate data by group
X_dev <- df %>% dplyr::filter(Status == "Developed") %>%
  dplyr::select(all_of(selected_vars)) %>% as.matrix()
X_devp <- df %>% dplyr::filter(Status == "Developing") %>%
  dplyr::select(all_of(selected_vars)) %>% as.matrix()

n1 <- nrow(X_dev)   # Developed sample size
n2 <- nrow(X_devp)  # Developing sample size

x_bar1 <- colMeans(X_dev)   # Developed mean vector
x_bar2 <- colMeans(X_devp)  # Developing mean vector

cat("Sample sizes: n1 =", n1, ", n2 =", n2, "\n\n")

# ----------------------------------------------------------------------------
# 5.1 Assumption 1: Multivariate Normality by Group
# ----------------------------------------------------------------------------
cat("--- Assumption 1: MVN by Group ---\n")

# Developed group
cat("\nDeveloped (n =", n1, "):\n")
tryCatch({
  mvn_dev <- MVN::mvn(as.data.frame(X_dev), mvn_test = "mardia")
  print(mvn_dev$multivariate_normality)
}, error = function(e) {
  cat("MVN test failed (singular covariance, likely constant variable).\n")
  cat("Checking univariate normality instead:\n")
  for (v in selected_vars) {
    tryCatch({
      sw_test <- shapiro.test(X_dev[, v])
      cat("  ", v, ": p =", format(sw_test$p.value, digits = 3), "\n")
    }, error = function(e2) {
      cat("  ", v, ": constant or NA values\n")
    })
  }
})

# Developing group
cat("\nDeveloping (n =", n2, "):\n")
tryCatch({
  mvn_devp <- MVN::mvn(as.data.frame(X_devp), mvn_test = "royston")
  print(mvn_devp$multivariate_normality)
}, error = function(e) {
  cat("Royston test failed. Using Mardia:\n")
  tryCatch({
    mvn_devp <- MVN::mvn(as.data.frame(X_devp), mvn_test = "mardia")
    print(mvn_devp$multivariate_normality)
  }, error = function(e2) {
    cat("MVN test failed:", conditionMessage(e2), "\n")
  })
})

# Univariate normality by group using rstatix
# Filter out constant variables (e.g., HIV_AIDS may be constant in Developed)
cat("\n--- Univariate Normality by Group (Shapiro-Wilk) ---\n")
df_for_test <- df %>%
  dplyr::select(Status, all_of(selected_vars))
univar_by_group <- df_for_test %>%
  pivot_longer(cols = -Status, names_to = "Variable", values_to = "Value") %>%
  group_by(Status, Variable) %>%
  # Filter out groups with constant values
  dplyr::filter(sd(Value) > 0) %>%
  shapiro_test(Value)
univar_by_group$Normal <- ifelse(univar_by_group$p > 0.05, "Yes", "No")
print(univar_by_group)

# ----------------------------------------------------------------------------
# 5.2 Assumption 2: Homogeneity of Covariance Matrices (Box's M Test)
# ----------------------------------------------------------------------------
# Box's M tests H₀: Σ₁ = Σ₂ = ... = Σₖ
# If rejected: use separate covariances or robust methods
# If not rejected: pooled covariance is appropriate
cat("\n--- Assumption 2: Box's M Test (Covariance Homogeneity) ---\n")
X_combined <- rbind(X_dev, X_devp)
groups <- factor(c(rep("Developed", n1), rep("Developing", n2)))
box_m <- biotools::boxM(X_combined, groups)
print(box_m)

if (is.infinite(box_m$statistic) || is.na(box_m$p.value)) {
  cat("\nBox's M test returned Inf/NA - possible singularity issue\n")
  cat("Using pooled covariance anyway, but interpret with caution\n")
  use_pooled <- TRUE
} else if (box_m$p.value < 0.05) {
  cat("\nBox's M significant (p < 0.05): Covariance matrices are NOT equal\n")
  cat("WARNING: Pooled covariance assumption violated!\n")
  cat("Consider: Nel-Van der Merwe test or separate covariances\n")
  cat("Proceeding with pooled Sp but results should be interpreted cautiously\n")
  use_pooled <- FALSE
} else {
  cat("\nBox's M not significant (p >= 0.05): Covariance matrices are equal\n")
  cat("Pooled covariance matrix (Sp) is appropriate\n")
  use_pooled <- TRUE
}

# ----------------------------------------------------------------------------
# 5.3 Two-Sample Hotelling's T² Calculation
# ----------------------------------------------------------------------------
# Pooled covariance: Sp = ((n1-1)S1 + (n2-1)S2) / (n1+n2-2)
# T² = (n1n2/(n1+n2)) × (x̄1 - x̄2)'Sp⁻¹(x̄1 - x̄2)
# Under H₀: ((n1+n2-p-1)/(p(n1+n2-2))) × T² ~ F(p, n1+n2-p-1)

cat("\nMean difference (Developed - Developing):\n")
print(round(x_bar1 - x_bar2, 3))

# Pooled covariance matrix
Sp <- ((n1-1)*cov(X_dev) + (n2-1)*cov(X_devp)) / (n1 + n2 - 2)
diff_means <- x_bar1 - x_bar2
T2_two <- (n1*n2/(n1+n2)) * t(diff_means) %*% solve(Sp) %*% diff_means
F_two <- ((n1 + n2 - p - 1) / (p * (n1 + n2 - 2))) * T2_two
p_two <- 1 - pf(F_two, p, n1 + n2 - p - 1)

cat("\nT² =", round(T2_two, 4), ", F =", round(F_two, 4),
    ", p =", format(p_two, scientific = TRUE), "\n")

if (p_two < 0.05) {
  cat("Reject H0: Developed and Developing mean vectors differ significantly\n")
}

# Verification using Hotelling package
cat("\n=== Hotelling Package Verification ===\n")
print(hotelling.test(X_dev, X_devp))

# ============================================================================
# SECTION 6: UNIVARIATE FOLLOW-UP TESTS
# ============================================================================
# After significant Hotelling's T², identify which variables
# contribute to the difference using univariate t-tests
# Apply Bonferroni correction: α* = α/p

cat("\n=== UNIVARIATE t-TESTS ===\n")
for (i in 1:p) {
  tt <- t.test(X_dev[,i], X_devp[,i], var.equal = TRUE)
  sig <- ifelse(tt$p.value < alpha/p, "*", "")  # Bonferroni-significant
  cat(selected_vars[i], ": t =", round(tt$statistic, 2),
      ", p =", format(tt$p.value, digits = 3), sig, "\n")
}

# ============================================================================
# SECTION 7: PERMUTATION TEST (ROBUST ALTERNATIVE)
# ============================================================================
# When MVN is violated, permutation-based tests provide robust inference
# This tests multivariate group differences without distributional assumptions

cat("\n=== PERMUTATION TEST (ROBUST ALTERNATIVE) ===\n")
cat("When MVN assumption is violated, permutation tests provide robust inference\n\n")

# Use vegan::adonis2 for permutation-based MANOVA (distance-based)
library(vegan)

# Create a distance matrix from the standardized data
X_scaled <- scale(X)
perm_result <- adonis2(X_scaled ~ Status, data = df, method = "euclidean",
                       permutations = 999)
cat("Permutation MANOVA (adonis2, 999 permutations):\n")
print(perm_result)

cat("\nComparison with parametric Hotelling's T²:\n")
cat("  Parametric p-value:   ", format(p_two, scientific = TRUE, digits = 3), "\n")
cat("  Permutation p-value:  ", format(perm_result$`Pr(>F)`[1], scientific = TRUE, digits = 3), "\n")

if (perm_result$`Pr(>F)`[1] < 0.05 && p_two < 0.05) {
  cat("  Both tests significant - results are ROBUST\n")
} else if (perm_result$`Pr(>F)`[1] >= 0.05 && p_two >= 0.05) {
  cat("  Both tests non-significant - results are ROBUST\n")
} else {
  cat("  Tests disagree - interpret with CAUTION\n")
}

# ============================================================================
# SECTION 8: OUTLIER SENSITIVITY ANALYSIS
# ============================================================================
# Examine whether conclusions change when multivariate outliers are removed
# Outliers identified via Mahalanobis distance with chi-squared cutoff

cat("\n=== OUTLIER SENSITIVITY ANALYSIS ===\n")

# Identify multivariate outliers
cutoff <- qchisq(0.975, df = p)
mahal_dist <- mahalanobis(X, colMeans(X), cov(X))
outlier_idx <- which(mahal_dist > cutoff)
outlier_countries <- df$Country[outlier_idx]

cat("Mahalanobis cutoff (chi-sq, df=", p, ", alpha=0.025):", round(cutoff, 2), "\n")
cat("Outliers identified:", length(outlier_idx), "countries\n")
if (length(outlier_idx) > 0 && length(outlier_idx) <= 15) {
  cat("Countries:", paste(outlier_countries, collapse = ", "), "\n")
} else if (length(outlier_idx) > 15) {
  cat("First 15:", paste(outlier_countries[1:15], collapse = ", "), "...\n")
}

# Re-run two-sample test WITHOUT outliers
if (length(outlier_idx) > 0) {
  df_no_out <- df[-outlier_idx, ]
  X_no_out_dev <- df_no_out %>% dplyr::filter(Status == "Developed") %>%
    dplyr::select(all_of(selected_vars)) %>% as.matrix()
  X_no_out_devp <- df_no_out %>% dplyr::filter(Status == "Developing") %>%
    dplyr::select(all_of(selected_vars)) %>% as.matrix()

  n1_no <- nrow(X_no_out_dev)
  n2_no <- nrow(X_no_out_devp)

  # Pooled covariance and T² without outliers
  Sp_no <- ((n1_no-1)*cov(X_no_out_dev) + (n2_no-1)*cov(X_no_out_devp)) / (n1_no + n2_no - 2)
  diff_no <- colMeans(X_no_out_dev) - colMeans(X_no_out_devp)
  T2_no <- (n1_no*n2_no/(n1_no+n2_no)) * t(diff_no) %*% solve(Sp_no) %*% diff_no
  F_no <- ((n1_no + n2_no - p - 1) / (p * (n1_no + n2_no - 2))) * T2_no
  p_no <- 1 - pf(F_no, p, n1_no + n2_no - p - 1)

  cat("\nTwo-sample Hotelling's T² comparison:\n")
  cat("  WITH outliers:    T² =", round(T2_two, 2), ", p =", format(p_two, scientific = TRUE, digits = 3), "\n")
  cat("  WITHOUT outliers: T² =", round(T2_no, 2), ", p =", format(p_no, scientific = TRUE, digits = 3), "\n")

  if (p_two < 0.05 && p_no < 0.05) {
    cat("  Conclusion: Results are ROBUST to outlier removal\n")
    outlier_status <- "Robust"
  } else if (p_two >= 0.05 && p_no >= 0.05) {
    cat("  Conclusion: Results are ROBUST (both non-significant)\n")
    outlier_status <- "Robust"
  } else {
    cat("  Conclusion: Results are SENSITIVE to outliers - interpret with caution\n")
    outlier_status <- "Sensitive"
  }
} else {
  cat("No outliers detected - sensitivity analysis not needed\n")
  outlier_status <- "N/A"
}

# ============================================================================
# SECTION 9: BONFERRONI CONFIDENCE INTERVALS BY STATUS
# ============================================================================
# Simultaneous confidence intervals for each group (Developed vs Developing)
# Adjusted for multiple comparisons across groups and variables

cat("\n=== BONFERRONI CONFIDENCE INTERVALS BY STATUS ===\n")

k_status <- 2  # 2 groups
p_status <- length(selected_vars)  # 6 variables
alpha_bonf_status <- alpha / (k_status * p_status)
t_crit_status <- qt(1 - alpha_bonf_status/2, df = n - k_status)

cat("Number of groups (k):", k_status, "\n")
cat("Number of variables (p):", p_status, "\n")
cat("Bonferroni alpha:", round(alpha_bonf_status, 5), "\n")
cat("Critical t-value:", round(t_crit_status, 3), "\n\n")

# Calculate CI for each Status and variable
ci_status_long <- df %>%
  group_by(Status) %>%
  dplyr::summarise(
    n = n(),
    across(all_of(selected_vars), list(mean = mean, sd = sd))
  ) %>%
  pivot_longer(
    cols = -c(Status, n),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(mean|sd)"
  ) %>%
  mutate(
    se = sd / sqrt(n),
    ci_lower = mean - t_crit_status * se,
    ci_upper = mean + t_crit_status * se
  )

# Plot 1: Bonferroni CI for Life Expectancy by Status
p_ci_status_life <- ggplot(ci_status_long %>% dplyr::filter(Variable == "Life_expectancy"),
                           aes(x = Status, y = mean, color = Status)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 1.2) +
  scale_color_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
  labs(title = "Bonferroni 95% CI: Life Expectancy by Status",
       subtitle = paste0("Adjusted for ", k_status * p_status, " comparisons (α = ", round(alpha_bonf_status, 5), ")"),
       x = "", y = "Life Expectancy (years)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

ggsave("figures/bonferroni_ci_status_life_expectancy.png", p_ci_status_life, width = 8, height = 6, dpi = 150)
cat("Saved: figures/bonferroni_ci_status_life_expectancy.png\n")

# Plot 2: Bonferroni CI for all variables by Status (faceted)
p_ci_status_all <- ggplot(ci_status_long, aes(x = Status, y = mean, color = Status)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 1) +
  scale_color_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
  facet_wrap(~Variable, scales = "free_y", ncol = 3) +
  labs(title = "Bonferroni 95% CI by Status (Developed vs Developing)",
       subtitle = paste0("α = ", round(alpha_bonf_status, 5), " (Bonferroni-adjusted)"),
       x = "", y = "Mean Value") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0))

ggsave("figures/bonferroni_ci_status_all_variables.png", p_ci_status_all, width = 12, height = 8, dpi = 150)
cat("Saved: figures/bonferroni_ci_status_all_variables.png\n")

# Save CI data
write.csv(ci_status_long, "figures/bonferroni_ci_status_data.csv", row.names = FALSE)

# Print CI table
cat("\nBonferroni CI Summary:\n")
ci_summary <- ci_status_long %>%
  mutate(CI = paste0("[", round(ci_lower, 2), ", ", round(ci_upper, 2), "]")) %>%
  dplyr::select(Status, Variable, mean, CI) %>%
  pivot_wider(names_from = Status, values_from = c(mean, CI))
print(ci_summary)

# ============================================================================
# SECTION 10: STATISTICAL VALIDITY SUMMARY
# ============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("STATISTICAL VALIDITY SUMMARY\n")
cat(strrep("=", 60), "\n")

# Determine MVN status
mvn_overall <- ifelse(mvn_royston$multivariate_normality$MVN == "YES" ||
                      grepl("Normal", mvn_royston$multivariate_normality$MVN),
                      "Met", "Violated")
cat("MVN Assumption (Overall):  ", mvn_overall, "\n")

# Box's M status
if (is.infinite(box_m$statistic) || is.na(box_m$p.value)) {
  boxm_status <- "Indeterminate (singularity)"
} else if (box_m$p.value < 0.05) {
  boxm_status <- "Violated (p < 0.05)"
} else {
  boxm_status <- "Met (p >= 0.05)"
}
cat("Covariance Homogeneity:    ", boxm_status, "\n")
cat("Outlier Sensitivity:       ", outlier_status, "\n")
cat("Robust Alternative:         Permutation MANOVA (vegan::adonis2)\n")
cat(strrep("=", 60), "\n")

cat("\n=== Mean Vector Inference Complete ===\n")
cat("Run 03_manova.R next.\n")

# ============================================================================
# STAT 467 - MULTIVARIATE ANALYSIS OF VARIANCE (MANOVA) BY CONTINENT
# ============================================================================
# Purpose: Test whether the 6 Continents differ across multiple health and
#          socioeconomic outcomes simultaneously
#
# Input:   data_country_level.csv (country-level aggregated data with Continent)
#
# Output:  - Profile plot (figures/manova_continent_profiles.png)
#          - Boxplots (figures/manova_continent_boxplots.png)
#          - Pairwise heatmap (figures/manova_pairwise_heatmap.png)
#          - Effect sizes (figures/manova_effect_sizes.png)
#          - Test statistics (figures/manova_test_statistics.csv)
#          - Univariate ANOVAs (figures/univariate_anovas.csv)
#          - Pairwise comparisons (figures/manova_pairwise_comparisons.csv)
#
# Methods:
#   - Four MANOVA test statistics (Pillai, Wilks, Hotelling-Lawley, Roy)
#   - Univariate follow-up ANOVAs with Bonferroni correction
#   - Post-hoc pairwise comparisons (Bonferroni-corrected)
#
# Key Concepts:
#   - MANOVA tests H0: mu_1 = mu_2 = ... = mu_6 (equal mean vectors)
#   - Accounts for correlations between dependent variables
#   - Controls Type I error when testing multiple outcomes
#   - Post-hoc tests identify WHICH groups differ
#
# Dependencies: tidyverse, MASS, car, biotools, MVN, ggpubr, rstatix
# ============================================================================

library(tidyverse)
library(MASS)
library(car)
library(biotools)
library(MVN)
library(ggpubr)
library(rstatix)

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# Load data
df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Continent <- as.factor(df$Continent)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n\n")
cat("=== CONTINENT DISTRIBUTION ===\n")
print(table(df$Continent))
cat("\n")

# Select dependent variables for MANOVA
# These are the health/socioeconomic outcomes we hypothesize differ by Continent
dv_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "Income_composition",
             "BMI", "GDP_log", "Diphtheria", "Polio")

# Prepare data with complete cases only
df_manova <- df %>% select(Country, Continent, Status, all_of(dv_vars)) %>% drop_na()
cat("Complete cases:", nrow(df_manova), "\n\n")

cat("Sample sizes by Continent:\n")
print(table(df_manova$Continent))
cat("\n")

# ============================================================================
# SECTION 1: MANOVA ASSUMPTIONS
# ============================================================================
# MANOVA requires:
#   1. Independent observations (satisfied by country-level aggregation)
#   2. Multivariate normality within each group
#   3. Homogeneity of covariance matrices (all Sigma_g equal)
#   4. Adequate sample size (n per group > p)

cat("=== ASSUMPTIONS ===\n")

# ----------------------------------------------------------------------------
# 1.1 Sample Size Check
# ----------------------------------------------------------------------------
n_per_group <- table(df_manova$Continent)
cat("Sample sizes by Continent:\n")
print(n_per_group)
cat("\nMinimum n:", min(n_per_group), " | Number of DVs (p):", length(dv_vars), "\n")
cat("n_min > p:", min(n_per_group), ">", length(dv_vars), "=", min(n_per_group) > length(dv_vars), "\n\n")

# ----------------------------------------------------------------------------
# 1.2 Multivariate Outlier Detection
# ----------------------------------------------------------------------------
cat("--- Outliers (Mahalanobis) ---\n")
X_outlier <- as.matrix(df_manova[, dv_vars])
center <- colMeans(X_outlier)
cov_mat <- cov(X_outlier)
mahal_dist <- mahalanobis(X_outlier, center, cov_mat)

cutoff <- qchisq(0.975, df = length(dv_vars))
outliers <- which(mahal_dist > cutoff)
cat("Chi-sq cutoff (df=", length(dv_vars), ", alpha=0.025):", round(cutoff, 2), "\n")
cat("Potential outliers:", length(outliers), "countries\n")
if (length(outliers) > 0 && length(outliers) <= 10) {
  cat("Countries:", paste(df_manova$Country[outliers], collapse = ", "), "\n")
}
cat("\n")

# ----------------------------------------------------------------------------
# 1.3 Multivariate Normality by Continent
# ----------------------------------------------------------------------------
cat("--- Multivariate Normality by Continent ---\n")
cat("Using Royston test (extension of Shapiro-Wilk)\n\n")

mvn_by_group <- data.frame(Group = character(), Test = character(),
                           Statistic = numeric(), p_value = character(),
                           MVN = character(), stringsAsFactors = FALSE)

for (cont in levels(df_manova$Continent)) {
  n_cont <- sum(df_manova$Continent == cont)
  cat(cont, " (n = ", n_cont, "):\n", sep = "")

  if (n_cont < 4) {
    cat("  Sample too small for MVN test\n")
    next
  }

  cont_data <- df_manova %>% dplyr::filter(Continent == cont) %>%
    dplyr::select(all_of(dv_vars))

  tryCatch({
    mvn_result <- MVN::mvn(cont_data, mvn_test = "royston")
    print(mvn_result$multivariate_normality)

    mvn_by_group <- rbind(mvn_by_group, data.frame(
      Group = cont,
      Test = "Royston",
      Statistic = as.numeric(mvn_result$multivariate_normality$Statistic),
      p_value = mvn_result$multivariate_normality$p.value,
      MVN = mvn_result$multivariate_normality$MVN
    ))
  }, error = function(e) {
    cat("  MVN test failed:", conditionMessage(e), "\n")
  })
  cat("\n")
}

# MVN Summary Table
cat("\n--- MVN SUMMARY TABLE ---\n")
if (nrow(mvn_by_group) > 0) {
  print(mvn_by_group)
  write.csv(mvn_by_group, "figures/manova_mvn_by_continent.csv", row.names = FALSE)

  n_violated <- sum(grepl("Not normal", mvn_by_group$MVN, ignore.case = TRUE))
  n_satisfied <- sum(grepl("Normal", mvn_by_group$MVN, ignore.case = TRUE) &
                     !grepl("Not normal", mvn_by_group$MVN, ignore.case = TRUE))

  cat("\n--- MVN ASSUMPTION STATUS ---\n")
  cat("Groups satisfying MVN:", n_satisfied, "/", nrow(mvn_by_group), "\n")
  cat("Groups violating MVN:", n_violated, "/", nrow(mvn_by_group), "\n\n")

  if (n_violated > 0) {
    cat("WARNING: MVN assumption is violated in", n_violated, "continent group(s).\n")
    cat("RECOMMENDATIONS:\n")
    cat("  1. Use Pillai's Trace (most robust to MVN violations)\n")
    cat("  2. Sample sizes are adequate (all n > p), which provides some protection\n")
    cat("  3. Consider permutation MANOVA as a non-parametric alternative\n")
    cat("  4. Interpret results with caution\n\n")

    # Flag for later use
    mvn_violated <- TRUE
  } else {
    cat("MVN assumption is satisfied in all groups.\n\n")
    mvn_violated <- FALSE
  }
} else {
  cat("No MVN results available.\n")
  mvn_violated <- TRUE
}

# ----------------------------------------------------------------------------
# 1.4 Box's M Test for Homogeneity of Covariance Matrices
# ----------------------------------------------------------------------------
cat("--- Box's M Test (Covariance Homogeneity) ---\n")
X_matrix <- as.matrix(df_manova[, dv_vars])

tryCatch({
  box_m <- biotools::boxM(X_matrix, df_manova$Continent)
  print(box_m)

  cat("\n--- Box's M Interpretation ---\n")
  if (is.infinite(box_m$statistic) || is.na(box_m$p.value)) {
    cat("Box's M returned Inf/NA - singularity in covariance matrix\n")
    cat("Recommendation: Use Pillai's Trace (most robust)\n")
  } else if (box_m$p.value < 0.05) {
    cat("Box's M significant (p < 0.05)\n")
    cat("Covariance matrices NOT equal across continents\n")
    cat("Recommendation: Use Pillai's Trace (most robust)\n")
  } else {
    cat("Box's M not significant - homogeneity assumption met\n")
  }
}, error = function(e) {
  cat("Box's M test failed:", conditionMessage(e), "\n")
  cat("Recommendation: Use Pillai's Trace\n")
})
cat("\n")

# ============================================================================
# SECTION 2: DESCRIPTIVE STATISTICS BY CONTINENT
# ============================================================================
cat("=== GROUP MEANS BY CONTINENT ===\n")
means_table <- df_manova %>%
  group_by(Continent) %>%
  summarise(
    n = n(),
    across(all_of(dv_vars), ~round(mean(.), 2))
  )
print(means_table, n = 6, width = Inf)
write.csv(means_table, "figures/manova_continent_means.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# 2.1 Bonferroni Confidence Intervals by Continent
# ----------------------------------------------------------------------------
# Simultaneous confidence intervals adjusted for multiple comparisons
# Alpha is divided by number of groups (k) and variables (p)

cat("\n=== BONFERRONI CONFIDENCE INTERVALS ===\n")

k <- length(levels(df_manova$Continent))  # 6 groups
p_vars <- length(dv_vars)                 # 8 variables
alpha <- 0.05
alpha_bonf <- alpha / (k * p_vars)        # Bonferroni-adjusted alpha

cat("Number of groups (k):", k, "\n")
cat("Number of variables (p):", p_vars, "\n")
cat("Bonferroni alpha:", round(alpha_bonf, 5), "\n")
cat("Critical t-value:", round(qt(1 - alpha_bonf/2, df = nrow(df_manova) - k), 3), "\n\n")

# Calculate CI for each continent and variable
ci_data <- df_manova %>%
  group_by(Continent) %>%
  summarise(
    n = n(),
    across(all_of(dv_vars), list(
      mean = ~mean(.),
      se = ~sd(.)/sqrt(n())
    ))
  )

# Create Bonferroni CI plot for key variables
t_crit <- qt(1 - alpha_bonf/2, df = nrow(df_manova) - k)

# Reshape for plotting
ci_long <- df_manova %>%
  group_by(Continent) %>%
  summarise(
    n = n(),
    across(all_of(dv_vars), list(mean = mean, sd = sd))
  ) %>%
  pivot_longer(
    cols = -c(Continent, n),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(mean|sd)"
  ) %>%
  mutate(
    se = sd / sqrt(n),
    ci_lower = mean - t_crit * se,
    ci_upper = mean + t_crit * se
  )

# Plot 1: Bonferroni CI for Life Expectancy
p_ci_life <- ggplot(ci_long %>% filter(Variable == "Life_expectancy"),
                    aes(x = reorder(Continent, mean), y = mean)) +
  geom_point(size = 4, color = "steelblue") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 1) +
  coord_flip() +
  labs(title = "Bonferroni 95% CI: Life Expectancy by Continent",
       subtitle = paste0("Adjusted for ", k * p_vars, " comparisons (α = ", round(alpha_bonf, 5), ")"),
       x = "", y = "Life Expectancy (years)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("figures/bonferroni_ci_life_expectancy.png", p_ci_life, width = 10, height = 6, dpi = 150)
cat("Saved: figures/bonferroni_ci_life_expectancy.png\n")

# Plot 2: Bonferroni CI for all DVs (faceted)
# Standardize variables for comparison
ci_long_std <- df_manova %>%
  mutate(across(all_of(dv_vars), scale)) %>%
  group_by(Continent) %>%
  summarise(
    n = n(),
    across(all_of(dv_vars), list(mean = mean, sd = sd))
  ) %>%
  pivot_longer(
    cols = -c(Continent, n),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(mean|sd)"
  ) %>%
  mutate(
    se = sd / sqrt(n),
    ci_lower = mean - t_crit * se,
    ci_upper = mean + t_crit * se
  )

p_ci_all <- ggplot(ci_long_std, aes(x = Continent, y = mean, color = Continent)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Variable, scales = "free_y", ncol = 4) +
  labs(title = "Bonferroni 95% CI by Continent (Standardized)",
       subtitle = paste0("Dashed line = overall mean; α = ", round(alpha_bonf, 5)),
       x = "", y = "Standardized Mean") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave("figures/bonferroni_ci_all_variables.png", p_ci_all, width = 14, height = 8, dpi = 150)
cat("Saved: figures/bonferroni_ci_all_variables.png\n")

# Save CI data
write.csv(ci_long, "figures/bonferroni_ci_data.csv", row.names = FALSE)

# ============================================================================
# SECTION 3: MANOVA TEST
# ============================================================================
# MANOVA tests H0: mu_Africa = mu_Asia = ... = mu_South_America
#
# Four test statistics:
#   1. Pillai's Trace: Most robust to assumption violations
#   2. Wilks' Lambda: Most common, likelihood ratio test
#   3. Hotelling-Lawley: Sum of eigenvalues
#   4. Roy's Largest Root: Based on first eigenvalue only

cat("\n=== MANOVA RESULTS ===\n")
dv_formula <- as.formula(paste("cbind(", paste(dv_vars, collapse = ", "), ") ~ Continent"))
manova_fit <- manova(dv_formula, data = df_manova)

cat("\n1. PILLAI'S TRACE:\n")
pillai <- summary(manova_fit, test = "Pillai")
print(pillai)

cat("\n2. WILKS' LAMBDA:\n")
wilks <- summary(manova_fit, test = "Wilks")
print(wilks)

cat("\n3. HOTELLING-LAWLEY:\n")
hotelling <- summary(manova_fit, test = "Hotelling-Lawley")
print(hotelling)

cat("\n4. ROY'S LARGEST ROOT:\n")
roy <- summary(manova_fit, test = "Roy")
print(roy)

# Create summary table
test_summary <- data.frame(
  Test = c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
  Statistic = round(c(pillai$stats[1,2], wilks$stats[1,2],
                      hotelling$stats[1,2], roy$stats[1,2]), 4),
  F = round(c(pillai$stats[1,3], wilks$stats[1,3],
              hotelling$stats[1,3], roy$stats[1,3]), 2),
  df1 = c(pillai$stats[1,4], wilks$stats[1,4],
          hotelling$stats[1,4], roy$stats[1,4]),
  df2 = c(pillai$stats[1,5], wilks$stats[1,5],
          hotelling$stats[1,5], roy$stats[1,5]),
  p_value = format(c(pillai$stats[1,6], wilks$stats[1,6],
                     hotelling$stats[1,6], roy$stats[1,6]),
                   scientific = TRUE, digits = 3)
)
cat("\n=== TEST SUMMARY ===\n")
print(test_summary)
write.csv(test_summary, "figures/manova_test_statistics.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# 3.1 Permutation MANOVA (Robust Alternative)
# ----------------------------------------------------------------------------
# When MVN is violated, permutation-based MANOVA provides a non-parametric
# alternative that doesn't assume multivariate normality.
# Uses vegan::adonis2 (PERMANOVA) based on Euclidean distances.

cat("\n=== PERMUTATION MANOVA (ROBUST ALTERNATIVE) ===\n")

if (mvn_violated) {
  cat("Since MVN assumption is violated, running permutation MANOVA...\n\n")
} else {
  cat("MVN satisfied, but running permutation MANOVA for comparison...\n\n")
}

# Check if vegan is available, if not install it
if (!requireNamespace("vegan", quietly = TRUE)) {
  cat("Installing vegan package for permutation MANOVA...\n")
  install.packages("vegan", repos = "https://cloud.r-project.org")
}
library(vegan)

# Permutation MANOVA using adonis2
# Uses Euclidean distance matrix and permutes group labels
set.seed(123)
perm_manova <- adonis2(X_matrix ~ df_manova$Continent,
                        method = "euclidean",
                        permutations = 999)

cat("PERMANOVA (Permutation-based MANOVA):\n")
cat("  Method: Euclidean distance\n")
cat("  Permutations: 999\n\n")
print(perm_manova)

# Extract results
perm_f <- perm_manova$F[1]
perm_r2 <- perm_manova$R2[1]
perm_p <- perm_manova$`Pr(>F)`[1]

cat("\n--- PERMANOVA Summary ---\n")
cat("F-statistic:", round(perm_f, 2), "\n")
cat("R² (variance explained):", round(perm_r2, 3), "\n")
cat("p-value:", format(perm_p, scientific = TRUE, digits = 3), "\n")

if (perm_p < 0.05) {
  cat("\nConclusion: Significant difference between continents (p < 0.05)\n")
  cat("The permutation test CONFIRMS the parametric MANOVA results.\n")
} else {
  cat("\nConclusion: No significant difference (p >= 0.05)\n")
}

# Add to test summary
test_summary <- rbind(test_summary, data.frame(
  Test = "PERMANOVA",
  Statistic = round(perm_r2, 4),
  F = round(perm_f, 2),
  df1 = perm_manova$Df[1],
  df2 = perm_manova$Df[2],
  p_value = format(perm_p, scientific = TRUE, digits = 3)
))
write.csv(test_summary, "figures/manova_test_statistics.csv", row.names = FALSE)

cat("\n")

# ============================================================================
# SECTION 4: UNIVARIATE FOLLOW-UP ANOVAS
# ============================================================================
cat("\n=== UNIVARIATE ANOVAs ===\n")
anova_summary <- summary.aov(manova_fit)

p <- length(dv_vars)
alpha_bonf <- 0.05 / p

univar_results <- data.frame(Variable = dv_vars, F = NA, p = NA, eta_sq = NA, Sig = NA)
for (i in 1:p) {
  aov_res <- anova_summary[[i]]
  univar_results$F[i] <- round(aov_res["Continent", "F value"], 2)
  univar_results$p[i] <- aov_res["Continent", "Pr(>F)"]
  univar_results$eta_sq[i] <- round(aov_res["Continent", "Sum Sq"] /
                                     sum(aov_res[, "Sum Sq"]), 3)
  univar_results$Sig[i] <- ifelse(aov_res["Continent", "Pr(>F)"] < alpha_bonf, "*", "")
}
print(univar_results)
write.csv(univar_results, "figures/univariate_anovas.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# 4.1 Kruskal-Wallis Tests (Non-parametric Alternative)
# ----------------------------------------------------------------------------
# When MVN is violated, Kruskal-Wallis provides a non-parametric alternative
# to one-way ANOVA for each dependent variable separately.
# Tests H0: All group medians are equal

cat("\n=== KRUSKAL-WALLIS TESTS (Non-parametric) ===\n")
cat("Alternative to univariate ANOVAs when normality is violated\n\n")

kw_results <- data.frame(Variable = dv_vars, Chi_sq = NA, df = NA, p_value = NA, Sig = NA)
for (i in 1:length(dv_vars)) {
  kw_test <- kruskal.test(df_manova[[dv_vars[i]]] ~ df_manova$Continent)
  kw_results$Chi_sq[i] <- round(kw_test$statistic, 2)
  kw_results$df[i] <- kw_test$parameter
  kw_results$p_value[i] <- kw_test$p.value
  kw_results$Sig[i] <- ifelse(kw_test$p.value < alpha_bonf, "*", "")
}

print(kw_results)
cat("* = Significant after Bonferroni correction (alpha =", round(alpha_bonf, 4), ")\n")
write.csv(kw_results, "figures/kruskal_wallis_tests.csv", row.names = FALSE)

# Compare parametric vs non-parametric results
cat("\n--- Parametric vs Non-parametric Comparison ---\n")
comparison <- data.frame(
  Variable = dv_vars,
  ANOVA_p = format(univar_results$p, scientific = TRUE, digits = 2),
  KW_p = format(kw_results$p_value, scientific = TRUE, digits = 2),
  Both_Sig = ifelse(univar_results$Sig == "*" & kw_results$Sig == "*", "Yes", "No")
)
print(comparison)

n_agree <- sum(comparison$Both_Sig == "Yes")
cat("\nAgreement:", n_agree, "/", length(dv_vars), "variables significant in both tests\n")

if (n_agree == length(dv_vars)) {
  cat("All parametric and non-parametric tests agree - results are ROBUST\n")
}

# ============================================================================
# SECTION 5: POST-HOC PAIRWISE COMPARISONS
# ============================================================================
# With 6 groups, we need to identify WHICH continents differ
# Using Bonferroni-corrected t-tests

cat("\n=== POST-HOC PAIRWISE COMPARISONS (Bonferroni) ===\n")

posthoc_list <- list()
for (var in dv_vars) {
  cat("\n---", var, "---\n")
  pairwise_result <- pairwise.t.test(df_manova[[var]], df_manova$Continent,
                                      p.adjust.method = "bonferroni")
  print(pairwise_result)

  # Store p-values for heatmap
  posthoc_list[[var]] <- pairwise_result$p.value
}

# Save first variable's pairwise comparison as example
write.csv(posthoc_list[[1]], "figures/manova_pairwise_life_expectancy.csv")

# Create summary of significant pairwise differences
cat("\n=== SIGNIFICANT PAIRWISE DIFFERENCES (p < 0.05) ===\n")
sig_pairs <- data.frame()
for (var in dv_vars) {
  pmat <- posthoc_list[[var]]
  for (i in 1:nrow(pmat)) {
    for (j in 1:ncol(pmat)) {
      if (!is.na(pmat[i, j]) && pmat[i, j] < 0.05) {
        sig_pairs <- rbind(sig_pairs, data.frame(
          Variable = var,
          Continent1 = rownames(pmat)[i],
          Continent2 = colnames(pmat)[j],
          p_value = round(pmat[i, j], 4)
        ))
      }
    }
  }
}
if (nrow(sig_pairs) > 0) {
  print(sig_pairs)
  write.csv(sig_pairs, "figures/manova_pairwise_significant.csv", row.names = FALSE)
} else {
  cat("No significant pairwise differences found\n")
}

# ============================================================================
# SECTION 6: VISUALIZATIONS
# ============================================================================

# Define continent colors
continent_colors <- c("Africa" = "#E41A1C", "Asia" = "#377EB8",
                      "Europe" = "#4DAF4A", "North America" = "#984EA3",
                      "South America" = "#FF7F00", "Oceania" = "#FFFF33")

# ----------------------------------------------------------------------------
# 6.1 Profile Plot (Standardized Means by Continent)
# ----------------------------------------------------------------------------
cat("\n=== Creating visualizations ===\n")

df_profile <- df_manova %>%
  mutate(across(all_of(dv_vars), ~scale(.)[,1])) %>%
  pivot_longer(cols = all_of(dv_vars), names_to = "Variable", values_to = "Value")

profile_summary <- df_profile %>%
  group_by(Continent, Variable) %>%
  summarise(Mean = mean(Value),
            SE = sd(Value)/sqrt(n()),
            .groups = "drop")

p_profile <- ggplot(profile_summary,
                    aes(x = Variable, y = Mean, group = Continent, color = Continent)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = Mean - 1.96*SE, ymax = Mean + 1.96*SE, fill = Continent),
              alpha = 0.1, color = NA) +
  scale_color_manual(values = continent_colors) +
  scale_fill_manual(values = continent_colors) +
  labs(title = "Mean Profile Plot by Continent (+/- 95% CI)",
       y = "Standardized Mean",
       x = "Variable") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"))

ggsave("figures/manova_continent_profiles.png", p_profile, width = 14, height = 8, dpi = 150)
cat("Saved: manova_continent_profiles.png\n")

# ----------------------------------------------------------------------------
# 6.2 Faceted Boxplots by Continent
# ----------------------------------------------------------------------------
df_boxplot <- df_manova %>%
  pivot_longer(cols = all_of(dv_vars), names_to = "Variable", values_to = "Value")

p_boxplots <- ggplot(df_boxplot, aes(x = Continent, y = Value, fill = Continent)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  facet_wrap(~Variable, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = continent_colors) +
  labs(title = "Distribution of Health/Socioeconomic Variables by Continent",
       x = "", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

ggsave("figures/manova_continent_boxplots.png", p_boxplots, width = 16, height = 10, dpi = 150)
cat("Saved: manova_continent_boxplots.png\n")

# ----------------------------------------------------------------------------
# 6.3 Pairwise Comparison Heatmap (Life Expectancy as example)
# ----------------------------------------------------------------------------
library(reshape2)

# Create heatmap for Life Expectancy pairwise comparisons
pmat_le <- posthoc_list[["Life_expectancy"]]
pmat_le_melt <- melt(pmat_le, na.rm = TRUE)
colnames(pmat_le_melt) <- c("Continent1", "Continent2", "p_value")

# Create significance labels
pmat_le_melt$sig <- ifelse(pmat_le_melt$p_value < 0.001, "***",
                           ifelse(pmat_le_melt$p_value < 0.01, "**",
                                  ifelse(pmat_le_melt$p_value < 0.05, "*", "")))

p_heatmap <- ggplot(pmat_le_melt, aes(x = Continent1, y = Continent2, fill = p_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 6, color = "white") +
  scale_fill_gradient2(low = "darkred", mid = "orange", high = "lightgray",
                       midpoint = 0.025, name = "p-value",
                       limits = c(0, 1)) +
  labs(title = "Pairwise Comparisons: Life Expectancy by Continent",
       subtitle = "*** p < 0.001, ** p < 0.01, * p < 0.05 (Bonferroni corrected)",
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

ggsave("figures/manova_pairwise_heatmap.png", p_heatmap, width = 10, height = 8, dpi = 150)
cat("Saved: manova_pairwise_heatmap.png\n")

# ----------------------------------------------------------------------------
# 6.4 Effect Size Plot
# ----------------------------------------------------------------------------
p_effect <- ggplot(univar_results, aes(x = reorder(Variable, eta_sq), y = eta_sq)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = c(0.01, 0.06, 0.14),
             linetype = "dashed", color = c("green", "orange", "red")) +
  geom_text(aes(label = paste0(round(eta_sq * 100, 1), "%")),
            hjust = -0.2, size = 3.5) +
  coord_flip() +
  labs(title = "Effect Sizes: Continent Differences (Eta-squared)",
       subtitle = "Green: small (0.01), Orange: medium (0.06), Red: large (0.14)",
       x = "", y = expression(eta^2)) +
  theme_minimal() +
  ylim(0, max(univar_results$eta_sq) * 1.2)

ggsave("figures/manova_effect_sizes.png", p_effect, width = 10, height = 6, dpi = 150)
cat("Saved: manova_effect_sizes.png\n")

# ----------------------------------------------------------------------------
# 6.5 Centroids Plot (2D Projection by Continent)
# ----------------------------------------------------------------------------
centroids <- df_manova %>%
  group_by(Continent) %>%
  summarise(Life_expectancy = mean(Life_expectancy),
            Schooling = mean(Schooling),
            .groups = "drop")

p_centroids <- ggplot(df_manova,
                      aes(x = Life_expectancy, y = Schooling, color = Continent)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_point(data = centroids, size = 6, shape = 18) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = continent_colors) +
  labs(title = "Continent Centroids with 95% Confidence Ellipses",
       subtitle = "Diamonds = continent centroids") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("figures/manova_centroids.png", p_centroids, width = 12, height = 8, dpi = 150)
cat("Saved: manova_centroids.png\n")

# ============================================================================
# SECTION 7: PERMUTATION MANOVA (ROBUST ALTERNATIVE)
# ============================================================================
cat("\n=== PERMUTATION MANOVA (ROBUST ALTERNATIVE) ===\n")
cat("Permutation tests do not require MVN or equal covariances\n\n")

library(vegan)

X_scaled <- scale(X_matrix)
set.seed(123)
perm_manova <- adonis2(X_scaled ~ Continent, data = df_manova,
                       method = "euclidean", permutations = 999)
cat("Permutation MANOVA (adonis2, 999 permutations):\n")
print(perm_manova)

# Compare with parametric results
pillai_p <- pillai$stats[1, 6]
perm_p <- perm_manova$`Pr(>F)`[1]

cat("\nComparison with parametric MANOVA:\n")
cat("  Parametric p-value (Pillai):  ", format(pillai_p, scientific = TRUE, digits = 3), "\n")
cat("  Permutation p-value:          ", format(perm_p, scientific = TRUE, digits = 3), "\n")

if (perm_p < 0.05 && pillai_p < 0.05) {
  cat("  Both tests significant - results are ROBUST\n")
} else if (perm_p >= 0.05 && pillai_p >= 0.05) {
  cat("  Both tests non-significant - results are ROBUST\n")
} else {
  cat("  Tests disagree - interpret with CAUTION\n")
}

# ============================================================================
# SECTION 8: SUMMARY
# ============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("MANOVA BY CONTINENT - SUMMARY\n")
cat(strrep("=", 60), "\n")
cat("Groups: 6 Continents\n")
cat("Dependent Variables:", length(dv_vars), "\n")
cat("Sample Size:", nrow(df_manova), "countries\n\n")

cat("Pillai's Trace:", round(pillai$stats[1, 2], 4), "\n")
cat("F-statistic:", round(pillai$stats[1, 3], 2), "\n")
cat("p-value:", format(pillai_p, scientific = TRUE, digits = 3), "\n\n")

cat("Conclusion: ")
if (pillai_p < 0.05) {
  cat("Significant multivariate differences exist between continents.\n")
  cat("Post-hoc tests reveal which specific pairs of continents differ.\n")
} else {
  cat("No significant multivariate differences between continents.\n")
}

cat(strrep("=", 60), "\n")
cat("\n=== MANOVA Complete ===\n")
cat("Run 04_pca_regression.R next.\n")

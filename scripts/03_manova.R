# STAT 467 - MANOVA
# Input: data_country_level.csv

library(tidyverse)
library(MASS)
library(car)
library(biotools)
library(MVN)
library(ggpubr)
library(rstatix)
library(heplots)

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n")
cat("Developed:", sum(df$Status == "Developed"), ", Developing:", sum(df$Status == "Developing"), "\n\n")

dv_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "Income_composition",
             "BMI", "HIV_AIDS", "GDP_log", "Diphtheria")

df_manova <- df %>% select(Country, Status, all_of(dv_vars)) %>% drop_na()
cat("Complete cases:", nrow(df_manova), "\n\n")

# Assumptions
cat("=== ASSUMPTIONS ===\n")

# Sample size
n_per_group <- table(df_manova$Status)
print(n_per_group)
cat("n_min > p:", min(n_per_group), ">", length(dv_vars), "\n\n")

# Multivariate outlier detection
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
} else if (length(outliers) > 10) {
  cat("First 10:", paste(df_manova$Country[outliers[1:10]], collapse = ", "), "...\n")
}
cat("\n")

# Multivariate normality by group
# Using Royston test (recommended for n < 5000, based on Shapiro-Wilk)
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
    mvn_royston <- MVN::mvn(grp_data, mvnTest = "royston")
    print(mvn_royston$multivariateNormality)

    mvn_by_group <- rbind(mvn_by_group, data.frame(
      Group = grp,
      Test = "Royston",
      Statistic = as.numeric(mvn_royston$multivariateNormality$Statistic),
      p_value = mvn_royston$multivariateNormality$`p value`,
      MVN = mvn_royston$multivariateNormality$MVN
    ))
  }, error = function(e) {
    cat("Royston test failed:", conditionMessage(e), "\n")
    cat("Trying Mardia test...\n")
    tryCatch({
      mvn_mardia <- MVN::mvn(grp_data, mvnTest = "mardia")
      print(mvn_mardia$multivariateNormality)
    }, error = function(e2) {
      cat("MVN test failed:", conditionMessage(e2), "\n")
    })
  })
}

# Univariate normality by group (Shapiro-Wilk)
cat("\n--- Univariate Normality by Group (Shapiro-Wilk) ---\n")
df_for_sw <- df_manova %>%
  dplyr::select(Status, all_of(dv_vars))
univar_by_group <- df_for_sw %>%
  pivot_longer(cols = -Status, names_to = "Variable", values_to = "Value") %>%
  group_by(Status, Variable) %>%
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

# Box's M test for homogeneity of covariance matrices
cat("--- Box's M Test (Covariance Homogeneity) ---\n")
X_matrix <- as.matrix(df_manova[, dv_vars])

# Using heplots::boxM for better output
box_m <- heplots::boxM(X_matrix, df_manova$Status)
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

# Descriptives
cat("=== GROUP MEANS ===\n")
means_table <- df_manova %>%
  group_by(Status) %>%
  summarise(across(all_of(dv_vars), ~round(mean(.), 2)))
print(means_table)

# MANOVA
cat("\n=== MANOVA RESULTS ===\n")
dv_formula <- as.formula(paste("cbind(", paste(dv_vars, collapse = ", "), ") ~ Status"))
manova_fit <- manova(dv_formula, data = df_manova)

cat("\n1. PILLAI'S TRACE:\n")
print(summary(manova_fit, test = "Pillai"))

cat("\n2. WILKS' LAMBDA:\n")
print(summary(manova_fit, test = "Wilks"))

cat("\n3. HOTELLING-LAWLEY:\n")
print(summary(manova_fit, test = "Hotelling-Lawley"))

cat("\n4. ROY'S LARGEST ROOT:\n")
print(summary(manova_fit, test = "Roy"))

# Summary table
pillai <- summary(manova_fit, test = "Pillai")
wilks <- summary(manova_fit, test = "Wilks")
hotelling <- summary(manova_fit, test = "Hotelling-Lawley")
roy <- summary(manova_fit, test = "Roy")

test_summary <- data.frame(
  Test = c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
  Statistic = round(c(pillai$stats[1,2], wilks$stats[1,2], hotelling$stats[1,2], roy$stats[1,2]), 4),
  F = round(c(pillai$stats[1,3], wilks$stats[1,3], hotelling$stats[1,3], roy$stats[1,3]), 2),
  p_value = format(c(pillai$stats[1,6], wilks$stats[1,6], hotelling$stats[1,6], roy$stats[1,6]), scientific = TRUE, digits = 3)
)
cat("\n=== TEST SUMMARY ===\n")
print(test_summary)
write.csv(test_summary, "figures/manova_test_statistics.csv", row.names = FALSE)

# Univariate ANOVAs
cat("\n=== UNIVARIATE ANOVAs ===\n")
anova_summary <- summary.aov(manova_fit)

p <- length(dv_vars)
alpha_bonf <- 0.05 / p

univar_results <- data.frame(Variable = dv_vars, F = NA, p = NA, eta_sq = NA, Sig = NA)
for (i in 1:p) {
  aov_res <- anova_summary[[i]]
  univar_results$F[i] <- round(aov_res["Status", "F value"], 2)
  univar_results$p[i] <- aov_res["Status", "Pr(>F)"]
  univar_results$eta_sq[i] <- round(aov_res["Status", "Sum Sq"] / sum(aov_res[, "Sum Sq"]), 3)
  univar_results$Sig[i] <- ifelse(aov_res["Status", "Pr(>F)"] < alpha_bonf, "*", "")
}
print(univar_results)
write.csv(univar_results, "figures/univariate_anovas.csv", row.names = FALSE)

# Profile plot
df_profile <- df_manova %>%
  mutate(across(all_of(dv_vars), ~scale(.)[,1])) %>%
  pivot_longer(cols = all_of(dv_vars), names_to = "Variable", values_to = "Value")

profile_summary <- df_profile %>%
  group_by(Status, Variable) %>%
  summarise(Mean = mean(Value), SE = sd(Value)/sqrt(n()), .groups = "drop")

p_profile <- ggplot(profile_summary, aes(x = Variable, y = Mean, group = Status, color = Status)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  scale_color_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
  labs(title = "Mean Profile Plot by Status", y = "Standardized Mean") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/manova_profile_plot.png", p_profile, width = 12, height = 7, dpi = 150)

# Effect size plot
p_effect <- ggplot(univar_results, aes(x = reorder(Variable, eta_sq), y = eta_sq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = c(0.01, 0.06, 0.14), linetype = "dashed", color = c("green", "orange", "red")) +
  coord_flip() +
  labs(title = "Effect Sizes (Eta-squared)", x = "", y = expression(eta^2)) +
  theme_minimal()
ggsave("figures/manova_effect_sizes.png", p_effect, width = 10, height = 6, dpi = 150)

# Centroids plot
centroids <- df_manova %>%
  group_by(Status) %>%
  summarise(Life_expectancy = mean(Life_expectancy), Schooling = mean(Schooling))

p_centroids <- ggplot(df_manova, aes(x = Life_expectancy, y = Schooling, color = Status)) +
  geom_point(alpha = 0.5) +
  geom_point(data = centroids, size = 5, shape = 18) +
  stat_ellipse(level = 0.95, linewidth = 1.2) +
  scale_color_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
  labs(title = "Group Centroids with 95% Ellipses") +
  theme_minimal()
ggsave("figures/manova_centroids.png", p_centroids, width = 10, height = 8, dpi = 150)

# Discriminant coefficients
# CRITICAL FIX: Calculate mean difference from RAW DATA, not rounded means_table
n1 <- sum(df_manova$Status == "Developed")
n2 <- sum(df_manova$Status == "Developing")

# Get group data
grp1_data <- df_manova %>% filter(Status == "Developed") %>% select(all_of(dv_vars))
grp2_data <- df_manova %>% filter(Status == "Developing") %>% select(all_of(dv_vars))

# Calculate covariance matrices
S1 <- cov(grp1_data)
S2 <- cov(grp2_data)

# Pooled covariance matrix
Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)

# Mean difference from RAW data (not rounded means_table)
x_bar1 <- colMeans(grp1_data)
x_bar2 <- colMeans(grp2_data)
mean_diff <- x_bar1 - x_bar2

# Raw discriminant coefficients: a = Sp^{-1} * (x_bar1 - x_bar2)
disc_coef <- solve(Sp) %*% mean_diff

# Standardized discriminant coefficients: a_std = a * sqrt(diag(Sp))
# This gives coefficients for standardized variables (mean=0, sd=1)
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

cat("\n=== MANOVA Complete ===\n")
cat("Run 04_pca_regression.R next.\n")




# STAT 467 - Mean Vector Inference
# Input: data_country_level.csv

library(tidyverse)
library(MASS)
library(Hotelling)
library(MVN)
library(ellipse)
library(rstatix)
library(heplots)

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n")
cat("Developed:", sum(df$Status == "Developed"), ", Developing:", sum(df$Status == "Developing"), "\n\n")

selected_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log", "BMI", "HIV_AIDS")

X <- df %>% select(all_of(selected_vars)) %>% as.matrix()
n <- nrow(X)
p <- ncol(X)

x_bar <- colMeans(X)
S <- cov(X)

cat("=== Sample Mean Vector ===\n")
print(round(x_bar, 3))

# ============================================================================
# ASSUMPTION CHECKING FOR HOTELLING'S T²
# ============================================================================

# 1. MULTIVARIATE NORMALITY (Overall)
cat("\n=== MULTIVARIATE NORMALITY CHECK ===\n")

# Use Royston test for n < 5000 (recommended by Güler, 2026)
mvn_royston <- MVN::mvn(as.data.frame(X), mvnTest = "royston")
cat("Royston's MVN Test (overall):\n")
print(mvn_royston$multivariateNormality)

# Mardia's test for skewness and kurtosis
mvn_mardia <- MVN::mvn(as.data.frame(X), mvnTest = "mardia")
cat("\nMardia's MVN Test (overall):\n")
print(mvn_mardia$multivariateNormality)

# 2. UNIVARIATE NORMALITY
cat("\n--- Univariate Normality (Shapiro-Wilk) ---\n")
univar_df <- as.data.frame(X) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  shapiro_test(Value)
univar_df$Normal <- ifelse(univar_df$p > 0.05, "Yes", "No")
print(univar_df)

# Q-Q plot using MVN package
png("figures/mean_vector_qq_plot.png", width = 800, height = 600, res = 120)
plot(mvn_mardia, diagnostic = "multivariate", type = "qq")
dev.off()

# One-sample Hotelling's T² test
cat("\n=== ONE-SAMPLE HOTELLING'S T² TEST ===\n")

mu0 <- c(Life_expectancy = 70, Adult_Mortality = 150, Schooling = 12,
         GDP_log = 8.5, BMI = 25, HIV_AIDS = 1.0)

cat("Hypothesized mean:\n")
print(mu0)

diff <- x_bar - mu0
T2 <- n * t(diff) %*% solve(S) %*% diff
F_stat <- ((n - p) / (p * (n - 1))) * T2
p_value <- 1 - pf(F_stat, p, n - p)

cat("\nT² =", round(T2, 4), ", F =", round(F_stat, 4), ", p =", format(p_value, scientific = TRUE), "\n")

if (p_value < 0.05) {
  cat("Reject H0: Mean vector differs from hypothesized values\n")
} else {
  cat("Fail to reject H0\n")
}

# Bonferroni confidence intervals
cat("\n=== BONFERRONI CONFIDENCE INTERVALS ===\n")
alpha <- 0.05
t_crit <- qt(1 - alpha/(2*p), df = n - 1)

ci_bonf <- data.frame(
  Variable = selected_vars,
  Mean = round(x_bar, 3),
  Lower = round(x_bar - t_crit * sqrt(diag(S)/n), 3),
  Upper = round(x_bar + t_crit * sqrt(diag(S)/n), 3),
  H0 = mu0
)
print(ci_bonf)

# Confidence ellipse plot
png("figures/confidence_ellipse.png", width = 900, height = 700, res = 120)
X_2d <- df %>% select(Life_expectancy, Schooling) %>% as.matrix()
x_bar_2d <- colMeans(X_2d)
S_2d <- cov(X_2d)
c2 <- (2 * (n - 1) / (n - 2)) * qf(0.95, 2, n - 2)

theta <- seq(0, 2*pi, length.out = 100)
eigen_S <- eigen(S_2d)
axes <- sqrt(c2/n) * sqrt(eigen_S$values)
rotation <- eigen_S$vectors
ellipse_pts <- t(rotation %*% rbind(axes[1]*cos(theta), axes[2]*sin(theta))) + matrix(x_bar_2d, nrow = 100, ncol = 2, byrow = TRUE)

plot(X_2d, pch = 19, col = "gray50", xlab = "Life Expectancy", ylab = "Schooling",
     main = "95% Confidence Ellipse for Mean Vector")
lines(ellipse_pts, col = "steelblue", lwd = 3)
points(x_bar_2d[1], x_bar_2d[2], pch = 19, col = "red", cex = 2)
points(mu0["Life_expectancy"], mu0["Schooling"], pch = 4, col = "darkgreen", cex = 2, lwd = 3)
legend("bottomright", c("Countries", "Sample Mean", "Hypothesized Mean", "95% CR"),
       col = c("gray50", "red", "darkgreen", "steelblue"), pch = c(19, 19, 4, NA), lty = c(NA, NA, NA, 1), lwd = c(NA, NA, 3, 3))
dev.off()

# ============================================================================
# TWO-SAMPLE HOTELLING'S T² TEST
# ============================================================================
cat("\n=== TWO-SAMPLE HOTELLING'S T² TEST ===\n")

X_dev <- df %>% dplyr::filter(Status == "Developed") %>%
  dplyr::select(all_of(selected_vars)) %>% as.matrix()
X_devp <- df %>% dplyr::filter(Status == "Developing") %>%
  dplyr::select(all_of(selected_vars)) %>% as.matrix()

n1 <- nrow(X_dev)
n2 <- nrow(X_devp)

x_bar1 <- colMeans(X_dev)
x_bar2 <- colMeans(X_devp)

cat("Sample sizes: n1 =", n1, ", n2 =", n2, "\n\n")

# ASSUMPTION 1: Multivariate Normality by Group
cat("--- Assumption 1: MVN by Group ---\n")

# Developed group
cat("\nDeveloped (n =", n1, "):\n")
if (n1 >= 7) {
  mvn_dev <- MVN::mvn(as.data.frame(X_dev), mvnTest = "royston")
  print(mvn_dev$multivariateNormality)
} else {
  cat("Sample too small for MVN test (n < 7). Checking univariate normality:\n")
  for (v in selected_vars) {
    sw_test <- shapiro.test(X_dev[, v])
    cat("  ", v, ": p =", format(sw_test$p.value, digits = 3), "\n")
  }
}

# Developing group
cat("\nDeveloping (n =", n2, "):\n")
mvn_devp <- MVN::mvn(as.data.frame(X_devp), mvnTest = "royston")
print(mvn_devp$multivariateNormality)

# Univariate normality by group (using rstatix)
cat("\n--- Univariate Normality by Group (Shapiro-Wilk) ---\n")
df_for_test <- df %>%
  dplyr::select(Status, all_of(selected_vars))
univar_by_group <- df_for_test %>%
  pivot_longer(cols = -Status, names_to = "Variable", values_to = "Value") %>%
  group_by(Status, Variable) %>%
  shapiro_test(Value)
univar_by_group$Normal <- ifelse(univar_by_group$p > 0.05, "Yes", "No")
print(univar_by_group)

# ASSUMPTION 2: Homogeneity of Covariance Matrices (Box's M test)
cat("\n--- Assumption 2: Box's M Test (Covariance Homogeneity) ---\n")
X_combined <- rbind(X_dev, X_devp)
groups <- factor(c(rep("Developed", n1), rep("Developing", n2)))
box_m <- heplots::boxM(X_combined, groups)
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

cat("\nMean difference (Developed - Developing):\n")
print(round(x_bar1 - x_bar2, 3))

# Pooled covariance matrix
Sp <- ((n1-1)*cov(X_dev) + (n2-1)*cov(X_devp)) / (n1 + n2 - 2)
diff_means <- x_bar1 - x_bar2
T2_two <- (n1*n2/(n1+n2)) * t(diff_means) %*% solve(Sp) %*% diff_means
F_two <- ((n1 + n2 - p - 1) / (p * (n1 + n2 - 2))) * T2_two
p_two <- 1 - pf(F_two, p, n1 + n2 - p - 1)

cat("\nT² =", round(T2_two, 4), ", F =", round(F_two, 4), ", p =", format(p_two, scientific = TRUE), "\n")

if (p_two < 0.05) {
  cat("Reject H0: Developed and Developing mean vectors differ significantly\n")
}

# Verification
cat("\n=== Hotelling Package Verification ===\n")
print(hotelling.test(X_dev, X_devp))

# Univariate follow-up
cat("\n=== UNIVARIATE t-TESTS ===\n")
for (i in 1:p) {
  tt <- t.test(X_dev[,i], X_devp[,i], var.equal = TRUE)
  sig <- ifelse(tt$p.value < alpha/p, "*", "")
  cat(selected_vars[i], ": t =", round(tt$statistic, 2), ", p =", format(tt$p.value, digits = 3), sig, "\n")
}

cat("\n=== Mean Vector Inference Complete ===\n")
cat("Run 03_manova.R next.\n")

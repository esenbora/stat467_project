# STAT 467 - Canonical Correlation Analysis (Bonus)
# Input: data_country_level.csv

library(tidyverse)
library(CCA)
library(corrplot)
library(gridExtra)

# Install CCP if not available
if (!require(CCP, quietly = TRUE)) {
  install.packages("CCP", repos = "https://cloud.r-project.org")
  library(CCP)
}

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n\n")

# Variable sets
socio_vars <- c("GDP_log", "Schooling", "Income_composition", "Total_expenditure", "Alcohol")
health_vars <- c("Life_expectancy", "Adult_Mortality", "infant_deaths", "HIV_AIDS", "BMI")

cat("Set X (Socioeconomic):", paste(socio_vars, collapse = ", "), "\n")
cat("Set Y (Health):", paste(health_vars, collapse = ", "), "\n")

df_cca <- df %>% select(Country, Status, all_of(socio_vars), all_of(health_vars)) %>% drop_na()
cat("Complete cases:", nrow(df_cca), "\n\n")

X <- as.matrix(df_cca[, socio_vars])
Y <- as.matrix(df_cca[, health_vars])
X_std <- scale(X)
Y_std <- scale(Y)

# Correlations
cat("=== CROSS-CORRELATIONS ===\n")
R_xy <- cor(X, Y)
print(round(R_xy, 3))

png("figures/cca_cross_correlations.png", width = 800, height = 600, res = 120)
corrplot(R_xy, method = "color", addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", col = colorRampPalette(c("#E94F37", "white", "#2E86AB"))(200),
         title = "Cross-Correlations", mar = c(0, 0, 2, 0))
dev.off()

# CCA
cat("\n=== CCA ===\n")
cca_result <- cc(X_std, Y_std)
n_cc <- min(ncol(X), ncol(Y))
rho <- cca_result$cor

cc_df <- data.frame(CV = paste0("CV", 1:n_cc), r = round(rho, 4), r_sq = round(rho^2, 4))
print(cc_df)

# Significance
cat("\n=== SIGNIFICANCE ===\n")
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)

wilks_test <- p.asym(rho, n, p, q, tstat = "Wilks")
print(wilks_test)

lambda1 <- prod(1 - rho^2)
chi_sq <- -(n - 1 - (p + q + 1)/2) * log(lambda1)
p_val <- 1 - pchisq(chi_sq, p * q)
cat("\nOverall Wilks:", round(lambda1, 4), ", Chi-sq:", round(chi_sq, 2), ", p:", format(p_val, scientific = TRUE), "\n")

# CRITICAL FIX: Sequential tests for each canonical correlation
# Tests H0: rho_i = rho_{i+1} = ... = rho_s = 0 (all remaining correlations are zero)
cat("\n=== SEQUENTIAL DIMENSION TESTS ===\n")
cat("Testing significance of each canonical correlation:\n\n")

seq_results <- data.frame(
  CV = paste0("CV", 1:n_cc),
  r = round(rho, 4),
  r_sq = round(rho^2, 4),
  Wilks_Lambda = NA,
  Chi_sq = NA,
  df = NA,
  p_value = NA,
  Significant = NA
)

for (i in 1:n_cc) {
  # Wilks' Lambda for dimensions i through s
  lambda_i <- prod(1 - rho[i:n_cc]^2)

  # Degrees of freedom for this test
  df_i <- (p - i + 1) * (q - i + 1)

  # Bartlett's chi-square approximation
  chi_i <- -(n - 1 - (p + q + 1)/2) * log(lambda_i)

  # P-value
  p_val_i <- 1 - pchisq(chi_i, df_i)

  seq_results$Wilks_Lambda[i] <- round(lambda_i, 4)
  seq_results$Chi_sq[i] <- round(chi_i, 2)
  seq_results$df[i] <- df_i
  seq_results$p_value[i] <- p_val_i
  seq_results$Significant[i] <- ifelse(p_val_i < 0.05, "*", "")
}

print(seq_results)
write.csv(seq_results, "figures/cca_sequential_tests.csv", row.names = FALSE)

# Count significant canonical correlations
n_sig <- sum(seq_results$p_value < 0.05)
cat("\nNumber of significant canonical correlations (alpha=0.05):", n_sig, "out of", n_cc, "\n")

# Coefficients
cat("\n=== COEFFICIENTS ===\n")
coef_X <- cca_result$xcoef
coef_Y <- cca_result$ycoef
colnames(coef_X) <- colnames(coef_Y) <- paste0("CV", 1:n_cc)
rownames(coef_X) <- socio_vars
rownames(coef_Y) <- health_vars

cat("X coefficients:\n")
print(round(coef_X, 4))
cat("\nY coefficients:\n")
print(round(coef_Y, 4))

write.csv(round(coef_X, 4), "figures/cca_coef_X.csv")
write.csv(round(coef_Y, 4), "figures/cca_coef_Y.csv")

# Canonical variates
U <- X_std %*% coef_X
V <- Y_std %*% coef_Y

# Loadings
loadings_X <- cor(X_std, U)
loadings_Y <- cor(Y_std, V)
colnames(loadings_X) <- paste0("U", 1:n_cc)
colnames(loadings_Y) <- paste0("V", 1:n_cc)

cat("\n=== LOADINGS ===\n")
cat("X on U:\n")
print(round(loadings_X, 3))
cat("\nY on V:\n")
print(round(loadings_Y, 3))

write.csv(round(loadings_X, 4), "figures/cca_loadings_X.csv")
write.csv(round(loadings_Y, 4), "figures/cca_loadings_Y.csv")

# LOADING INTERPRETATION: Identify important variables (|loading| >= 0.4)
cat("\n=== LOADING INTERPRETATION ===\n")
cat("Variables with |loading| >= 0.4 are considered important contributors\n\n")

loading_threshold <- 0.4

for (cv_idx in 1:min(n_sig, 3)) {  # Focus on significant CVs, max 3
  cat("--- CV", cv_idx, "(r =", round(rho[cv_idx], 3), ") ---\n")

  # Important X variables
  x_important <- which(abs(loadings_X[, cv_idx]) >= loading_threshold)
  if (length(x_important) > 0) {
    cat("Socioeconomic (X) variables:\n")
    for (idx in x_important) {
      cat("  ", rownames(loadings_X)[idx], ": ", round(loadings_X[idx, cv_idx], 3), "\n", sep = "")
    }
  } else {
    cat("Socioeconomic (X): No variables meet threshold\n")
  }

  # Important Y variables
  y_important <- which(abs(loadings_Y[, cv_idx]) >= loading_threshold)
  if (length(y_important) > 0) {
    cat("Health (Y) variables:\n")
    for (idx in y_important) {
      cat("  ", rownames(loadings_Y)[idx], ": ", round(loadings_Y[idx, cv_idx], 3), "\n", sep = "")
    }
  } else {
    cat("Health (Y): No variables meet threshold\n")
  }
  cat("\n")
}

# Redundancy Analysis
# Redundancy = (variance in own set explained by own CV) * (shared variance via canonical r^2)
var_X_U <- colMeans(loadings_X^2)  # Avg variance in X explained by each U
var_Y_V <- colMeans(loadings_Y^2)  # Avg variance in Y explained by each V
rd_X_Y <- var_X_U * rho^2  # Variance in X explained by Y's canonical variates
rd_Y_X <- var_Y_V * rho^2  # Variance in Y explained by X's canonical variates

cat("\n=== REDUNDANCY ANALYSIS ===\n")
cat("Redundancy = (own variance extracted) × (canonical r²)\n\n")
cat("Total variance in Socioeconomic (X) explained by Health (Y):", round(sum(rd_X_Y) * 100, 1), "%\n")
cat("Total variance in Health (Y) explained by Socioeconomic (X):", round(sum(rd_Y_X) * 100, 1), "%\n")

redundancy <- data.frame(CV = paste0("CV", 1:n_cc), r = round(rho, 3), r_sq = round(rho^2, 3),
                         Rd_X_Y = round(rd_X_Y, 3), Rd_Y_X = round(rd_Y_X, 3))
print(redundancy)
write.csv(redundancy, "figures/cca_redundancy.csv", row.names = FALSE)

# Plots
cv_data <- data.frame(Country = df_cca$Country, Status = df_cca$Status,
                      U1 = U[,1], V1 = V[,1], U2 = U[,2], V2 = V[,2])

p_cv1 <- ggplot(cv_data, aes(x = U1, y = V1, color = Status)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = paste("CV1 (r =", round(rho[1], 3), ")"), x = "U1 (Socio)", y = "V1 (Health)") +
  theme_minimal()
ggsave("figures/cca_cv1_plot.png", p_cv1, width = 10, height = 8, dpi = 150)

p_cv2 <- ggplot(cv_data, aes(x = U2, y = V2, color = Status)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = paste("CV2 (r =", round(rho[2], 3), ")"), x = "U2", y = "V2") +
  theme_minimal()
ggsave("figures/cca_cv2_plot.png", p_cv2, width = 10, height = 8, dpi = 150)

# Loading plot
load_X <- data.frame(Variable = rownames(loadings_X), CV1 = loadings_X[,1], CV2 = loadings_X[,2], Set = "Socio")
load_Y <- data.frame(Variable = rownames(loadings_Y), CV1 = loadings_Y[,1], CV2 = loadings_Y[,2], Set = "Health")
loadings_all <- bind_rows(load_X, load_Y)

p_load <- ggplot(loadings_all, aes(x = CV1, y = CV2, color = Set, label = Variable)) +
  geom_point(size = 3) + geom_text(hjust = -0.1, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Socio" = "#2E86AB", "Health" = "#E94F37")) +
  labs(title = "Canonical Loadings") + theme_minimal() +
  coord_fixed() + xlim(-1.1, 1.1) + ylim(-1.1, 1.1)
ggsave("figures/cca_loadings_plot.png", p_load, width = 10, height = 10, dpi = 150)

p_bar <- ggplot(cc_df, aes(x = CV, y = r)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(r, 3)), vjust = -0.5) +
  ylim(0, 1) + labs(title = "Canonical Correlations") + theme_minimal()
ggsave("figures/cca_correlations_bar.png", p_bar, width = 8, height = 6, dpi = 150)

cat("\n=== CCA Complete ===\n")
cat("All 8 methods implemented!\n")

# STAT 467 - Factor Analysis
# Input: data_country_level.csv

library(tidyverse)
library(psych)
library(GPArotation)
library(corrplot)
library(gridExtra)
library(nFactors)

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n\n")

fa_vars <- c("Adult_Mortality", "infant_deaths", "Alcohol", "Hepatitis_B", "BMI",
             "under_five_deaths", "Polio", "Total_expenditure", "Diphtheria", "HIV_AIDS",
             "thinness_10_19_years", "thinness_5_9_years", "Income_composition", "Schooling", "GDP_log")

X_fa <- df %>% select(all_of(fa_vars)) %>% drop_na()
df_complete <- df %>% select(Country, Status, all_of(fa_vars)) %>% drop_na()

cat("Variables:", length(fa_vars), ", Obs:", nrow(X_fa), "\n\n")

R <- cor(X_fa)

# KMO
kmo <- KMO(R)
cat("=== KMO ===\n")
cat("Overall MSA:", round(kmo$MSA, 3), "\n")

# Bartlett
bartlett <- cortest.bartlett(R, n = nrow(X_fa))
cat("\n=== Bartlett's Test ===\n")
cat("Chi-sq:", round(bartlett$chisq, 2), ", df:", bartlett$df, ", p:", format(bartlett$p.value, scientific = TRUE), "\n")

# Eigenvalues
eigenvalues <- eigen(R)$values
n_kaiser <- sum(eigenvalues > 1)
cat("\nKaiser criterion:", n_kaiser, "factors\n")

# Parallel analysis
cat("\n=== Parallel Analysis ===\n")
png("figures/fa_parallel_analysis.png", width = 900, height = 600, res = 120)
pa <- fa.parallel(X_fa, fa = "fa", n.iter = 100, main = "Parallel Analysis")
dev.off()
cat("Parallel analysis suggests:", pa$nfact, "factors\n")

n_factors <- max(pa$nfact, 4)
cat("Using", n_factors, "factors\n")

# Factor Analysis - Varimax
cat("\n=== FA WITH VARIMAX ===\n")
fa_varimax <- fa(X_fa, nfactors = n_factors, rotate = "varimax", fm = "ml", scores = "regression")
print(fa_varimax, cut = 0.3, sort = TRUE)

# Factor Analysis - Promax
cat("\n=== FA WITH PROMAX ===\n")
fa_promax <- fa(X_fa, nfactors = n_factors, rotate = "promax", fm = "ml", scores = "regression")
print(fa_promax, cut = 0.3, sort = TRUE)

# Loadings table
loadings_mat <- fa_varimax$loadings
class(loadings_mat) <- "matrix"
loadings_df <- as.data.frame(loadings_mat)
colnames(loadings_df) <- paste0("F", 1:n_factors)
loadings_df$Variable <- rownames(loadings_df)
loadings_df$Communality <- fa_varimax$communality

cat("\n=== LOADINGS ===\n")
loadings_print <- loadings_df[, c("Variable", paste0("F", 1:n_factors), "Communality")]
loadings_print[, -1] <- round(loadings_print[, -1], 3)
print(loadings_print)
write.csv(loadings_df, "figures/fa_loadings.csv", row.names = FALSE)

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

# Factor interpretation
cat("\n=== FACTOR INTERPRETATION ===\n")
cat("F1: Economic Development (GDP, Schooling, Income)\n")
cat("F2: Healthcare Access (Diphtheria, Polio, Hep B)\n")
cat("F3: Mortality Burden (Adult/infant/under-5 deaths)\n")
cat("F4: Nutritional Status (BMI, thinness, Alcohol)\n")

# Factor correlations (Promax)
if (!is.null(fa_promax$Phi)) {
  cat("\n=== FACTOR CORRELATIONS (Promax) ===\n")
  print(round(fa_promax$Phi, 3))

  png("figures/fa_factor_correlations.png", width = 600, height = 500, res = 120)
  corrplot(fa_promax$Phi, method = "color", type = "lower", addCoef.col = "black",
           tl.col = "black", title = "Factor Correlations", mar = c(0, 0, 2, 0))
  dev.off()
}

# Factor scores
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

ggsave("figures/fa_factor_scores_plot.png", grid.arrange(p_f1f2, p_f1f3, ncol = 2), width = 14, height = 6, dpi = 150)

# Model fit
cat("\n=== MODEL FIT ===\n")
cat("TLI:", round(fa_varimax$TLI, 3), "\n")
cat("RMSEA:", round(fa_varimax$RMSEA[1], 3), "\n")
cat("BIC:", round(fa_varimax$BIC, 2), "\n")

cat("\n=== Factor Analysis Complete ===\n")
cat("Run 06_classification.R next.\n")

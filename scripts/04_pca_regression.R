# STAT 467 - PCA and PC Regression
# Input: data_country_level.csv

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(gridExtra)
library(caret)
library(car)
library(psych)

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n\n")

pca_vars <- c("Adult_Mortality", "infant_deaths", "Alcohol", "Hepatitis_B", "BMI",
              "under_five_deaths", "Polio", "Total_expenditure", "Diphtheria", "HIV_AIDS",
              "thinness_10_19_years", "thinness_5_9_years", "Income_composition",
              "Schooling", "GDP_log", "Population_log")

X_pca <- df %>% select(all_of(pca_vars)) %>% drop_na()
df_complete <- df %>% select(Country, Status, Life_expectancy, all_of(pca_vars)) %>% drop_na()

cat("Variables:", length(pca_vars), "\n")
cat("Complete obs:", nrow(X_pca), "\n\n")

# KMO and Bartlett
cor_matrix <- cor(X_pca)
kmo <- KMO(cor_matrix)
cat("KMO:", round(kmo$MSA, 3), "\n")

bartlett <- cortest.bartlett(cor_matrix, n = nrow(X_pca))
cat("Bartlett p-value:", format(bartlett$p.value, scientific = TRUE), "\n\n")

# PCA
pca_result <- prcomp(X_pca, center = TRUE, scale. = TRUE)

eigenvalues <- pca_result$sdev^2
var_exp <- eigenvalues / sum(eigenvalues)
cum_var <- cumsum(var_exp)

pca_summary <- data.frame(
  PC = paste0("PC", 1:length(eigenvalues)),
  Eigenvalue = round(eigenvalues, 4),
  Var_Pct = round(var_exp * 100, 2),
  Cum_Pct = round(cum_var * 100, 2)
)
cat("=== EIGENVALUES ===\n")
print(pca_summary)
write.csv(pca_summary, "figures/pca_eigenvalues.csv", row.names = FALSE)

n_kaiser <- sum(eigenvalues > 1)
n_70var <- which(cum_var >= 0.70)[1]
n_80var <- which(cum_var >= 0.80)[1]

cat("\nKaiser criterion:", n_kaiser, "components\n")
cat("70% variance:", n_70var, "components\n")
cat("80% variance:", n_80var, "components\n")

n_components <- max(n_kaiser, n_70var)

# Scree plot
p_scree <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50),
                    barfill = "steelblue", linecolor = "red") +
  geom_hline(yintercept = 100/length(eigenvalues), linetype = "dashed") +
  labs(title = "Scree Plot") + theme_minimal()
ggsave("figures/pca_screeplot.png", p_scree, width = 10, height = 6, dpi = 150)

# Cumulative variance plot
p_cumvar <- ggplot(pca_summary, aes(x = 1:nrow(pca_summary), y = Cum_Pct)) +
  geom_line(color = "steelblue", linewidth = 1.2) + geom_point(color = "steelblue", size = 3) +
  geom_hline(yintercept = c(70, 80), linetype = "dashed", color = c("orange", "red")) +
  labs(title = "Cumulative Variance", x = "Components", y = "Cumulative %") + theme_minimal()
ggsave("figures/pca_cumulative_variance.png", p_cumvar, width = 10, height = 6, dpi = 150)

# Loadings
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

# Biplot
p_biplot <- fviz_pca_biplot(pca_result, geom.ind = "point", col.ind = df_complete$Status,
                             palette = c("#2E86AB", "#E94F37"), addEllipses = TRUE,
                             ellipse.level = 0.95, repel = TRUE, legend.title = "Status") +
  labs(title = "PCA Biplot") + theme_minimal()
ggsave("figures/pca_biplot.png", p_biplot, width = 14, height = 10, dpi = 150)

# PC scores plot
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

# Contributions
p_contrib1 <- fviz_contrib(pca_result, choice = "var", axes = 1, top = 10) + theme_minimal()
p_contrib2 <- fviz_contrib(pca_result, choice = "var", axes = 2, top = 10) + theme_minimal()
ggsave("figures/pca_contributions.png", grid.arrange(p_contrib1, p_contrib2, ncol = 2), width = 14, height = 6, dpi = 150)

# PC Regression
cat("\n=== PC REGRESSION ===\n")
pc_data <- as.data.frame(pca_result$x)
pc_data$Life_expectancy <- df_complete$Life_expectancy
pc_data$Status <- df_complete$Status

pcr_m3 <- lm(Life_expectancy ~ PC1 + PC2 + PC3, data = pc_data)
pcr_m5 <- lm(Life_expectancy ~ PC1 + PC2 + PC3 + PC4 + PC5, data = pc_data)

cat("Model (3 PCs): R² =", round(summary(pcr_m3)$r.squared, 4), "\n")
cat("Model (5 PCs): R² =", round(summary(pcr_m5)$r.squared, 4), "\n")

cat("\n=== 5-PC MODEL SUMMARY ===\n")
print(summary(pcr_m5))

write.csv(data.frame(
  Term = names(coef(pcr_m5)),
  Coefficient = round(coef(pcr_m5), 4),
  p_value = summary(pcr_m5)$coefficients[, 4]
), "figures/pcr_coefficients.csv", row.names = FALSE)

# Cross-validation
cat("\n=== CROSS-VALIDATION ===\n")
set.seed(123)
train_control <- trainControl(method = "cv", number = 5)

cv_results <- data.frame(n_PCs = 1:10, RMSE = NA, R2 = NA)
for (k in 1:10) {
  formula_k <- as.formula(paste("Life_expectancy ~", paste(paste0("PC", 1:k), collapse = " + ")))
  model_cv <- train(formula_k, data = pc_data, method = "lm", trControl = train_control)
  cv_results$RMSE[k] <- model_cv$results$RMSE
  cv_results$R2[k] <- model_cv$results$Rsquared
}
print(cv_results)
write.csv(cv_results, "figures/pcr_cv_results.csv", row.names = FALSE)

p_cv <- ggplot(cv_results, aes(x = n_PCs, y = R2)) +
  geom_line(color = "steelblue", linewidth = 1.2) + geom_point(color = "steelblue", size = 3) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "CV R² by Number of PCs", x = "PCs", y = "R²") + theme_minimal()
ggsave("figures/pcr_cv_r2.png", p_cv, width = 10, height = 6, dpi = 150)

cat("Optimal PCs:", cv_results$n_PCs[which.max(cv_results$R2)], "\n")

# Standard regression comparison
std_model <- lm(Life_expectancy ~ ., data = df_complete %>% select(Life_expectancy, all_of(pca_vars)))
cat("\nStandard Regression R²:", round(summary(std_model)$r.squared, 4), "\n")

vif_vals <- car::vif(std_model)
cat("Variables with VIF > 5:", sum(vif_vals > 5), "\n")

# Diagnostics
p_fitted <- ggplot(pc_data, aes(x = fitted(pcr_m5), y = Life_expectancy)) +
  geom_point(aes(color = Status), alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = "Fitted vs Actual (5-PC Model)") + theme_minimal()

p_resid <- ggplot(data.frame(Fitted = fitted(pcr_m5), Residuals = residuals(pcr_m5)),
                  aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residual Plot") + theme_minimal()

ggsave("figures/pcr_diagnostics.png", grid.arrange(p_fitted, p_resid, ncol = 2), width = 14, height = 6, dpi = 150)

cat("\n=== PCA and PC Regression Complete ===\n")
cat("Run 05_factor_analysis.R next.\n")

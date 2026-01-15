# STAT 467 - Exploratory Data Analysis
# Input: data_country_level.csv

library(tidyverse)
library(corrplot)
library(ggcorrplot)
library(GGally)
library(gridExtra)
library(MVN)

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries,", ncol(df), "variables\n")
cat("Status:", table(df$Status), "\n\n")

key_vars <- c("Life_expectancy", "Adult_Mortality", "infant_deaths", "Alcohol",
              "Hepatitis_B", "BMI", "under_five_deaths", "Polio", "Total_expenditure",
              "Diphtheria", "HIV_AIDS", "thinness_10_19_years", "thinness_5_9_years",
              "Income_composition", "Schooling", "GDP_log", "Population_log")

# Descriptive statistics
cat("=== Descriptive Statistics ===\n")
print(summary(df[, key_vars]))

# Summary by Status
summary_by_status <- df %>%
  group_by(Status) %>%
  summarise(n = n(),
            Life_exp = mean(Life_expectancy),
            Mortality = mean(Adult_Mortality),
            Schooling = mean(Schooling),
            GDP = mean(GDP_log))
print(summary_by_status)

# Correlation matrix
cor_matrix <- cor(df[, key_vars], use = "pairwise.complete.obs")
write.csv(cor_matrix, "figures/correlation_matrix.csv")

png("figures/correlation_heatmap.png", width = 1400, height = 1200, res = 120)
corrplot(cor_matrix, method = "color", type = "lower", order = "hclust",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.5)
dev.off()

# High correlations
high_cor <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("\nHigh correlations (|r| > 0.7):\n")
  for (i in 1:nrow(high_cor)) {
    if (high_cor[i,1] < high_cor[i,2]) {
      cat(rownames(cor_matrix)[high_cor[i,1]], "-",
          colnames(cor_matrix)[high_cor[i,2]], ":",
          round(cor_matrix[high_cor[i,1], high_cor[i,2]], 3), "\n")
    }
  }
}

# Boxplots by Status
selected_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
                   "Income_composition", "HIV_AIDS", "BMI", "Alcohol")

p_list <- lapply(selected_vars, function(var) {
  ggplot(df, aes_string(x = "Status", y = var, fill = "Status")) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
    labs(title = var) +
    theme_minimal() +
    theme(legend.position = "none")
})

p_boxplots <- do.call(grid.arrange, c(p_list, ncol = 4))
ggsave("figures/boxplots_by_status.png", p_boxplots, width = 14, height = 8, dpi = 150)

# Scatter plots - Life expectancy relationships
p1 <- ggplot(df, aes(x = Schooling, y = Life_expectancy, color = Status)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()

p2 <- ggplot(df, aes(x = GDP_log, y = Life_expectancy, color = Status)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()

p3 <- ggplot(df, aes(x = Adult_Mortality, y = Life_expectancy, color = Status)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()

p4 <- ggplot(df, aes(x = Income_composition, y = Life_expectancy, color = Status)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()

ggsave("figures/scatter_life_expectancy.png", grid.arrange(p1, p2, p3, p4, ncol = 2),
       width = 12, height = 10, dpi = 150)

# Pairs plot
pairs_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log", "Income_composition", "HIV_AIDS")
p_pairs <- ggpairs(df, columns = pairs_vars, mapping = aes(color = Status, alpha = 0.5),
                   upper = list(continuous = wrap("cor", size = 3)),
                   lower = list(continuous = wrap("points", alpha = 0.5, size = 1))) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  theme_minimal()
ggsave("figures/pairs_plot.png", p_pairs, width = 14, height = 12, dpi = 150)

# Multivariate normality check
norm_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log", "Income_composition", "BMI", "HIV_AIDS")
mvn_result <- MVN::mvn(df[, norm_vars], mvn_test = "mardia")
cat("\n=== Multivariate Normality (Mardia) ===\n")
print(mvn_result$multivariateNormality)

# Q-Q plot
mahal_dist <- mahalanobis(df[, norm_vars], colMeans(df[, norm_vars]), cov(df[, norm_vars]))
png("figures/mvn_qq_plot.png", width = 800, height = 600, res = 120)
qqplot(qchisq(ppoints(nrow(df)), df = length(norm_vars)), mahal_dist,
       main = "Q-Q Plot: Mahalanobis Distances", xlab = "Chi-squared Quantiles", ylab = "Mahalanobis Distances")
abline(0, 1, col = "red", lwd = 2)
dev.off()

cat("\n=== EDA Complete ===\n")
cat("Figures saved. Run 02_mean_vector_inference.R next.\n")

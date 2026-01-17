# STAT 467 - Data Preprocessing
# Input: data.csv | Output: data_cleaned.csv, data_country_level.csv

library(tidyverse)
library(mice)
library(VIM)

# Load data (use population-corrected version)
df_raw <- read.csv("data_population_corrected.csv", stringsAsFactors = FALSE,
                   check.names = FALSE)

# Clean column names
colnames(df_raw) <- trimws(colnames(df_raw))
colnames(df_raw) <- gsub("\\s+", "_", colnames(df_raw))

df_raw <- df_raw %>%
  rename(
    HIV_AIDS = `HIV/AIDS`,
    thinness_10_19_years = `thinness_1-19_years`,
    thinness_5_9_years = `thinness_5-9_years`,
    under_five_deaths = `under-five_deaths`,
    Income_composition = Income_composition_of_resources
  )

cat("Raw data:", dim(df_raw)[1], "x", dim(df_raw)[2], "\n")

# Missing value analysis
missing_summary <- df_raw %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing") %>%
  mutate(Percent = round(Missing / nrow(df_raw) * 100, 2)) %>%
  filter(Missing > 0) %>%
  arrange(desc(Percent))

print(missing_summary)

# Visualize missing data
png("figures/missing_data_pattern.png", width = 1200, height = 800)
VIM::aggr(df_raw, col = c('navyblue', 'red'), numbers = TRUE, sortVars = TRUE,
          cex.axis = 0.7, gap = 3)
dev.off()

# Remove rows missing Life expectancy
df_clean <- df_raw %>% filter(!is.na(Life_expectancy))
cat("After removing missing Life_expectancy:", nrow(df_clean), "rows\n")

# MICE imputation for numeric columns
# NOTE: Multiple Imputation by Chained Equations (MICE)
# - m=5: Creates 5 imputed datasets
# - maxit=10: 10 iterations for convergence
# - method='pmm': Predictive Mean Matching (preserves distribution)
#
# IMPORTANT: We use only the FIRST imputed dataset (complete(mice_model, 1))
# Proper multiple imputation would require:
#   1. Running analyses on all m=5 datasets
#   2. Combining results using Rubin's rules
# For this educational project, single imputation is used for simplicity.
# This may underestimate standard errors - interpret CIs with caution.

numeric_cols <- df_clean %>%
  select(where(is.numeric), -Year) %>%
  colnames()

set.seed(123)
df_for_imputation <- df_clean[, numeric_cols]

cat("\n=== MICE IMPUTATION ===\n")
cat("Method: Predictive Mean Matching (pmm)\n")
cat("Number of imputations (m):", 5, "\n")
cat("Max iterations:", 10, "\n")

mice_model <- mice(df_for_imputation, m = 5, maxit = 10, method = 'pmm',
                   seed = 123, printFlag = FALSE)

# Using first imputed dataset (see note above about Rubin's rules)
df_imputed <- complete(mice_model, 1)
df_clean[, numeric_cols] <- df_imputed

cat("NOTE: Using single imputed dataset (1 of 5) for simplicity\n")
cat("Imputation complete. Remaining NA:", sum(is.na(df_clean[, numeric_cols])), "\n")

# Log transformations for skewed variables
df_clean <- df_clean %>%
  mutate(
    GDP_log = log1p(GDP),
    Population_log = log1p(Population),
    Measles_log = log1p(Measles),
    percentage_expenditure_log = log1p(percentage_expenditure),
    infant_deaths_log = log1p(infant_deaths),
    under_five_deaths_log = log1p(under_five_deaths)
  )

df_clean$Status <- as.factor(df_clean$Status)

# Aggregate to country level
numeric_vars <- df_clean %>% select(where(is.numeric), -Year) %>% colnames()

df_country <- df_clean %>%
  group_by(Country) %>%
  summarise(
    Status = first(Status),
    across(all_of(numeric_vars), ~mean(., na.rm = TRUE)),
    n_years = n()
  ) %>%
  ungroup()

cat("Country-level data:", nrow(df_country), "countries\n")
print(table(df_country$Status))

# Outlier detection (Mahalanobis)
key_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log", "Income_composition", "BMI", "HIV_AIDS")
X_outlier <- df_country[, key_vars] %>% as.matrix()
mahal_dist <- mahalanobis(X_outlier, colMeans(X_outlier), cov(X_outlier))
df_country$Mahalanobis_Distance <- mahal_dist

chi_sq_crit <- qchisq(0.999, df = length(key_vars))
outliers <- df_country %>% filter(Mahalanobis_Distance > chi_sq_crit)
cat("Potential outliers:", nrow(outliers), "\n")

# Save datasets
write.csv(df_clean, "data_cleaned.csv", row.names = FALSE)
write.csv(df_country, "data_country_level.csv", row.names = FALSE)

cat("\nSaved: data_cleaned.csv (", nrow(df_clean), " rows)\n")
cat("Saved: data_country_level.csv (", nrow(df_country), " rows)\n")
cat("\nPreprocessing complete. Run 01_eda.R next.\n")

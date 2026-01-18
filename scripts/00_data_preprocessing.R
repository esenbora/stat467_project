# ============================================================================
# STAT 467 - DATA PREPROCESSING
# ============================================================================
# Purpose: Clean and prepare WHO Life Expectancy data for multivariate analysis
#
# Input:   data_population_corrected.csv (raw data with fixed population values)
#          - Original data had implausible population values fixed via World Bank API
#          - See 00a_population_fix.R for details on population correction
#
# Output:  data_cleaned.csv (2938 rows, year-level observations after imputation)
#          data_country_level.csv (193 rows, country means across 2000-2015)
#
# Methods:
#   - MICE imputation (Predictive Mean Matching)
#   - Log transformation for skewed variables
#   - Aggregation to country-level means
#   - Mahalanobis distance for outlier detection
#
# Dependencies: tidyverse, mice, VIM
# ============================================================================

library(tidyverse)
library(mice)
library(VIM)

# ----------------------------------------------------------------------------
# STEP 1: Load and Clean Raw Data
# ----------------------------------------------------------------------------
# Load the population-corrected dataset
# Note: Population values were corrected using World Bank API (see 00a_population_fix.R)
df_raw <- read.csv("data_population_corrected.csv", stringsAsFactors = FALSE,
                   check.names = FALSE)

# Standardize column names: trim whitespace and replace spaces with underscores
colnames(df_raw) <- trimws(colnames(df_raw))
colnames(df_raw) <- gsub("\\s+", "_", colnames(df_raw))

# Rename problematic column names for easier R handling
# - Replace slashes and hyphens with underscores
# - Shorten overly long names
df_raw <- df_raw %>%
  rename(
    HIV_AIDS = `HIV/AIDS`,
    thinness_10_19_years = `thinness_1-19_years`,
    thinness_5_9_years = `thinness_5-9_years`,
    under_five_deaths = `under-five_deaths`,
    Income_composition = Income_composition_of_resources
  )

cat("Raw data:", dim(df_raw)[1], "x", dim(df_raw)[2], "\n")

# ----------------------------------------------------------------------------
# STEP 2: Missing Value Analysis
# ----------------------------------------------------------------------------
# Calculate missing values per variable
# This helps identify which variables need imputation and their severity
missing_summary <- df_raw %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing") %>%
  mutate(Percent = round(Missing / nrow(df_raw) * 100, 2)) %>%
  filter(Missing > 0) %>%
  arrange(desc(Percent))

print(missing_summary)

# Visualize missing data patterns using VIM package
# Blue = present, Red = missing
# This helps identify if missingness is random (MCAR) or patterned
png("figures/missing_data_pattern.png", width = 1200, height = 800)
VIM::aggr(df_raw, col = c('navyblue', 'red'), numbers = TRUE, sortVars = TRUE,
          cex.axis = 0.7, gap = 3)
dev.off()

# ----------------------------------------------------------------------------
# STEP 3: Handle Missing Life Expectancy (Response Variable)
# ----------------------------------------------------------------------------
# Remove rows where Life_expectancy is missing
# Rationale: This is our primary response variable - imputing it would be circular
df_clean <- df_raw %>% filter(!is.na(Life_expectancy))
cat("After removing missing Life_expectancy:", nrow(df_clean), "rows\n")

# ----------------------------------------------------------------------------
# STEP 4: Multiple Imputation by Chained Equations (MICE)
# ----------------------------------------------------------------------------
# MICE Algorithm Overview:
#   1. Fill missing values with initial guesses (e.g., column means)
#   2. For each variable with missing values:
#      a. Set missing values back to NA
#      b. Regress this variable on all other variables
#      c. Replace missing values with predictions from the model
#   3. Repeat for multiple iterations (maxit) until convergence
#   4. Create m imputed datasets
#
# Parameters:
#   - m=5: Creates 5 imputed datasets (standard in literature)
#   - maxit=10: 10 iterations for algorithm convergence
#   - method='pmm': Predictive Mean Matching
#
# Why Predictive Mean Matching (PMM)?
#   - Preserves the original distribution of observed values
#   - Imputed values are always plausible (drawn from observed data)
#   - Robust to model misspecification
#   - Better for non-normal data than regression imputation
#
# IMPORTANT STATISTICAL NOTE:
#   We use only the FIRST imputed dataset (complete(mice_model, 1)).
#   Proper multiple imputation would require:
#     1. Running all analyses on all m=5 datasets
#     2. Combining results using Rubin's rules (pooling estimates and SE)
#   For this educational project, single imputation is used for simplicity.
#   Consequence: Standard errors may be underestimated - interpret CIs with caution.

numeric_cols <- df_clean %>%
  select(where(is.numeric), -Year) %>%
  colnames()

set.seed(123)  # For reproducibility
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

# ----------------------------------------------------------------------------
# STEP 5: Log Transformations for Skewed Variables
# ----------------------------------------------------------------------------
# Rationale: Many health/economic variables are highly right-skewed
# Log transformation helps:
#   - Normalize distributions for parametric tests (MVN assumption)
#   - Stabilize variance (homoscedasticity)
#   - Make relationships more linear
#
# We use log1p(x) = log(1 + x) to handle zero values
# Variables transformed: GDP, Population, Measles, expenditure, deaths
df_clean <- df_clean %>%
  mutate(
    GDP_log = log1p(GDP),
    Population_log = log1p(Population),
    Measles_log = log1p(Measles),
    percentage_expenditure_log = log1p(percentage_expenditure),
    infant_deaths_log = log1p(infant_deaths),
    under_five_deaths_log = log1p(under_five_deaths)
  )

# Convert Status to factor for downstream analyses
df_clean$Status <- as.factor(df_clean$Status)

# ----------------------------------------------------------------------------
# STEP 6: Aggregate to Country-Level Data
# ----------------------------------------------------------------------------
# Rationale for aggregation:
#   - Original data has 16 years per country (2000-2015) = repeated measures
#   - Multivariate methods assume independent observations
#   - Country-level means provide stable estimates and satisfy independence
#   - n=193 countries provides adequate sample size for multivariate analysis
#
# Aggregation method: Mean across all years for each country
numeric_vars <- df_clean %>% select(where(is.numeric), -Year) %>% colnames()

df_country <- df_clean %>%
  group_by(Country) %>%
  summarise(
    Status = first(Status),  # Status is constant within country
    across(all_of(numeric_vars), ~mean(., na.rm = TRUE)),
    n_years = n()  # Track how many years of data per country
  ) %>%
  ungroup()

cat("Country-level data:", nrow(df_country), "countries\n")
print(table(df_country$Status))

# ----------------------------------------------------------------------------
# STEP 7: Multivariate Outlier Detection (Mahalanobis Distance)
# ----------------------------------------------------------------------------
# Mahalanobis distance measures how far each observation is from the centroid,
# accounting for correlations between variables.
#
# Formula: D² = (x - μ)' Σ⁻¹ (x - μ)
#
# Under multivariate normality, D² follows a chi-squared distribution with p df.
# We use alpha = 0.001 (very conservative) to flag extreme outliers.
# These are flagged but NOT removed - removal requires domain justification.
key_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log", "Income_composition", "BMI", "HIV_AIDS")
X_outlier <- df_country[, key_vars] %>% as.matrix()
mahal_dist <- mahalanobis(X_outlier, colMeans(X_outlier), cov(X_outlier))
df_country$Mahalanobis_Distance <- mahal_dist

# Chi-squared critical value at alpha = 0.001
chi_sq_crit <- qchisq(0.999, df = length(key_vars))
outliers <- df_country %>% filter(Mahalanobis_Distance > chi_sq_crit)
cat("Potential outliers:", nrow(outliers), "\n")

# ----------------------------------------------------------------------------
# STEP 8: Save Processed Datasets
# ----------------------------------------------------------------------------
write.csv(df_clean, "data_cleaned.csv", row.names = FALSE)
write.csv(df_country, "data_country_level.csv", row.names = FALSE)

cat("\nSaved: data_cleaned.csv (", nrow(df_clean), " rows)\n")
cat("Saved: data_country_level.csv (", nrow(df_country), " rows)\n")
cat("\nPreprocessing complete. Run 01_eda.R next.\n")

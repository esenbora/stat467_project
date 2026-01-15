# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

STAT 467 Multivariate Data Analysis term project analyzing the WHO Life Expectancy Dataset (2000-2015). Implements 8 multivariate statistical methods to examine relationships between socioeconomic/health factors across 193 countries, with focus on Developed vs. Developing country comparisons.

## Running the Analysis

Scripts must be executed sequentially in R/RStudio (from the project root directory):

```r
source("scripts/utils.R")                      # Load utility functions first
source("scripts/00a_population_fix.R")         # Fix Population data (World Bank API) - run once
source("scripts/00_data_preprocessing.R")      # Creates data_cleaned.csv, data_country_level.csv
source("scripts/01_eda.R")
source("scripts/02_mean_vector_inference.R")
source("scripts/03_manova.R")
source("scripts/04_pca_regression.R")
source("scripts/05_factor_analysis.R")
source("scripts/06_classification.R")
source("scripts/07_clustering.R")
source("scripts/08_cca.R")
```

To re-run a specific analysis, you can source from that script onward (earlier scripts create required datasets).

## Architecture

**Data Flow Pipeline:**
```
data.csv → 00a_population_fix.R → data_population_corrected.csv
                                              ↓
                         00_data_preprocessing.R → data_cleaned.csv + data_country_level.csv
                                              ↓
                              01-08 analysis scripts → figures/*.png + figures/*.csv
```

**Key Design Decisions:**
- Raw data (2,938 time-series rows) is aggregated to country-level means (193 rows) for multivariate analysis
- MICE imputation (m=5, maxit=10, PMM) handles missing values before aggregation
- Log transforms applied to highly skewed variables (GDP, Population, Measles, expenditure, deaths)
- All scripts use `set.seed(123)` for reproducibility

## Utils Module (scripts/utils.R)

Reusable functions available after sourcing `utils.R`:

| Function | Purpose |
|----------|---------|
| `calc_mahalanobis(X)` | Mahalanobis distance from centroid |
| `identify_outliers(X, alpha)` | Chi-squared based outlier detection |
| `hotelling_one_sample(X, mu0)` | One-sample Hotelling's T² test |
| `hotelling_two_sample(X1, X2)` | Two-sample Hotelling's T² test |
| `bonferroni_ci(X, alpha)` | Confidence intervals with Bonferroni correction |
| `create_cor_heatmap(cor_mat, ...)` | Correlation matrix visualization |
| `create_mvn_qq(X, ...)` | MVN Q-Q plot with chi-squared quantiles |
| `find_optimal_k(X, max_k)` | WSS and silhouette for optimal clustering k |
| `cluster_profile(data, cluster_var, vars)` | Summary statistics by cluster |
| `compare_classifiers(actual, predictions)` | Model comparison (accuracy, kappa) |
| `standardize(X)` / `unstandardize(X_std)` | Scaling with preserved attributes |
| `load_packages(packages)` | Auto-install and load dependencies |
| `print_section(title)` / `print_subsection(title)` | Console output formatting |

## Key Variable Groups

Used consistently across analysis scripts:

- **Response:** `Life.expectancy`
- **Health indicators:** `Adult.Mortality`, `infant.deaths`, `HIV.AIDS`, `Hepatitis.B`, `Polio`, `Diphtheria`, `Measles`, `BMI`
- **Socioeconomic:** `GDP`, `Schooling`, `Income.composition.of.resources`, `percentage.expenditure`, `Total.expenditure`
- **Grouping factor:** `Status` (Developed/Developing)

## Required R Packages

Core packages (auto-installed via `load_packages()`):
- tidyverse, dplyr, ggplot2, corrplot, gridExtra, factoextra
- MASS, cluster, psych, car, MVN, caret
- mice, Hotelling, GPArotation, CCA, biotools, dendextend, mclust, pROC, nFactors
- WDI (World Bank data API for population correction)

## Output

All figures and result tables are saved to `figures/` directory as PNG and CSV files.

# STAT 467 - Multivariate Data Analysis Term Project

## Dataset
WHO Life Expectancy Dataset (2000-2015)
- 193 countries, 22 variables
- Aggregated to country-level means

## Methods Implemented
1. **Data Preprocessing** - MICE imputation, log transforms
2. **EDA** - Correlation analysis, MVN testing
3. **Mean Vector Inference** - Hotelling's T² (one/two-sample)
4. **MANOVA** - Pillai, Wilks, Hotelling-Lawley, Roy
5. **PCA & Regression** - Principal component regression
6. **Factor Analysis** - Varimax/Promax rotation
7. **Classification** - LDA, QDA with ROC analysis
8. **Clustering** - Hierarchical, K-means
9. **CCA** - Canonical correlation analysis

## Usage
```r
# Run scripts in order
setwd("scripts")
source("utils.R")
source("00_data_preprocessing.R")  # Creates cleaned data
source("01_eda.R")
# ... continue through 08_cca.R
```

## Structure
```
├── data.csv                 # Original dataset
├── scripts/
│   ├── utils.R              # Helper functions
│   ├── 00_data_preprocessing.R
│   ├── 01_eda.R
│   ├── 02_mean_vector_inference.R
│   ├── 03_manova.R
│   ├── 04_pca_regression.R
│   ├── 05_factor_analysis.R
│   ├── 06_classification.R
│   ├── 07_clustering.R
│   └── 08_cca.R
└── figures/                 # Output directory
```

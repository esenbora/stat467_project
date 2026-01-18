# ============================================================================
# STAT 467 - CLASSIFICATION (LDA/QDA)
# ============================================================================
# Purpose: Classify countries as Developed/Developing using discriminant analysis
#          Compare Linear (LDA) and Quadratic (QDA) discriminant functions
#
# Input:   data_country_level.csv (country-level aggregated data)
#
# Output:  - ROC curves (figures/classification_roc_curves.png)
#          - Confusion matrices (figures/classification_confusion_matrices.png)
#          - LDA projection (figures/lda_projection.png)
#          - Partition plots (figures/classification_partition_*.png)
#          - Variable importance (figures/classification_var_importance.png)
#          - LDA coefficients (figures/lda_coefficients.csv)
#          - Misclassified countries (figures/lda_misclassified.csv)
#
# Methods:
#   1. Linear Discriminant Analysis (LDA)
#   2. Quadratic Discriminant Analysis (QDA)
#   3. 5-fold cross-validation
#   4. ROC/AUC evaluation
#
# Key Concepts:
#   LDA vs QDA:
#     - LDA assumes equal covariance matrices (Σ₁ = Σ₂)
#       → Uses pooled covariance, produces linear decision boundaries
#       → More stable with small samples, fewer parameters to estimate
#     - QDA allows different covariances per group
#       → Produces quadratic (curved) decision boundaries
#       → Requires more data per group (estimates pq(p+1)/2 parameters per group)
#
#   Classification Rule (Bayes' Theorem):
#     Assign x to group k that maximizes: P(k|x) = P(x|k) × P(k) / P(x)
#     Where P(k) is the prior probability and P(x|k) is the likelihood
#
#   Assumptions:
#     1. Multivariate normality within each group
#     2. Equal covariances across groups (LDA only)
#     3. Independence of observations
#
# Prior Probabilities:
#   - Equal priors (0.5, 0.5): Assume equal group sizes
#   - Observed priors: Use sample proportions (recommended for imbalanced data)
#   - Affects decision boundary location
#
# Dependencies: tidyverse, MASS, caret, pROC, gridExtra, biotools, MVN, klaR
# ============================================================================

library(tidyverse)
library(MASS)
library(caret)
library(pROC)
library(gridExtra)
library(biotools)  # For Box's M test
library(MVN)       # For multivariate normality
library(klaR)      # For partition plots (partimat)

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

df <- read.csv("data_country_level.csv", stringsAsFactors = FALSE)
df$Status <- as.factor(df$Status)

cat("Data:", nrow(df), "countries\n")
print(table(df$Status))

# Select predictor variables
# These capture health outcomes, socioeconomic status, and development indicators
# Note: HIV_AIDS excluded as it is constant (0.1) for all Developed countries,
# causing singular covariance matrix
pred_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
               "Income_composition", "BMI", "Diphtheria", "thinness_10_19_years")

df_class <- df %>% select(Country, Status, all_of(pred_vars)) %>% drop_na()
cat("Complete cases:", nrow(df_class), "\n\n")

# ============================================================================
# SECTION 1: TRAIN-TEST SPLIT
# ============================================================================
# Stratified split maintains class proportions in both sets
# 70% training / 30% testing is a common split for moderate sample sizes

set.seed(123)
train_idx <- createDataPartition(df_class$Status, p = 0.7, list = FALSE)
train_data <- df_class[train_idx, ]
test_data <- df_class[-train_idx, ]

cat("Train:", nrow(train_data), ", Test:", nrow(test_data), "\n\n")

# ============================================================================
# SECTION 2: LDA ASSUMPTION TESTS
# ============================================================================
# LDA requires:
#   1. Multivariate normality within each group
#   2. Homogeneous covariance matrices (Σ₁ = Σ₂)
# QDA only requires assumption 1 (allows heterogeneous covariances)

cat("=== LDA ASSUMPTION TESTS ===\n")

# ----------------------------------------------------------------------------
# 2.1 Box's M Test for Equality of Covariance Matrices
# ----------------------------------------------------------------------------
# Tests H₀: Σ₁ = Σ₂ = ... = Σₖ
# If rejected (p < 0.05): Covariances differ, LDA assumption violated
# Note: Box's M is sensitive to departures from MVN
#       With large samples, even minor differences can be significant

cat("\n1. Box's M Test (Covariance Homogeneity):\n")
X_class_matrix <- as.matrix(df_class[, pred_vars])
box_m_result <- biotools::boxM(X_class_matrix, df_class$Status)
print(box_m_result)

if (is.infinite(box_m_result$statistic) || is.na(box_m_result$p.value)) {
  cat("WARNING: Box's M test returned Inf/NA - possible singularity\n")
  cat("Recommendation: LDA may still be appropriate, but interpret with caution\n")
} else if (box_m_result$p.value < 0.05) {
  cat("WARNING: Covariance matrices are NOT equal (p < 0.05)\n")
  cat("Recommendation: Consider QDA or interpret LDA results cautiously\n")
} else {
  cat("Covariance homogeneity assumption is met (p >= 0.05)\n")
}

# ----------------------------------------------------------------------------
# 2.2 Multivariate Normality by Group
# ----------------------------------------------------------------------------
# MVN is required for both LDA and QDA
# Using Mardia's test which evaluates multivariate skewness and kurtosis
# Note: With small group sizes, MVN tests have low power

cat("\n2. Multivariate Normality by Group (Mardia's Test):\n")
for (grp in levels(df_class$Status)) {
  cat("\n  ", grp, ":\n", sep = "")
  grp_data <- df_class %>% dplyr::filter(Status == grp) %>%
    dplyr::select(all_of(pred_vars))
  tryCatch({
    mvn_result <- MVN::mvn(grp_data, mvn_test = "mardia")
    cat("    Skewness p =", format(mvn_result$multivariate_normality$p.value[1],
                                   digits = 3), "\n")
    cat("    Kurtosis p =", format(mvn_result$multivariate_normality$p.value[2],
                                   digits = 3), "\n")
    # New MVN package uses "✓ Normal" / "✗ Not Normal" instead of "YES"/"NO"
    if (!grepl("Normal", mvn_result$multivariate_normality$MVN[1]) ||
        !grepl("Normal", mvn_result$multivariate_normality$MVN[2])) {
      cat("    WARNING: MVN may be violated for this group\n")
    } else {
      cat("    MVN assumption met\n")
    }
  }, error = function(e) {
    cat("    MVN test failed:", conditionMessage(e), "\n")
  })
}

# ----------------------------------------------------------------------------
# 2.3 Prior Probabilities
# ----------------------------------------------------------------------------
# Priors affect the decision boundary:
#   - Higher prior for a class shifts boundary toward that class
#   - Using observed proportions is generally recommended for imbalanced data
#   - Equal priors (0.5, 0.5) assume no prior knowledge of class membership

class_props <- prop.table(table(train_data$Status))
cat("\n3. Class Proportions (for prior probabilities):\n")
print(round(class_props, 3))

# ============================================================================
# SECTION 3: LINEAR DISCRIMINANT ANALYSIS (LDA)
# ============================================================================
# LDA finds linear combinations of predictors that maximize between-class
# to within-class variance ratio (Fisher's criterion)
#
# Discriminant function: LD₁ = w₁x₁ + w₂x₂ + ... + wₚxₚ
# Where w = Sp⁻¹(x̄₁ - x̄₂) and Sp is the pooled covariance matrix
#
# Classification: Assign to group with higher posterior probability

cat("\n=== LDA ===\n")
# LDA uses training class proportions as priors by default
# (more realistic than equal priors)
lda_model <- lda(Status ~ ., data = train_data[, c("Status", pred_vars)])
print(lda_model)

# Predict on test set
lda_pred <- predict(lda_model, test_data)

# Confusion matrix with various performance metrics
lda_cm <- confusionMatrix(lda_pred$class, test_data$Status)
print(lda_cm)

# Save discriminant coefficients (LD1)
# These show each variable's contribution to the linear discriminant function
# Larger |coefficient| = more discriminating power (after standardization)
lda_coef <- data.frame(Variable = rownames(lda_model$scaling), LD1 = round(lda_model$scaling[,1], 4))
lda_coef <- lda_coef %>% arrange(desc(abs(LD1)))
write.csv(lda_coef, "figures/lda_coefficients.csv", row.names = FALSE)

# ============================================================================
# SECTION 4: QUADRATIC DISCRIMINANT ANALYSIS (QDA)
# ============================================================================
# QDA estimates separate covariance matrices for each group
# This allows quadratic (curved) decision boundaries
#
# Advantage: More flexible when group covariances truly differ
# Disadvantage: Requires more parameters → needs larger sample per group
#              With p predictors and k groups: estimates k × p(p+1)/2 parameters
#
# May fail with rank deficiency when n_group < p (too few observations per group)

cat("\n=== QDA ===\n")
qda_success <- FALSE
tryCatch({
  # QDA uses training class proportions as priors by default
  qda_model <- qda(Status ~ ., data = train_data[, c("Status", pred_vars)])
  qda_pred <- predict(qda_model, test_data)
  qda_cm <- confusionMatrix(qda_pred$class, test_data$Status)
  print(qda_cm)
  qda_success <- TRUE
}, error = function(e) {
  cat("QDA failed (rank deficiency - group size too small for # predictors)\n")
  cat("Skipping QDA, using LDA only\n")
})

# ============================================================================
# SECTION 5: CROSS-VALIDATION
# ============================================================================
# K-fold CV provides unbiased estimate of generalization performance
# Using ROC (area under curve) as metric for imbalanced binary classification
# ROC is preferred over accuracy when class sizes are unequal

cat("\n=== CROSS-VALIDATION ===\n")
cv_ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

lda_cv <- train(Status ~ ., data = df_class[, c("Status", pred_vars)], method = "lda", trControl = cv_ctrl, metric = "ROC")
cat("LDA CV ROC:", round(lda_cv$results$ROC, 4), "\n")

if (qda_success) {
  qda_cv <- train(Status ~ ., data = df_class[, c("Status", pred_vars)], method = "qda", trControl = cv_ctrl, metric = "ROC")
  cat("QDA CV ROC:", round(qda_cv$results$ROC, 4), "\n")
}

# ============================================================================
# SECTION 6: ROC CURVES AND AUC
# ============================================================================
# ROC (Receiver Operating Characteristic) curve plots:
#   - Sensitivity (True Positive Rate) on y-axis
#   - 1 - Specificity (False Positive Rate) on x-axis
#
# AUC (Area Under Curve) interpretation:
#   - 1.0: Perfect discrimination
#   - 0.9-1.0: Excellent
#   - 0.8-0.9: Good
#   - 0.7-0.8: Fair
#   - 0.5: No discrimination (random guessing)
#   - < 0.5: Worse than random (model inverted)

roc_lda <- roc(test_data$Status, lda_pred$posterior[, "Developed"], levels = c("Developing", "Developed"))
cat("\nLDA AUC:", round(auc(roc_lda), 4), "\n")

png("figures/classification_roc_curves.png", width = 900, height = 700, res = 120)
if (qda_success) {
  roc_qda <- roc(test_data$Status, qda_pred$posterior[, "Developed"], levels = c("Developing", "Developed"))
  cat("QDA AUC:", round(auc(roc_qda), 4), "\n")
  plot(roc_lda, col = "#2E86AB", lwd = 2, main = "ROC Curves: LDA vs QDA")
  plot(roc_qda, col = "#E94F37", lwd = 2, add = TRUE)
  abline(a = 0, b = 1, lty = 2, col = "gray")
  legend("bottomright", c(paste("LDA (AUC =", round(auc(roc_lda), 3), ")"),
                          paste("QDA (AUC =", round(auc(roc_qda), 3), ")")),
         col = c("#2E86AB", "#E94F37"), lwd = 2)
} else {
  plot(roc_lda, col = "#2E86AB", lwd = 2, main = "ROC Curve: LDA")
  abline(a = 0, b = 1, lty = 2, col = "gray")
  legend("bottomright", paste("LDA (AUC =", round(auc(roc_lda), 3), ")"), col = "#2E86AB", lwd = 2)
}
dev.off()

# ============================================================================
# SECTION 7: CONFUSION MATRIX VISUALIZATION
# ============================================================================
# Confusion matrix shows:
#   - True Positives (TP): Correctly classified as positive
#   - True Negatives (TN): Correctly classified as negative
#   - False Positives (FP): Incorrectly classified as positive (Type I error)
#   - False Negatives (FN): Incorrectly classified as negative (Type II error)

plot_cm <- function(cm, title) {
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Reference", "Freq")
  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile() + geom_text(aes(label = Freq), size = 8, color = "white") +
    scale_fill_gradient(low = "#E94F37", high = "#2E86AB") +
    labs(title = title) + theme_minimal() + theme(legend.position = "none")
}

p_lda <- plot_cm(lda_cm, paste("LDA\nAcc:", round(lda_cm$overall["Accuracy"], 3)))
if (qda_success) {
  p_qda <- plot_cm(qda_cm, paste("QDA\nAcc:", round(qda_cm$overall["Accuracy"], 3)))
  ggsave("figures/classification_confusion_matrices.png", grid.arrange(p_lda, p_qda, ncol = 2), width = 12, height = 5, dpi = 150)
} else {
  ggsave("figures/classification_confusion_matrices.png", p_lda, width = 6, height = 5, dpi = 150)
}

# ============================================================================
# SECTION 8: LDA PROJECTION
# ============================================================================
# Projects all observations onto the first linear discriminant (LD1)
# Shows how well the discriminant function separates the two groups
# Good separation: minimal overlap between density curves

lda_scores <- predict(lda_model, df_class)$x
df_lda <- data.frame(Country = df_class$Country, Status = df_class$Status, LD1 = lda_scores[,1])

p_proj <- ggplot(df_lda, aes(x = LD1, fill = Status)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
  labs(title = "LDA Projection", x = "LD1") + theme_minimal()
ggsave("figures/lda_projection.png", p_proj, width = 10, height = 6, dpi = 150)

# 2D classification plot showing decision space
p_2d <- ggplot(df_class, aes(x = Life_expectancy, y = Schooling, color = Status)) +
  geom_point(alpha = 0.7, size = 2.5) + stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = "Classification Space") + theme_minimal()
ggsave("figures/classification_2d.png", p_2d, width = 10, height = 8, dpi = 150)

# ============================================================================
# SECTION 9: PARTITION PLOTS
# ============================================================================
# Partition plots show classification boundaries for pairs of variables
# Useful for visualizing how LDA/QDA divide the feature space
# Note: Only uses 2 variables at a time (projection of full decision boundary)

cat("\n=== PARTITION PLOTS ===\n")
cat("Creating LDA partition plots for key variable pairs...\n")

# Select two most important variables for visualization
top2_vars <- c("Life_expectancy", "Schooling")

png("figures/classification_partition_lda.png", width = 800, height = 700, res = 120)
partimat(Status ~ Life_expectancy + Schooling,
         data = df_class, method = "lda",
         main = "LDA Partition Plot")
dev.off()

png("figures/classification_partition_qda.png", width = 800, height = 700, res = 120)
tryCatch({
  partimat(Status ~ Life_expectancy + Schooling,
           data = df_class, method = "qda",
           main = "QDA Partition Plot")
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, "QDA partition plot failed\n(possible singularity)", cex = 1.2)
})
dev.off()

cat("Saved: figures/classification_partition_lda.png\n")
cat("Saved: figures/classification_partition_qda.png\n")

# ============================================================================
# SECTION 10: MISCLASSIFICATION ANALYSIS
# ============================================================================
# Examine which countries were misclassified in the TEST SET
# (Using test set for honest error estimate, not train set)
# Posterior probability shows model's confidence in its prediction

cat("\n=== MISCLASSIFIED (Test Set) ===\n")
test_misclass <- test_data %>%
  mutate(Predicted = lda_pred$class,
         Prob_Developed = round(lda_pred$posterior[, "Developed"], 3)) %>%
  dplyr::filter(Status != Predicted) %>%
  dplyr::select(Country, Status, Predicted, Prob_Developed)

if (nrow(test_misclass) > 0) {
  print(test_misclass)
} else {
  cat("No misclassifications in test set\n")
}
write.csv(test_misclass, "figures/lda_misclassified.csv", row.names = FALSE)

# ============================================================================
# SECTION 11: VARIABLE IMPORTANCE
# ============================================================================
# Standardized discriminant coefficients: |LD1 coefficient| × variable SD
# This heuristic accounts for different variable scales
# Higher values indicate stronger discrimination between groups

cat("\n=== VARIABLE IMPORTANCE ===\n")
cat("(Using |LD1 coefficient| × variable SD - a common heuristic)\n\n")
sd_vars <- apply(df_class[, pred_vars], 2, sd)
var_imp <- data.frame(
  Variable = names(sd_vars),
  LD1_Coef = round(lda_model$scaling[, 1], 4),
  SD = round(sd_vars, 4),
  Importance = round(abs(lda_model$scaling[, 1]) * sd_vars, 4)
) %>%
  arrange(desc(Importance))

print(var_imp)

p_imp <- ggplot(var_imp, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Variable Importance (LDA)", x = "") + theme_minimal()
ggsave("figures/classification_var_importance.png", p_imp, width = 10, height = 6, dpi = 150)

# ============================================================================
# SECTION 12: OUTLIER SENSITIVITY ANALYSIS
# ============================================================================
# Examine whether classification accuracy is sensitive to outliers
# This is important because outliers can distort discriminant coefficients

cat("\n=== OUTLIER SENSITIVITY ANALYSIS ===\n")

# Identify multivariate outliers using Mahalanobis distance
p_vars <- length(pred_vars)
cutoff_mahal <- qchisq(0.975, df = p_vars)
mahal_dist <- mahalanobis(X_class_matrix, colMeans(X_class_matrix), cov(X_class_matrix))
outlier_idx <- which(mahal_dist > cutoff_mahal)

cat("Mahalanobis cutoff (chi-sq, df=", p_vars, ", alpha=0.025):", round(cutoff_mahal, 2), "\n")
cat("Outliers identified:", length(outlier_idx), "countries\n")

if (length(outlier_idx) > 0 && length(outlier_idx) <= 15) {
  cat("Countries:", paste(df_class$Country[outlier_idx], collapse = ", "), "\n")
} else if (length(outlier_idx) > 15) {
  cat("First 15:", paste(df_class$Country[outlier_idx[1:15]], collapse = ", "), "...\n")
}

# Re-run classification WITHOUT outliers
if (length(outlier_idx) > 0 && length(outlier_idx) < nrow(df_class) * 0.2) {
  df_no_outliers <- df_class[-outlier_idx, ]

  # New train-test split without outliers
  set.seed(123)
  train_idx_no <- createDataPartition(df_no_outliers$Status, p = 0.7, list = FALSE)
  train_no <- df_no_outliers[train_idx_no, ]
  test_no <- df_no_outliers[-train_idx_no, ]

  # Retrain LDA
  lda_no <- lda(Status ~ ., data = train_no[, c("Status", pred_vars)])
  lda_pred_no <- predict(lda_no, test_no)
  lda_acc_no <- mean(lda_pred_no$class == test_no$Status)

  # ROC for model without outliers
  roc_lda_no <- roc(test_no$Status, lda_pred_no$posterior[, "Developed"],
                    levels = c("Developing", "Developed"), quiet = TRUE)
  auc_no <- auc(roc_lda_no)

  cat("\nLDA Performance Comparison:\n")
  cat("  WITH outliers:    Accuracy =", round(lda_cm$overall["Accuracy"], 4),
      ", AUC =", round(auc(roc_lda), 4), "\n")
  cat("  WITHOUT outliers: Accuracy =", round(lda_acc_no, 4),
      ", AUC =", round(auc_no, 4), "\n")

  # Check sensitivity
  acc_diff <- abs(lda_cm$overall["Accuracy"] - lda_acc_no)
  auc_diff <- abs(auc(roc_lda) - auc_no)

  if (acc_diff < 0.05 && auc_diff < 0.05) {
    cat("  Conclusion: Results are ROBUST (accuracy/AUC change < 0.05)\n")
    outlier_status <- "Robust"
  } else {
    cat("  Conclusion: Results are SENSITIVE to outliers\n")
    cat("  (accuracy/AUC changed by", round(max(acc_diff, auc_diff), 3), ")\n")
    outlier_status <- "Sensitive"
  }

  # QDA sensitivity (if it succeeded)
  if (qda_success) {
    tryCatch({
      qda_no <- qda(Status ~ ., data = train_no[, c("Status", pred_vars)])
      qda_pred_no <- predict(qda_no, test_no)
      qda_acc_no <- mean(qda_pred_no$class == test_no$Status)
      roc_qda_no <- roc(test_no$Status, qda_pred_no$posterior[, "Developed"],
                        levels = c("Developing", "Developed"), quiet = TRUE)

      cat("\nQDA Performance Comparison:\n")
      cat("  WITH outliers:    Accuracy =", round(qda_cm$overall["Accuracy"], 4),
          ", AUC =", round(auc(roc_qda), 4), "\n")
      cat("  WITHOUT outliers: Accuracy =", round(qda_acc_no, 4),
          ", AUC =", round(auc(roc_qda_no), 4), "\n")
    }, error = function(e) {
      cat("\nQDA without outliers failed - likely rank deficiency\n")
    })
  }

} else if (length(outlier_idx) == 0) {
  cat("No outliers detected - sensitivity analysis not needed\n")
  outlier_status <- "N/A"
} else {
  cat("Too many outliers (>20%) - sensitivity analysis not meaningful\n")
  outlier_status <- "N/A"
}

# ============================================================================
# SECTION 13: STATISTICAL VALIDITY SUMMARY
# ============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("CLASSIFICATION VALIDITY SUMMARY\n")
cat(strrep("=", 60), "\n")

# Box's M status
if (is.infinite(box_m_result$statistic) || is.na(box_m_result$p.value)) {
  boxm_status <- "Indeterminate (singularity)"
} else if (box_m_result$p.value < 0.05) {
  boxm_status <- "Violated - use QDA"
} else {
  boxm_status <- "Met - LDA appropriate"
}
cat("Covariance Homogeneity:    ", boxm_status, "\n")
cat("Outlier Sensitivity:       ", outlier_status, "\n")
cat("LDA Test Accuracy:         ", round(lda_cm$overall["Accuracy"], 4), "\n")
cat("LDA Test AUC:              ", round(auc(roc_lda), 4), "\n")
if (qda_success) {
  cat("QDA Test Accuracy:         ", round(qda_cm$overall["Accuracy"], 4), "\n")
  cat("QDA Test AUC:              ", round(auc(roc_qda), 4), "\n")
}
cat(strrep("=", 60), "\n")

cat("\n=== Classification Complete ===\n")
cat("Run 07_clustering.R next.\n")

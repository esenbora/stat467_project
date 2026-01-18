# STAT 467 - Classification (LDA/QDA)
# Input: data_country_level.csv

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

pred_vars <- c("Life_expectancy", "Adult_Mortality", "Schooling", "GDP_log",
               "Income_composition", "BMI", "HIV_AIDS", "Diphtheria", "thinness_10_19_years")

df_class <- df %>% select(Country, Status, all_of(pred_vars)) %>% drop_na()
cat("Complete cases:", nrow(df_class), "\n\n")

# Train-test split
set.seed(123)
train_idx <- createDataPartition(df_class$Status, p = 0.7, list = FALSE)
train_data <- df_class[train_idx, ]
test_data <- df_class[-train_idx, ]

cat("Train:", nrow(train_data), ", Test:", nrow(test_data), "\n\n")

# LDA ASSUMPTION TESTS
cat("=== LDA ASSUMPTION TESTS ===\n")

# 1. Box's M Test for equality of covariance matrices
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

# 2. Multivariate Normality by Group
cat("\n2. Multivariate Normality by Group (Mardia's Test):\n")
for (grp in levels(df_class$Status)) {
  cat("\n  ", grp, ":\n", sep = "")
  grp_data <- df_class %>% dplyr::filter(Status == grp) %>%
    dplyr::select(all_of(pred_vars))
  tryCatch({
    mvn_result <- MVN::mvn(grp_data, mvnTest = "mardia")
    cat("    Skewness p =", format(mvn_result$multivariateNormality$`p value`[1],
                                   digits = 3), "\n")
    cat("    Kurtosis p =", format(mvn_result$multivariateNormality$`p value`[2],
                                   digits = 3), "\n")
    if (mvn_result$multivariateNormality$Result[1] == "NO" ||
        mvn_result$multivariateNormality$Result[2] == "NO") {
      cat("    WARNING: MVN may be violated for this group\n")
    } else {
      cat("    MVN assumption met\n")
    }
  }, error = function(e) {
    cat("    MVN test failed:", conditionMessage(e), "\n")
  })
}

# Calculate class proportions for priors
class_props <- prop.table(table(train_data$Status))
cat("\n3. Class Proportions (for prior probabilities):\n")
print(round(class_props, 3))

# LDA
cat("\n=== LDA ===\n")
# Use observed class proportions as priors (more realistic than equal priors)
lda_model <- lda(Status ~ ., data = train_data[, c("Status", pred_vars)],
                 prior = class_props)
print(lda_model)

lda_pred <- predict(lda_model, test_data)
lda_cm <- confusionMatrix(lda_pred$class, test_data$Status)
print(lda_cm)

lda_coef <- data.frame(Variable = rownames(lda_model$scaling), LD1 = round(lda_model$scaling[,1], 4))
lda_coef <- lda_coef %>% arrange(desc(abs(LD1)))
write.csv(lda_coef, "figures/lda_coefficients.csv", row.names = FALSE)

# QDA
cat("\n=== QDA ===\n")
qda_success <- FALSE
tryCatch({
  # Use same priors as LDA for consistency
  qda_model <- qda(Status ~ ., data = train_data[, c("Status", pred_vars)],
                   prior = class_props)
  qda_pred <- predict(qda_model, test_data)
  qda_cm <- confusionMatrix(qda_pred$class, test_data$Status)
  print(qda_cm)
  qda_success <- TRUE
}, error = function(e) {
  cat("QDA failed (rank deficiency - group size too small for # predictors)\n")
  cat("Skipping QDA, using LDA only\n")
})

# Cross-validation
cat("\n=== CROSS-VALIDATION ===\n")
cv_ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

lda_cv <- train(Status ~ ., data = df_class[, c("Status", pred_vars)], method = "lda", trControl = cv_ctrl, metric = "ROC")
cat("LDA CV ROC:", round(lda_cv$results$ROC, 4), "\n")

if (qda_success) {
  qda_cv <- train(Status ~ ., data = df_class[, c("Status", pred_vars)], method = "qda", trControl = cv_ctrl, metric = "ROC")
  cat("QDA CV ROC:", round(qda_cv$results$ROC, 4), "\n")
}

# ROC curves
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

# Confusion matrix plots
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

# LDA projection
lda_scores <- predict(lda_model, df_class)$x
df_lda <- data.frame(Country = df_class$Country, Status = df_class$Status, LD1 = lda_scores[,1])

p_proj <- ggplot(df_lda, aes(x = LD1, fill = Status)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Developed" = "#2E86AB", "Developing" = "#E94F37")) +
  labs(title = "LDA Projection", x = "LD1") + theme_minimal()
ggsave("figures/lda_projection.png", p_proj, width = 10, height = 6, dpi = 150)

# 2D classification plot
p_2d <- ggplot(df_class, aes(x = Life_expectancy, y = Schooling, color = Status)) +
  geom_point(alpha = 0.7, size = 2.5) + stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("#2E86AB", "#E94F37")) +
  labs(title = "Classification Space") + theme_minimal()
ggsave("figures/classification_2d.png", p_2d, width = 10, height = 8, dpi = 150)

# Partition plot using klaR::partimat
# Shows classification boundaries for variable pairs
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

# Misclassified countries (using TEST SET for honest error estimate)
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

# Variable importance (standardized discriminant coefficients)
# Note: This is |coefficient| × SD, a common heuristic for importance
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

cat("\n=== Classification Complete ===\n")
cat("Run 07_clustering.R next.\n")


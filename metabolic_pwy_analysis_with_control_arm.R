
###################################################################################
## Author: Peter Mangani
## Date: September 2023
## Purpose : Script for analyzing MetaCyc metabolic pathways obtained from PICRUSt2
## Project: Msc Bioinformatics Project
#################################################################################

## initialization
rm(list=ls()) ## Don't execute this line if there are  objects in workspace that you need later.

library(tidyverse)
library(caret)
library(stringr)
library(ggplot2)
library(glmnet)
library(reshape2)
library(dplyr) 

## Find the path to the directory where the current script is located
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

############### Input metadata and MetaCyC Pathway abundance table #################

## Import ASV table
pwy <- read.table(paste0(current_dir,"/Analysis/path_abun_unstrat.tsv"),sep="\t",header=T,row.names=1, stringsAsFactors=TRUE)

## Import Sample Metadata
metadata <- read.table(paste0(current_dir,"/Analysis/sample_metadata.tsv"),sep="\t",header=T,row.names=1, stringsAsFactors=TRUE)

## Identify common samples and ensure identical ordering
common_samples <- intersect(rownames(metadata), colnames(pwy))
pwy <- pwy[,common_samples]
metadata <- metadata[common_samples,]

pwy[is.na(pwy)] <- 0 ## replace NAs with zero
features <- data.frame(pwy)
features <- t(features)

metadata[is.na(metadata)] <- 0 ## replace NAs with 0
set.seed(1000)

## Convert labels to double for LASSO
label_laz <- as.double(metadata$LAZ)

##################################################################################
####     Perform LASSO Regression with Cross Validation
##################################################################################

cv_model <- cv.glmnet(features, label_laz , alpha = 1, nfolds = 10)
plot(cv_model) ## Plot the cross validation model

best_lambda <- cv_model$lambda.min ## Extract best lambda value
lasso_model <- glmnet(features, label_laz, alpha = 1, lambda = best_lambda)

## Obtain coefficients and isolate non-zero pathways
selected_feature <- coef(lasso_model)
nonzero_coef_pathways <- which(selected_feature != 0, arr.ind = TRUE)

selected_feature_indices <- unique(nonzero_coef_pathways[, "row"])
selected_feature_indices <- selected_feature_indices - 1 

X_test_selected <- features[, selected_feature_indices]

##################################################################################
## Perform Univariate (Wilcoxon) Analysis on Selected Features (Stunted vs Healthy)
##################################################################################

healthy_row_names <- row.names(metadata %>% filter(LAZ >= -2))
stunted_row_names <- row.names(metadata %>% filter(LAZ < -2))

stunted_pathways <- as.data.frame(X_test_selected[row.names(X_test_selected) %in% stunted_row_names,])
healthy_pathways <- as.data.frame(X_test_selected[row.names(X_test_selected) %in% healthy_row_names,])
w_null_expected <- (nrow(healthy_pathways) * nrow(stunted_pathways)) / 2

prob_value <- vector("numeric", length = ncol(healthy_pathways)) 
pathways <- vector("character", length = ncol(healthy_pathways)) 
t_value <- vector("character", length = ncol(healthy_pathways)) 
counter <- 1 

for (variable in colnames(healthy_pathways)) {
  healthy_values <- healthy_pathways[variable]    
  stunted_values <- stunted_pathways[variable]
  
  test_result <- wilcox.test(healthy_values[,1], stunted_values[,1])
  
  prob_value[counter] <- test_result$p.value
  pathways[counter] <- variable
  t_value[counter] <- test_result$statistic
  counter <- counter + 1 
}

adjusted_p_values <- p.adjust(prob_value , method = "BH")

univariate_results <- data.frame(
  pathway = pathways, 
  p_value = prob_value, 
  adj_pvalue = adjusted_p_values, 
  W_Statistic = t_value,
  Expected_Null_W = w_null_expected,
  stringsAsFactors = FALSE
  
  )
print(univariate_results)
# 
#    pathway       p_value     adj_pvalue       W_Statistic   Expected_Null_W
# 1  POLYAMSYN-PWY 0.003630298 0.03630298       39942         34843.5
# 2       PWY-4722 0.018335777 0.09167888       38063         34843.5
# 3       PWY-6470 0.081033206 0.18633410       31785         34843.5
# 4       PWY-6749 0.232645474 0.38774246       36936         34843.5
# 5       PWY-6901 0.471159743 0.47115974       36107         34843.5
# 6       PWY-7031 0.402449663 0.44716629     33402.5         34843.5
# 7       PWY-7332 0.346233892 0.43279236       36470         34843.5
# 8       PWY-7373 0.289444008 0.41349144       33006         34843.5
# 9       PWY-7374 0.093167051 0.18633410       37512         34843.5
# 10     PWY0-1533 0.046581406 0.15527135       38332         34843.5


# Add a clear interpretation column showing the biological direction
univariate_results$Biological_Direction <- ifelse(
  univariate_results$W_Statistic > univariate_results$Expected_Null_W,
  "Reduced in Stunted Cohort",
  "Enriched in Stunted Cohort"
)

# Print the final consolidated table to the console for easy comparison
print("--- Wilcoxon Rank-Sum Results (All 10 Pathways) ---")
print(univariate_results, row.names = FALSE)

# pathway       p_value       adj_pvalue    W_Statistic Expected_Null_W   Biological_Direction
# POLYAMSYN-PWY 0.003630298   0.03630298       39942         34843.5      Reduced in Stunted Cohort
# PWY-4722      0.018335777   0.09167888       38063         34843.5      Reduced in Stunted Cohort
# PWY-6470      0.081033206   0.18633410       31785         34843.5      Enriched in Stunted Cohort
# PWY-6749      0.232645474   0.38774246       36936         34843.5      Reduced in Stunted Cohort
# PWY-6901      0.471159743   0.47115974       36107         34843.5      Reduced in Stunted Cohort
# PWY-7031      0.402449663   0.44716629     33402.5         34843.5      Enriched in Stunted Cohort
# PWY-7332      0.346233892   0.43279236       36470         34843.5      Reduced in Stunted Cohort
# PWY-7373      0.289444008   0.41349144       33006         34843.5      Enriched in Stunted Cohort
# PWY-7374      0.093167051   0.18633410       37512         34843.5      Reduced in Stunted Cohort
# PWY0-1533     0.046581406   0.15527135       38332         34843.5      Reduced in Stunted Cohort


##################################################################################
## Rigorous Multivariable Logistic Modeling & Interaction Q-Values
##################################################################################

# Establish the clean binary clinical groups 
pwy_names_in_df <- colnames(X_test_selected)
linkage_df <- data.frame(
  SampleID = rownames(metadata),
  LAZ = label_laz,
  InterventionCode = as.factor(metadata$InterventionCode),
  Polyamines = X_test_selected[, pwy_names_in_df[stringr::str_detect(colnames(X_test_selected), "POLYAMSYN")]],
  PWY_4722 = X_test_selected[, pwy_names_in_df[stringr::str_detect(colnames(X_test_selected), "PWY.4722")]]
)

linkage_df$StuntingStatus <- factor(
  ifelse(linkage_df$LAZ < -2, "Stunted", "Non-Stunted"),
  levels = c("Non-Stunted", "Stunted")
)

# Reference group set to Arm A 
linkage_df$InterventionCode <- relevel(linkage_df$InterventionCode, ref = "A")

# Fit Logistic Interaction Models
logit_polyam <- glm(StuntingStatus ~ Polyamines * InterventionCode, data = linkage_df, family = binomial)
logit_4722   <- glm(StuntingStatus ~ PWY_4722 * InterventionCode, data = linkage_df, family = binomial)

# Extract and compile interaction terms for Multiple-Testing Correction
coef_polyam <- as.data.frame(summary(logit_polyam)$coefficients)
coef_4722   <- as.data.frame(summary(logit_4722)$coefficients)

interaction_data <- rbind(
  data.frame(Pathway = "POLYAMSYN-PWY", Term = "Polyamines:InterventionCodeB", coef_polyam["Polyamines:InterventionCodeB", ]),
  data.frame(Pathway = "POLYAMSYN-PWY", Term = "Polyamines:InterventionCodeC", coef_polyam["Polyamines:InterventionCodeC", ]),
  data.frame(Pathway = "PWY-4722", Term = "PWY_4722:InterventionCodeB", coef_4722["PWY_4722:InterventionCodeB", ]),
  data.frame(Pathway = "PWY-4722", Term = "PWY_4722:InterventionCodeC", coef_4722["PWY_4722:InterventionCodeC", ])
)

colnames(interaction_data)[3:6] <- c("Estimate_LogOdds", "Std_Error", "Z_Value", "P_Value")
interaction_data$Q_Value <- p.adjust(interaction_data$P_Value, method = "BH")

print("--- Final Cohort Interaction P-Values and Q-Values ---")
print(interaction_data[, c("Pathway", "Term", "Estimate_LogOdds", "P_Value", "Q_Value")])

 
# Pathway          Term                            Estimate_LogOdds  P_Value    Q_Value
# POLYAMSYN-PWY    Polyamines:InterventionCodeB    -0.0001929111     0.5329693  0.6441544
# POLYAMSYN-PWY    Polyamines:InterventionCodeC     0.0001348981     0.6441544  0.6441544
# PWY-4722         PWY_4722:InterventionCodeB     0.0179151478       0.4831167  0.6441544
# PWY-4722         PWY_4722:InterventionCodeC     0.0327221198       0.2859895  0.6441544


##################################################################################
## Publication-Ready Grouped Visualization
##################################################################################

cohort_plot_df <- linkage_df %>%
  select(InterventionCode, StuntingStatus, Polyamines, PWY_4722) %>%
  melt(id.vars = c("InterventionCode", "StuntingStatus"))

levels(cohort_plot_df$variable) <- c(
  "Polyamines" = "Polyamine Biosynthesis I",
  "PWY_4722" = "Creatinine Degradation II"
)

final_plot <- ggplot(cohort_plot_df, aes(x = InterventionCode, y = value, fill = StuntingStatus)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.6, color = "#2c3e50",
               position = position_dodge(width = 0.75)) +
  geom_point(aes(color = StuntingStatus), 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
             alpha = 0.35, size = 1.2) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_manual(values = c("Non-Stunted" = "#5dade2", "Stunted" = "#ec7063")) +
  scale_color_manual(values = c("Non-Stunted" = "#2e86c1", "Stunted" = "#cb4335")) +
  labs(
    title = "Functional Pathway Abundances Stratified by Trial Arm and Stunting Status",
    x = "Blinded Intervention Code",
    y = "Relative Functional Abundance",
    fill = "Legend",
    color = "Legend"
  ) +
  theme_classic(base_size = 12) + 
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 15)),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_blank(), 
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "gray95") 
  )

# Output plot to screen
print(final_plot)

# Automatically save plot in your directory
ggsave(
  filename = paste0(current_dir, "/Analysis/functional_pathway_abundances_30months.tiff"),
  plot = final_plot,
  device = "tiff",
  dpi = 300,
  width = 8.5,
  height = 5.5,
  units = "in",
  compression = "lzw"
)


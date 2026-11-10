
## ==============================================================================
## Microbiome Classification and Differential Abundance Analysis
## Author: Peter Mangani
## Date: September 2023
## Purpose : Script for taxonomic data
## Project: Msc Bioinformatics Project
# # ==============================================================================

## Initialization
rm(list = ls()) # Uncomment if you need to clear workspace

library(tidyverse)
library(caret)
library(randomForest)
library(pROC)
library(stringr)
library(missForest)

# path detection for both RStudio and command-line execution
current_dir <- tryCatch({
  dirname(rstudioapi::getSourceEditorContext()$path)
}, error = function(e) {
  getwd()
})

############### Input metadata, ASV and Taxonomy table #################

## import ASV table
asv <- read.table(file.path(current_dir, "/Analysis/table_from_biom.tsv"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = TRUE)
dim(asv) ## [1] 4377  529

## Import Sample Metadata
metadata <- read.table(file.path(current_dir, "/Analysis/sample_metadata.tsv"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = TRUE)
dim(metadata) ## [1] 536  13

## import taxonomy
taxonomy <- read.table(file.path(current_dir, "/Analysis/taxa_metadata.tsv"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = TRUE)
dim(taxonomy) ## [1] 5689  3

## perform inner join on asv with taxonomy 
asv_taxa <- inner_join(asv, taxonomy, by = "asv")
dim(asv_taxa) ## [1] 4377  531


## Binning & Genus-Level Aggregation
tax_table <- str_split_fixed(asv_taxa$Taxon, ";", n = 8)

asv_taxa <- data.frame(genus = tax_table[, 6], asv_taxa)

asv_taxa <- asv_taxa[!(asv_taxa$genus == "g__"), ] 
asv_taxa <- asv_taxa[!(asv_taxa$genus == ""), ] 

txonID <- which(colnames(asv_taxa) == "Taxon") 
asv_taxa <- asv_taxa[, -txonID] 

confID <- which(colnames(asv_taxa) == "Confidence") 
asv_taxa <- asv_taxa[, -confID] 

row.names(asv_taxa) <- asv_taxa[, 2]

otu <- asv_taxa[, -2]
dim(otu) ## [1] 2338  530

otu <- aggregate(otu[, -1], by = list(otu$genus), sum)
dim(otu) ## [1] 329  530

rownames(otu) <- otu[, 1]
otu <- otu[, -1]
dim(otu) ## [1] 329  529

y <- t(otu)

## Identify common samples between metadata and otu table
common_samples <- intersect(rownames(metadata), colnames(otu))
length(common_samples) # 528

otu <- otu[, common_samples]
metadata <- metadata[common_samples, ]

all(rownames(metadata) == colnames(otu)) # TRUE


## Normalization 
#calculate relative abundance of each genus by dividing its value by the total reads per sample
otu <- sweep(otu, 2, colSums(otu), "/")
otu[is.na(otu)] <- 0 

## Prepare for training
x <- data.matrix(otu)
x <- t(x)

metadata[is.na(metadata)] <- 0 
metadata$Identifier <- as.factor(metadata$Identifier)


################ Train and Test ###############

set.seed(1000)
train_index <- createDataPartition(metadata$Identifier, p = 0.8, list = FALSE) 

x.train <- x[train_index, ] 
y.train <- metadata$Identifier[train_index]

x.test <- x[-train_index, ]
y.test <- metadata$Identifier[-train_index]

train_control <- trainControl(
  method = "cv",
  number = 5, 
  summaryFunction = twoClassSummary, 
  classProbs = TRUE 
)

set.seed(1000)
rf_train <- train(
  x = x.train, 
  y = as.factor(y.train),
  method = 'rf',
  metric = "ROC", 
  trControl = train_control
)
rf_train

rf_test <- predict(rf_train, x.test) 
conf_matrix <- confusionMatrix(rf_test, y.test)
conf_matrix$table
conf_matrix$byClass

rf_test <- predict(rf_train, x.test, type = "prob")
rf_test <- rf_test[, 1]


############ Plot ROC curve ############
rf <- roc(y.test, rf_test) 
auc <- rf$auc
cat("Area under the curve:", auc, "\n")

pdf(file.path(current_dir, "ROC_curve_taxonomy.pdf"))
plot(rf, col = "blue", legacy.axes = TRUE, main = "ROC Curve Plot - Bacteria taxa abundances", 
     xlab = "False Positive Rate", ylab = "True Positive Rate", print.auc = TRUE)
legend("bottomright", legend = c("Random Forest", "Reference Line(Random Chance)"), col = c("blue", "grey"), lwd = 2)
dev.off()


########## Feature Importance & Differential Abundance ###########
impVars <- varImp(rf_train) 
ImpMeasure <- data.frame(impVars$importance)
ImpMeasure <- ImpMeasure[order(ImpMeasure$Overall, decreasing = TRUE), , drop = FALSE]
ImpMeasure_top30 <- ImpMeasure[1:nrow(ImpMeasure), , drop = FALSE]

healthy_row_names <- row.names(metadata %>% filter(LAZ >= -2))
stunted_row_names <- row.names(metadata %>% filter(LAZ < -2))

stunted_taxa <- as.data.frame(x[row.names(x) %in% stunted_row_names, ])
healthy_taxa <- as.data.frame(x[row.names(x) %in% healthy_row_names, ])

p_values <- vector("numeric", length = ncol(healthy_taxa))
counter <- 1

for (variable in colnames(healthy_taxa)) {
  healthy_values <- healthy_taxa[variable]
  stunted_values <- stunted_taxa[variable]  
  
  test_result <- wilcox.test(healthy_values[, 1], stunted_values[, 1])
  p_values[counter] <- test_result$p.value
  counter <- counter + 1 
}
names(p_values) <- colnames(healthy_taxa)

adjusted_p_values <- p.adjust(p_values, method = "BH")

ggplot(ImpMeasure_top30, aes(x = reorder(rownames(ImpMeasure_top30), Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "blue", width = 0.5) +
  coord_flip() +
  labs(x = "Taxa", y = "Importance", title = "Feature Importance Plot")

stunted_taxa1 <- as.data.frame(y[row.names(y) %in% stunted_row_names, ])
healthy_taxa1 <- as.data.frame(y[row.names(y) %in% healthy_row_names, ])

p_taxa_stunted <- vector("numeric", length = nrow(ImpMeasure_top30))
p_taxa_healthy <- vector("numeric", length = nrow(ImpMeasure_top30))
counter1 <- 1

for (name in row.names(ImpMeasure_top30)){
  stunted_val <- colSums(stunted_taxa1[name], na.rm = FALSE, dims = 1)
  healthy_val <- colSums(healthy_taxa1[name], na.rm = FALSE, dims = 1)
  
  p_taxa_stunted[counter1] <- stunted_val
  p_taxa_healthy[counter1] <- healthy_val
  counter1 <- counter1 + 1
}

results <- data.frame(
  Overall = ImpMeasure_top30$Overall,
  non_stunted_count = p_taxa_healthy,
  stunted_count = p_taxa_stunted,
  p_value = p_values,
  adjusted_p_value = adjusted_p_values,
  row.names = row.names(ImpMeasure_top30)
)

print(results)

results_filtered <- results[results$p_value <= 0.05, ]
print(results_filtered)

results_adj_filtered <- results[results$adjusted_p_value <= 0.05, ]
print(results_adj_filtered)

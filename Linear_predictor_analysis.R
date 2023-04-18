####Script for linear predictor analysis

#Read in top 20 genes from caret model
top_20_genes_knnProfile <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/top_20_genes_knnProfile.txt", quote="", comment.char="")

#Need to run setup in caret_models_script.r for "x"
exp_data_training <- t(x)
dim(exp_data_training)


# extract Limma stats for the top 20 genes
tables_top_genes <- na.omit(Limma_DE_genes_Tum_met_vs_Tum_no_met[match(top_20_genes_knnProfile$V1, Limma_DE_genes_Tum_met_vs_Tum_no_met$Gene),])
rownames(tables_top_genes) <- tables_top_genes$Gene

###############################################################################

##### Testing the linear predictor in both training and testing datasets

###############################################################################


##### testing in training dataset first
exp_data_training_topgenes <- exp_data_training[match(top_20_genes_knnProfile$V1, rownames(exp_data_training)),]
linear_pred_training <- vector()

for(i in 1:length(colnames(exp_data_training_topgenes))) {
  linear_pred_training[i] = 0
}

#work through columns (patients), and work through each row (gene)
for(i in 1:length(colnames(exp_data_training_topgenes))) {
  for(j in 1:length(rownames(exp_data_training_topgenes))) {
    linear_pred_training[i] <- linear_pred_training[i] + tables_top_genes[j,2] * exp_data_training_topgenes[j,i]
  }
}

names(linear_pred_training) <- colnames(exp_data_training_topgenes)

##Supplemental Fig 3 (Training set)
#Boxplot for tum_no_met vs tum_met
cmpr <- list(c("N", "Y"))

predictor_training <- data.frame(linear_predictor= linear_pred_training, group=training$group)
p <- ggplot(predictor_training, aes(x=group, y=linear_predictor)) + geom_boxplot()
p + geom_jitter(shape=16, position=position_jitter(0.2))
p + stat_compare_means(comparisons = cmpr, tip.length=0.01,label = "p.format") + geom_jitter(shape=16, position=position_jitter(0.2)) + labs(title="Training set",x="Group", y = "Linear predictors")



###############################################################################

##### next, testing in testing dataset
exp_data_testing <- testing[,c(1:1274)]
exp_data_testing <- t(exp_data_testing)

exp_data_testing_topgenes <- exp_data_testing[match(top_20_genes_knnProfile$V1, rownames(exp_data_testing)),]
linear_pred_testing <- vector()
for(i in 1:length(colnames(exp_data_testing_topgenes))) {
  linear_pred_testing[i] = 0
}

#work through columns (patients), and work through each row (gene)
for(i in 1:length(colnames(exp_data_testing_topgenes))) {
  for(j in 1:length(rownames(exp_data_testing_topgenes))) {
    linear_pred_testing[i] <- linear_pred_testing[i] + tables_top_genes[j,2] * exp_data_testing_topgenes[j,i]
  }
}

names(linear_pred_testing) <- colnames(exp_data_testing_topgenes)

#Supplemental fig 3 (Testing set)
#Boxplot for tum_no_met vs tum_met

cmpr <- list(c("N", "Y"))

predictor_testing <- data.frame(linear_predictor=linear_pred_testing, group=testing$group)
p <- ggplot(predictor_testing, aes(x=group, y=linear_predictor)) + geom_boxplot()
p + geom_jitter(shape=16, position=position_jitter(0.2))
p + stat_compare_means(comparisons = cmpr, tip.length=0.01,label = "p.format") + geom_jitter(shape=16, position=position_jitter(0.2)) + labs(title="Testing set",x="Group", y = "Linear predictors")


###############################################################################

## ROC plots to show model performance

## Figure 3 

library(pROC)


roc1 <- roc(testing$group, linear_pred_testing, percent=TRUE)

plot(roc1)

roc1 <- roc(training$group, linear_pred_training, percent=TRUE)
ci(roc1)
roc2 <- roc(testing$group, linear_pred_testing, percent=TRUE)
ci(roc2)

auc(roc1)
auc(roc2)

#Manually added the AUC values from above outputs
plot(roc1, col="red", main="ROC curves of linear predictor", xlim=c(100,0))
plot(roc2, add=TRUE, col="blue")
legend("bottomright", legend=c("Training, AUC=0.85 (0.80-0.91)", "Testing, AUC=0.88 (0.78-0.99)"),
       col=c("red", "blue"), lty=1, cex=0.7)

###############################################################################








###############################################################################
## Supplemental Fig 4
## Boxplots of linear predictors across the four clinical groups in our cohort

#List of all IDs (immuno-suppressed excluded)
Tum_no_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Tum_no_met_IDs.txt", quote="\"", comment.char="")
Tum_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Tum_met_IDs.txt", quote="\"", comment.char="")
Norm_with_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Norm_with_met_IDs.txt", quote="\"", comment.char="")
Norm_of_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Norm_of_met_IDs.txt", quote="\"", comment.char="")
Norm_no_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Norm_no_met_IDs.txt", quote="\"", comment.char="")
Met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Met_IDs.txt", quote="\"", comment.char="")
exp_data.normalised_no_IS_merged_clean.all_532 <- read.delim("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/exp_data.normalised_no_IS_merged_clean.all_532.txt", row.names=1, quote="")

sanofi_NT <- rbind(Norm_no_met_IDs, Norm_with_met_IDs, Tum_met_IDs, Tum_no_met_IDs)

#Extract all expression data for Norm/Tumour samples (immuno-suppressed excluded)
NT_exp_data <- exp_data.normalised_no_IS_merged_clean.all_532[,which(colnames(exp_data.normalised_no_IS_merged_clean.all_532) %in% sanofi_NT$V1)]

## Import top 20 genes from Caret script
top_20_genes_knnProfile <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/top_20_genes_knnProfile.txt", quote="", comment.char="")

#Extract expression data for these 20 genes only
NT_exp_data <- NT_exp_data[match(top_20_genes_knnProfile$V1, rownames(NT_exp_data)),]


o <- data.frame(rep(c("Norm_no_met", "Norm_with_met", "Tum_met", "Tum_no_met"), c(132, 71, 84, 146)))
colnames(o) <- "Group"


linear_pred_all <- vector()
for(i in 1:length(colnames(NT_exp_data))) {
  linear_pred_all[i] = 0
}
for(i in 1:length(colnames(NT_exp_data))) {
  for(j in 1:length(rownames(NT_exp_data))) {
    linear_pred_all[i] <- linear_pred_all[i] + tables_top_genes[j,2] * NT_exp_data[j,i]
  }
}


library(ggpubr)

#For comparisons on boxplots
cmpr <- list(c("Norm_with_met","Tum_no_met"), c("Norm_no_met","Tum_no_met"), c("Tum_no_met", "Tum_met"))


##### ggplot2 boxplot of linear predictor (20 genes) across all samples
temp_data <- data.frame(linear_pred=linear_pred_all, group=o$Group)

p <- ggplot(temp_data, aes(reorder(group, linear_pred), y=linear_pred)) + geom_boxplot()
p <- ggplot(temp_data, aes(reorder(group, linear_pred), y=linear_pred)) + geom_boxplot(width=0.6) + ylim(NA, 200)
p + geom_jitter(shape=16, position=position_jitter(0.2))
p + stat_compare_means(comparisons = cmpr, tip.length=0.01,label = "p.format") + geom_jitter(shape=16, position=position_jitter(0.2)) + labs(title="Linear predictors of disease progression across all samples",x="Group", y = "Linear predictors")





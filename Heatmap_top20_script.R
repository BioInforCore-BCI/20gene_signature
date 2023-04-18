## Script to produce heatmap of top 20 genes from KNN model across all primary tumour samples
## Supplemental fig 2


library(pheatmap)
library(methods)

#Read in top 20 genes from predicted model 
#Obtained from running script "Caret_Models_script.R"
top_20_genes_knnProfile <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/top_20_genes_knnProfile.txt", quote="", comment.char="")

#Load normalized expression matrix for all samples
exp_data.normalised_no_IS_merged_clean.all_532 <- read.delim("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/exp_data.normalised_no_IS_merged_clean.all_532.txt", row.names=1, quote="")

#Read in IDs for tumour samples
Tum_no_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Tum_no_met_IDs.txt", quote="\"", comment.char="")
Tum_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Tum_met_IDs.txt", quote="\"", comment.char="")
Tumour_IDs <- rbind(Tum_no_met_IDs, Tum_met_IDs)

#Extract expression data from tumour samples for top 20 genes only
top_20_KNN<-exp_data.normalised_no_IS_merged_clean.all_532[which(rownames(exp_data.normalised_no_IS_merged_clean.all_532) %in% top_20_genes_knnProfile$V1),]
top_20_KNN<-t(top_20_KNN[,which(colnames(top_20_KNN) %in% Tumour_IDs$V1)])

#Add column with sample type
group_tumour <- as.data.frame(rep(c("Tum_no_met", "Tum_met"), c(146, 84)))
colnames(group_tumour) <- "group"
rownames(group_tumour) <- rownames(top_20_KNN)

## Supplemental fig 2 heatmap
pheatmap(t(top_20_KNN), annotation_col = group_tumour, scale="row", show_colnames = F, clustering_method = "average")





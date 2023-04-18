library(ggfortify)
library(edgeR)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(msigdbr)

################### PREPARATION ###################

#Import all data
Tum_no_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Tum_no_met_IDs.txt", quote="\"", comment.char="")
Tum_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Tum_met_IDs.txt", quote="\"", comment.char="")
Norm_with_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Norm_with_met_IDs.txt", quote="\"", comment.char="")
Norm_of_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Norm_of_met_IDs.txt", quote="\"", comment.char="")
Norm_no_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Norm_no_met_IDs.txt", quote="\"", comment.char="")
Met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Met_IDs.txt", quote="\"", comment.char="")
exp_data.normalised_no_IS_merged_clean.all_532 <- read.delim("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/exp_data.normalised_no_IS_merged_clean.all_532.txt", row.names=1, quote="")


#List of all IDs
sanofi_all <- rbind(Met_IDs, Norm_no_met_IDs, Norm_of_met_IDs, Norm_with_met_IDs, Tum_met_IDs, Tum_no_met_IDs)

#Extract all expression data for all patients (immuno-suppressed excluded)
All_exp_data<-t(exp_data.normalised_no_IS_merged_clean.all_532[,which(colnames(exp_data.normalised_no_IS_merged_clean.all_532) %in% sanofi_all$V1)])

#Numbers per group for PCA labelling
z <- data.frame(rep(c("Met", "Norm_no_met", "Norm_of_met", "Norm_with_met", "Tum_met", "Tum_no_met"), c(73, 132, 7, 71, 84, 146)))
colnames(z) <- "Type"


###############################################################################

##PCA plot of all samples (not included in final paper)

#Separate numeric version of this matrix for compatibility with prcomp
All_exp_data_numeric <- All_exp_data

#Add labels for PCA
All_exp_data <- cbind(All_exp_data, z)

#Run PCA and plot with colours by type
pca_all <- prcomp(All_exp_data_numeric)
autoplot(pca_all, data = All_exp_data, colour = 'Type')



###############################################################################

### Tum-met versus Tum_no_met

################### Differential expression analysis ###################

#Extract all expression data for Tum_met and Tum_no_met patients
Tum_met<-exp_data.normalised_no_IS_merged_clean.all_532[,which(colnames(exp_data.normalised_no_IS_merged_clean.all_532) %in% Tum_met_IDs$V1)]
Tum_no_met<-exp_data.normalised_no_IS_merged_clean.all_532[,which(colnames(exp_data.normalised_no_IS_merged_clean.all_532) %in% Tum_no_met_IDs$V1)]

#Merge into one data frame
Tum_met_vs_tum_no_met <- cbind(Tum_met, Tum_no_met)

#Factor of sample types
factors_tum_met_vs_tum_no_met <- factor(rep(c("Tum_met", "Tum_no_met"), c(84,146)))

#Create design matrix
mm.tm_vs_tnm <- model.matrix(~factors_tum_met_vs_tum_no_met-1)
colnames(mm.tm_vs_tnm) <- levels(factors_tum_met_vs_tum_no_met)

#Set up contrasts for pairwise comparisons
tm_vs_tnm.contrast.mat <- makeContrasts(
  Diff = Tum_met - Tum_no_met,
  levels = colnames(mm.tm_vs_tnm))

#lmFit fits a linear model using weighted least squares for each gene
tm_vs_tnm_fit <- lmFit(Tum_met_vs_tum_no_met, mm.tm_vs_tnm)

#Estimate contrast for each gene
tm_vs_tnm_contrast <- contrasts.fit(tm_vs_tnm_fit, tm_vs_tnm.contrast.mat)

#Empirical Bayes smoothing of standard errors
tm_vs_tnm_fit <- eBayes(tm_vs_tnm_contrast)

#Show top DE genes
top.table <- topTable(tm_vs_tnm_fit, n = Inf, p.value = 0.05, sort.by = "P")

#Save limma output to table
write.table(top.table, file = "Limma_DE_genes_Tum_met_vs_Tum_no_met.tsv", quote = FALSE, sep = "\t")

#Table of all genes, including p>0.05
#Need all genes for GSEA later
top.table_unsig <- topTable(tm_vs_tnm_fit, n = Inf, sort.by = "P")

#Subset all DE genes with p value <0.05 and abs(log fold) > 1
DE_genes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > 1))

DE_tm_vs_tnm <- subset(Tum_met_vs_tum_no_met, rownames(Tum_met_vs_tum_no_met) %in% DE_genes)



##PCA plot (Not in final paper)

#Numbers per group for PCA labelling
y <- data.frame(rep(c("Tum_met", "Tum_no_met"), c(84, 146)))
colnames(y) <- "Type"

DE_tm_vs_tnm <- cbind(t(DE_tm_vs_tnm), y)

#PCA plot of DE genes
pca_tm_tnm <- prcomp(DE_tm_vs_tnm[1:229],)
autoplot(pca_tm_tnm, data = DE_tm_vs_tnm, colour = 'Type')

#################################################################




################### ClusterProfiler analysis ###################

######## Setup

organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#Convert gene symbols to entrez IDs
entrez <- bitr(DE_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb = organism)
entrez <- entrez$ENTREZID


######### Over representation analysis

### Reactome

reactomeEnrich <- enrichPathway(gene = entrez, pvalueCutoff = 0.05, readable = TRUE)

#Bar plot using qscore
mutate(reactomeEnrich, qscore = -log(p.adjust, base=10)) %>%
  arrange(desc(qscore)) %>%
  barplot(x="qscore")

#Add qscore to result table
test <- mutate(reactomeEnrich, qscore = -log(p.adjust, base=10))
reactomeEnrichTable <- test@result

##Manually produced supplemental fig 1A from this table in Excel
write.table(reactomeEnrichTable, file = "Over_representation_tumourvstumourmet.tsv", quote = FALSE, sep = "\t")




######## Cell markers

#Setup
#Manually downloaded these, should be the same as ones used on CLusterProfiler
Human_cell_markers <- read.delim("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Human_cell_markers.txt")


cells <- Human_cell_markers %>%
  dplyr::select(cellName, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', ')) %>%
  tidyr::unnest(cols = c(geneID))

cell_markers <- enricher(entrez, TERM2GENE = cells)

##Supplemental fig 1B
#Barplot of over-represented cellular signatures for DE genes (with qscore)
mutate(cell_markers, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")




######### GSEA analysis

##Generate ranked gene list

#Use all genes, not just significant
DE_genes <- top.table_unsig[3]
#Bring rownames (gene symbols) into a column
DE_genes <- cbind(DE_genes, rownames(DE_genes))
#Rename column to allow merge
colnames(DE_genes)[2] <- "SYMBOL"
#Convert gene symbols to entrez IDs
entrez_temp <- bitr(rownames(DE_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb = organism)
#Merge entrez output with original table
DE_genes <- merge(DE_genes, entrez_temp, by = "SYMBOL")
#Extract entrez id and t-statistic columns only
DE_genes <- DE_genes[,c(3,2)]
#Convert to list format
geneList <- DE_genes[,2]
names(geneList) <- as.character(DE_genes[,1])
geneList <- sort(geneList, decreasing = TRUE)


#Load MSigDB datasets, C2 = curated gene sets
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene)

#GSEA against MSigDB curated gene set
msigdb_gsea <- GSEA(geneList, TERM2GENE = m_df)

gsea_result <- msigdb_gsea@result
rownames(gsea_result) <- 1:nrow(gsea_result)

write.table(gsea_result, file = "Tumour_GSEA_results_Reactome.tsv", quote = FALSE, sep = "\t")


#Load MSigDB datasets, C2 = curated gene sets
m_df_KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)

#GSEA against MSigDB curated gene set
msigdb_gsea_KEGG <- GSEA(geneList, TERM2GENE = m_df_KEGG)

gsea_result_KEGG <- msigdb_gsea_KEGG@result
rownames(gsea_result_KEGG) <- 1:nrow(gsea_result_KEGG)

write.table(gsea_result_KEGG, file = "Tumour_GSEA_results_KEGG.tsv", quote = FALSE, sep = "\t")


### Manually produced figure 1 using outputs from KEGG and REACTOME GSEA analysis in Excel


#Plots of GSEA results

#Harcoded for Upregulated manually chosen plots
p1 <- gseaplot(msigdb_gsea, geneSetID = 85, by = "runningScore", title = gsea_result$Description[85])
p2 <- gseaplot(msigdb_gsea, geneSetID = 86, by = "runningScore", title = gsea_result$Description[86])
p3 <- gseaplot(msigdb_gsea, geneSetID = 89, by = "runningScore", title = gsea_result$Description[89])
p4 <- gseaplot(msigdb_gsea, geneSetID = 190, by = "runningScore", title = gsea_result$Description[190])
p5 <- gseaplot(msigdb_gsea, geneSetID = 154, by = "runningScore", title = gsea_result$Description[154])
p6 <- gseaplot(msigdb_gsea, geneSetID = 94, by = "runningScore", title = gsea_result$Description[94])

par(mfrow = c(2,3))

p1 + p2 + p3 + p4 + p5 + p6


#####

#Harcoded for Downregulated manually chosen plots
p1 <- gseaplot(msigdb_gsea, geneSetID = 84, by = "runningScore", title = gsea_result$Description[84])
p2 <- gseaplot(msigdb_gsea, geneSetID = 95, by = "runningScore", title = gsea_result$Description[95])
p3 <- gseaplot(msigdb_gsea, geneSetID = 91, by = "runningScore", title = gsea_result$Description[91])
p4 <- gseaplot(msigdb_gsea, geneSetID = 413, by = "runningScore", title = gsea_result$Description[413])
p5 <- gseaplot(msigdb_gsea, geneSetID = 194, by = "runningScore", title = gsea_result$Description[194])
p6 <- gseaplot(msigdb_gsea, geneSetID = 354, by = "runningScore", title = gsea_result$Description[354])

par(mfrow = c(2,3))

p1 + p2 + p3 + p4 + p5 + p6


#### Cell markers

#GSEA against cell markers
gsea_cell_markers <- GSEA(geneList, TERM2GENE = cells)

gsea_cell_markers_table <- gsea_cell_markers@result
rownames(gsea_cell_markers_table) <- 1:nrow(gsea_cell_markers_table)

#Plot of GSEA results
p1 <- gseaplot(gsea_cell_markers, geneSetID = 2, by = "runningScore", title = gsea_cell_markers_table$Description[2])
p2 <- gseaplot(gsea_cell_markers, geneSetID = 3, by = "runningScore", title = gsea_cell_markers_table$Description[3])
p3 <- gseaplot(gsea_cell_markers, geneSetID = 5, by = "runningScore", title = gsea_cell_markers_table$Description[5])
p4 <- gseaplot(gsea_cell_markers, geneSetID = 1, by = "runningScore", title = gsea_cell_markers_table$Description[1])
p5 <- gseaplot(gsea_cell_markers, geneSetID = 9, by = "runningScore", title = gsea_cell_markers_table$Description[9])
p6 <- gseaplot(gsea_cell_markers, geneSetID = 8, by = "runningScore", title = gsea_cell_markers_table$Description[8])

par(mfrow = c(2,3))

p1 + p2 + p3 + p4 + p5 + p6





###############################################################################


## Norm_with_met versus Norm_no_met

################### Differential expression analysis ###################

#Extract all expression data for Norm_met and Norm_no_met patients
Norm_met<-exp_data.normalised_no_IS_merged_clean.all_532[,which(colnames(exp_data.normalised_no_IS_merged_clean.all_532) %in% Norm_with_met_IDs$V1)]
Norm_no_met<-exp_data.normalised_no_IS_merged_clean.all_532[,which(colnames(exp_data.normalised_no_IS_merged_clean.all_532) %in% Norm_no_met_IDs$V1)]
Norm_of_met<-exp_data.normalised_no_IS_merged_clean.all_532[,which(colnames(exp_data.normalised_no_IS_merged_clean.all_532) %in% Norm_of_met_IDs$V1)]

#Merge into one data frame
Norm_met_vs_norm_no_met <- cbind(Norm_met, Norm_no_met)

#Factor of sample types
factors_norm_met_vs_norm_no_met <- factor(rep(c("Norm_met", "Norm_no_met"), c(71,132)))

#Create design matrix
mm.nm_vs_nnm <- model.matrix(~factors_norm_met_vs_norm_no_met-1)
colnames(mm.nm_vs_nnm) <- levels(factors_norm_met_vs_norm_no_met)


#Set up contrasts for pairwise comparisons
nm_vs_nnm.contrast.mat <- makeContrasts(
  Diff = Norm_met - Norm_no_met,
  levels = colnames(mm.nm_vs_nnm))


#lmFit fits a linear model using weighted least squares for each gene
nm_vs_nnm_fit <- lmFit(Norm_met_vs_norm_no_met, mm.nm_vs_nnm)

#Estimate contrast for each gene
nm_vs_nnm_contrast <- contrasts.fit(nm_vs_nnm_fit, nm_vs_nnm.contrast.mat)

#Empirical Bayes smoothing of standard errors
nm_vs_nnm_fit <- eBayes(nm_vs_nnm_contrast)

#Show top DE genes
top.table_norm <- topTable(nm_vs_nnm_fit, n = Inf, p.value = 0.05, sort.by = "P")

#Save limma output to table
write.table(top.table_norm, file = "Limma_DE_genes_Norm_met_vs_Norm_no_met.tsv", quote = FALSE, sep = "\t")

#Table of all genes, including p<0.05
top.table_norm_unsig <- topTable(nm_vs_nnm_fit, n = Inf, sort.by = "P")

#Subset all DE genes with p value <0.05 and abs(log fold) > 1
DE_genes_norm <- rownames(subset(top.table_norm, top.table_norm$adj.P.Val < 0.05 & abs(top.table_norm$logFC) > 1))

DE_nm_vs_nnm <- subset(Norm_met_vs_norm_no_met, rownames(Norm_met_vs_norm_no_met) %in% DE_genes_norm)


##PCA plot

#Bring in Norm_of_met
All_norm <- cbind(Norm_met, Norm_no_met, Norm_of_met)
DE_all_norm <- subset(All_norm, rownames(All_norm) %in% DE_genes_norm)


#Numbers per group for PCA labelling
w <- data.frame(rep(c("Norm_with_met", "Norm_no_met", "Norm_of_met"), c(71, 132, 7)))
colnames(w) <- "Type"

DE_all_norm <- cbind(t(DE_all_norm), w)

#PCA plot of DE genes
pca_norm <- prcomp(DE_all_norm[1:198],)
autoplot(pca_norm, data = DE_all_norm, colour = 'Type')

#################################################################




################### ClusterProfiler analysis ###################

## Results from this section were used to produce supplemental Table II


organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#Convert gene symbols to entrez IDs
entrez_norm <- bitr(DE_genes_norm, fromType="SYMBOL", toType="ENTREZID", OrgDb = organism)
entrez_norm <- entrez_norm$ENTREZID


######### Over representation analysis

### Reactome

reactomeEnrich_norm <- enrichPathway(gene = entrez_norm, pvalueCutoff = 0.05, readable = TRUE)


#Bar plot using qscore
mutate(reactomeEnrich_norm, qscore = -log(p.adjust, base=10)) %>%
  arrange(desc(qscore)) %>%
  barplot(x="qscore")

#Add qscore to result table
test <- mutate(reactomeEnrich_norm, qscore = -log(p.adjust, base=10))
reactomeEnrichTable_norm <- test@result

write.table(reactomeEnrichTable_norm, file = "reactomEnrichTable_norm.tsv", sep = "\t", quote = FALSE)

######### GSEA analysis

#Use all genes, not just significant
DE_genes <- top.table_norm_unsig[3]
#Bring rownames (gene symbols) into a column
DE_genes <- cbind(DE_genes, rownames(DE_genes))
#Rename column to allow merge
colnames(DE_genes)[2] <- "SYMBOL"
#Convert gene symbols to entrez IDs
entrez_temp <- bitr(rownames(DE_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb = organism)
#Merge entrez output with original table
DE_genes <- merge(DE_genes, entrez_temp, by = "SYMBOL")
#Extract entrez id and t-statistic columns only
DE_genes <- DE_genes[,c(3,2)]
#Convert to list format
geneList_norm <- DE_genes[,2]
names(geneList_norm) <- as.character(DE_genes[,1])
geneList_norm <- sort(geneList_norm, decreasing = TRUE)


#Load MSigDB datasets, C2 = curated gene sets
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene)

#GSEA against MSigDB curated gene set
msigdb_gsea_norm <- GSEA(geneList_norm, TERM2GENE = m_df)

gsea_result_norm <- msigdb_gsea_norm@result
rownames(gsea_result_norm) <- 1:nrow(gsea_result_norm)

write.table(gsea_result_norm, file = "gsea_result_norm_REACTOME.tsv", quote = FALSE, sep = "\t")


#Load MSigDB datasets, C2 = curated gene sets
m_df_KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)

#GSEA against MSigDB curated gene set
msigdb_gsea_norm <- GSEA(geneList_norm, TERM2GENE = m_df_KEGG)

gsea_result_norm <- msigdb_gsea_norm@result
rownames(gsea_result_norm) <- 1:nrow(gsea_result_norm)

write.table(gsea_result_norm, file = "gsea_result_norm_KEGG.tsv", quote = FALSE, sep = "\t")


#Plots of GSEA results

#Harcoded for Upregulated manually chosen plots
p1 <- gseaplot(msigdb_gsea_norm, geneSetID = 38, by = "runningScore", title = gsea_result_norm$Description[38])
p2 <- gseaplot(msigdb_gsea_norm, geneSetID = 40, by = "runningScore", title = gsea_result_norm$Description[40])
p3 <- gseaplot(msigdb_gsea_norm, geneSetID = 43, by = "runningScore", title = gsea_result_norm$Description[43])
p4 <- gseaplot(msigdb_gsea_norm, geneSetID = 50, by = "runningScore", title = gsea_result_norm$Description[50])
p5 <- gseaplot(msigdb_gsea_norm, geneSetID = 51, by = "runningScore", title = gsea_result_norm$Description[51])
p6 <- gseaplot(msigdb_gsea_norm, geneSetID = 41, by = "runningScore", title = gsea_result_norm$Description[41])
p7 <- gseaplot(msigdb_gsea_norm, geneSetID = 181, by = "runningScore", title = gsea_result_norm$Description[181])


par(mfrow = c(2,3))

p1 + p2 + p3 + p4 + p5 + p6

p7






##Script for testing models using Caret (on tumour samples with and without met)

######## Set up for Caret model training 

library(caret)

set.seed(998)

#Load in limma output table for tumour with met vs tumour without met
#Outputted from "limma_DE_analysis.R"
Limma_DE_genes_Tum_met_vs_Tum_no_met <- read.delim("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Limma_DE_genes_Tum_met_vs_Tum_no_met.tsv", quote="")

#Subset all DE genes with p value <0.05 and abs(log fold) > 1
DE_genes <- subset(Limma_DE_genes_Tum_met_vs_Tum_no_met, Limma_DE_genes_Tum_met_vs_Tum_no_met$adj.P.Val < 0.05 & abs(Limma_DE_genes_Tum_met_vs_Tum_no_met$logFC) > 1)
DE_genes <- as.data.frame(DE_genes$Gene)

#Load in normalized expression matrix for all samples and extract tumour samples
exp_data.normalised_no_IS_merged_clean.all_532 <- read.delim("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/exp_data.normalised_no_IS_merged_clean.all_532.txt", row.names=1, quote="")
Tum_no_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Tum_no_met_IDs.txt", quote="\"", comment.char="")
Tum_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Tum_met_IDs.txt", quote="\"", comment.char="")
Tumour_IDs <- rbind(Tum_no_met_IDs, Tum_met_IDs)
Tumour_exp_data<-exp_data.normalised_no_IS_merged_clean.all_532[which(rownames(exp_data.normalised_no_IS_merged_clean.all_532) %in% DE_genes$`DE_genes$Gene`),]
Tumour_exp_data<-t(Tumour_exp_data[,which(colnames(Tumour_exp_data) %in% Tumour_IDs$V1)])


#Set up for tumours with and without met 
group_tumour <- as.data.frame(rep(c("N", "Y"), c(146, 84)))
colnames(group_tumour) <- "group"
#Add this data to tumour data frame
Tumour_exp_data <- cbind(Tumour_exp_data, group_tumour)


####
#Avoid rerunning below line 
#This will randomly create new training/testing sets
inTraining <- createDataPartition(Tumour_exp_data$group, p = .75, list = FALSE)
####


#Extract testing and training data sets for Caret
training <- Tumour_exp_data[inTraining,]
testing  <- Tumour_exp_data[-inTraining,]


x <- training[,c(1:1274)]
x <- as.data.frame(x)
y <- as.vector(training$group)
y <- factor(y, levels=c('N', 'Y'))





###############################################################################

#SVM ML

ctrl <- rfeControl(functions = caretFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

svmProfile <- rfe(x, y,
                  sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),
                  rfeControl = ctrl,
                  ## pass options to train()
                  method = "svmRadial")



##### Use top genes in the training and testing data to perform SVM 

## Tested a variety of numbers of top genes, but not shown here to reduce script length

top_genes <- svmProfile$fit$ptype[1:40]
top_genes <- t(top_genes)


training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[1275,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[41] <- "group"
training_t_sel <- as.data.frame(training_t_sel)
training_top_genes <- t(training_t_sel)
training_top_genes <- as.data.frame(training_top_genes)
training_top_genes1 <- as.data.frame(apply(training_top_genes, 2, as.numeric))
training_top_genes1$group <- training_top_genes$group
training_top_genes <- training_top_genes1
rm(training_top_genes1)

testing_t <- t(testing)
testing_t_sel <- testing_t[match(top_genes$V1, rownames(testing_t)),]
head(testing_t_sel)
testing_t_sel <- rbind(testing_t_sel, testing_t[1275,])
rownames(testing_t_sel)
rownames(testing_t_sel)[41] <- "group"
testing_t_sel <- as.data.frame(testing_t_sel)
testing_top_genes <- t(testing_t_sel)
testing_top_genes <- as.data.frame(testing_top_genes)
testing_top_genes1 <- as.data.frame(apply(testing_top_genes, 2, as.numeric))
testing_top_genes1$group <- testing_top_genes$group
testing_top_genes <- testing_top_genes1
rm(testing_top_genes1)


##### perform a SVM Radial now based on the selected training matrix
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary,
                           savePredictions = TRUE)

set.seed(500)

svmFit1 <- train(group ~ ., data = training_top_genes, 
                 method = "svmRadial", 
                 trControl = fitControl, 
                 preProc = c("center", "scale"),
                 tuneLength = 8,
                 metric = "ROC")



predict(svmFit1, newdata = testing_top_genes)
test_prediction <- predict(svmFit1, newdata = testing_top_genes)

test_set <- data.frame(obs=factor(testing_top_genes$group), 
                       N=predict(svmFit1, newdata = testing_top_genes, type = "prob")[,1], 
                       Y=predict(svmFit1, newdata = testing_top_genes, type = "prob")[,2], 
                       pred=predict(svmFit1, testing_top_genes))


confusionMatrix(data = test_set$pred, reference = test_set$obs)
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")



##############################################################################

## Random Forest ML


set.seed(825)
rfProfile <- rfe(x, y, sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),rfeControl = rfeControl(functions = rfFuncs))



##### Use top genes in the training and testing data to perform random forest

## Tested a variety of numbers of top genes, but not shown here to reduce script length
top_genes <- rfProfile$fit$forest$ncat[1:60]
top_genes$V1 <- rownames(top_genes)


training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[1275,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[61] <- "group"
training_t_sel <- as.data.frame(training_t_sel)
training_top_genes <- t(training_t_sel)
training_top_genes <- as.data.frame(training_top_genes)
training_top_genes1 <- as.data.frame(apply(training_top_genes, 2, as.numeric))
training_top_genes1$group <- training_top_genes$group
training_top_genes <- training_top_genes1
rm(training_top_genes1)

testing_t <- t(testing)
testing_t_sel <- testing_t[match(top_genes$V1, rownames(testing_t)),]
head(testing_t_sel)
testing_t_sel <- rbind(testing_t_sel, testing_t[1275,])
rownames(testing_t_sel)
rownames(testing_t_sel)[61] <- "group"
testing_t_sel <- as.data.frame(testing_t_sel)
testing_top_genes <- t(testing_t_sel)
testing_top_genes <- as.data.frame(testing_top_genes)
testing_top_genes1 <- as.data.frame(apply(testing_top_genes, 2, as.numeric))
testing_top_genes1$group <- testing_top_genes$group
testing_top_genes <- testing_top_genes1
rm(testing_top_genes1)



##### perform a random forest now based on the selected training matrix
set.seed(825)
rfFit1 <- train(group ~ ., data = training_top_genes, 
                method = "rf", 
                trControl = fitControl, 
                tuneLength = 8,
                metric = "ROC")

predict(rfFit1, newdata = testing_top_genes)


test_set <- data.frame(obs=factor(testing_top_genes$group), N=predict(rfFit1, newdata = testing_top_genes, type = "prob")[,1], 
                       Y=predict(rfFit1, newdata = testing_top_genes, type = "prob")[,2], pred=predict(rfFit1, testing_top_genes))


confusionMatrix(data = test_set$pred, reference = test_set$obs)
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")


##############################################################################

## Lda ML


set.seed(825)
ldaProfile <- rfe(x, y,
                  sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),
                  rfeControl = rfeControl(functions = ldaFuncs, method = "repeatedcv", repeats = 5, verbose = FALSE))


##### Use top genes in the training and testing data to perform Lda 

## Tested a variety of numbers of top genes, but not shown here to reduce script length

top_genes <- rownames(ldaProfile$fit$scaling)[1:10]
top_genes$V1 <- top_genes$x

training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[1275,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[11] <- "group"
training_t_sel <- as.data.frame(training_t_sel)
training_top_genes <- t(training_t_sel)
training_top_genes <- as.data.frame(training_top_genes)
training_top_genes1 <- as.data.frame(apply(training_top_genes, 2, as.numeric))
training_top_genes1$group <- training_top_genes$group
training_top_genes <- training_top_genes1
rm(training_top_genes1)

testing_t <- t(testing)
testing_t_sel <- testing_t[match(top_genes$V1, rownames(testing_t)),]
head(testing_t_sel)
testing_t_sel <- rbind(testing_t_sel, testing_t[1275,])
rownames(testing_t_sel)
rownames(testing_t_sel)[11] <- "group"
testing_t_sel <- as.data.frame(testing_t_sel)
testing_top_genes <- t(testing_t_sel)
testing_top_genes <- as.data.frame(testing_top_genes)
testing_top_genes1 <- as.data.frame(apply(testing_top_genes, 2, as.numeric))
testing_top_genes1$group <- testing_top_genes$group
testing_top_genes <- testing_top_genes1
rm(testing_top_genes1)


##### perform a lda now based on selected training matrix
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

set.seed(825)
ldaFit1 <- train(group ~ ., data = training_top_genes, 
                 method = "lda", 
                 trControl = fitControl, 
                 tuneLength = 8,
                 metric = "ROC")


predict(ldaFit1, newdata = testing_top_genes)


###### for lda
test_set <- data.frame(obs=factor(testing_top_genes$group), N=predict(ldaFit1, newdata = testing_top_genes, type = "prob")[,1], Y=predict(ldaFit1, newdata = testing_top_genes, type = "prob")[,2], pred=predict(ldaFit1, testing_top_genes))

confusionMatrix(data = test_set$pred, reference = test_set$obs)
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")


##############################################################################

## PAM ML

set.seed(825)

ctrl <- rfeControl(functions = caretFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

pamProfile <- rfe(x, y,
                  sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),
                  rfeControl = ctrl,
                  ## pass options to train()
                  method = "pam")


##### Use top genes in the training and testing data to perform PAM 

## Tested a variety of numbers of top genes, but not shown here to reduce script length
top_genes <- pamProfile$fit$ptype[1:30]
top_genes <- t(top_genes)


training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[1275,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[31] <- "group"
training_t_sel <- as.data.frame(training_t_sel)
training_top_genes <- t(training_t_sel)
training_top_genes <- as.data.frame(training_top_genes)
training_top_genes1 <- as.data.frame(apply(training_top_genes, 2, as.numeric))
training_top_genes1$group <- training_top_genes$group
training_top_genes <- training_top_genes1
rm(training_top_genes1)

testing_t <- t(testing)
testing_t_sel <- testing_t[match(top_genes$V1, rownames(testing_t)),]
head(testing_t_sel)
testing_t_sel <- rbind(testing_t_sel, testing_t[1275,])
rownames(testing_t_sel)
rownames(testing_t_sel)[31] <- "group"
testing_t_sel <- as.data.frame(testing_t_sel)
testing_top_genes <- t(testing_t_sel)
testing_top_genes <- as.data.frame(testing_top_genes)
testing_top_genes1 <- as.data.frame(apply(testing_top_genes, 2, as.numeric))
testing_top_genes1$group <- testing_top_genes$group
testing_top_genes <- testing_top_genes1
rm(testing_top_genes1)


##### perform a PAM now based on selected training matrix
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

set.seed(825)
pamFit1 <- train(group ~ ., data = training_top_genes, 
                 method = "pam", 
                 trControl = fitControl, 
                 tuneLength = 8,
                 metric = "ROC")

predict(pamFit1, newdata = testing_top_genes)


##### for PAM
test_set <- data.frame(obs=factor(testing_top_genes$group), N=predict(pamFit1, newdata = testing_top_genes, type = "prob")[,1], Y=predict(pamFit1, newdata = testing_top_genes, type = "prob")[,2], pred=predict(pamFit1, testing_top_genes))

confusionMatrix(data = test_set$pred, reference = test_set$obs)
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")


##############################################################################

#K-nearest neighbours

set.seed(825)
library(nnet)

ctrl <- rfeControl(functions = caretFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

knnProfile <- rfe(x, y,
                  sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),
                  rfeControl = ctrl,
                  ## pass options to train()
                  method = "knn")


##### Use top genes in the training and testing data to perform KNN

## Tested a variety of numbers of top genes, but not shown here to reduce script length

top_genes <- knnProfile$fit$ptype[1:20]
top_genes <- t(top_genes)


training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[1275,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[21] <- "group"
training_t_sel <- as.data.frame(training_t_sel)
training_top_genes <- t(training_t_sel)
training_top_genes <- as.data.frame(training_top_genes)
training_top_genes1 <- as.data.frame(apply(training_top_genes, 2, as.numeric))
training_top_genes1$group <- training_top_genes$group
training_top_genes <- training_top_genes1
rm(training_top_genes1)

testing_t <- t(testing)
testing_t_sel <- testing_t[match(top_genes$V1, rownames(testing_t)),]
head(testing_t_sel)
testing_t_sel <- rbind(testing_t_sel, testing_t[1275,])
rownames(testing_t_sel)
rownames(testing_t_sel)[21] <- "group"
testing_t_sel <- as.data.frame(testing_t_sel)
testing_top_genes <- t(testing_t_sel)
testing_top_genes <- as.data.frame(testing_top_genes)
testing_top_genes1 <- as.data.frame(apply(testing_top_genes, 2, as.numeric))
testing_top_genes1$group <- testing_top_genes$group
testing_top_genes <- testing_top_genes1
rm(testing_top_genes1)


##### perform a KNN now based on selected training matrix

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

set.seed(825)
knnFit1 <- train(group ~ ., data = training_top_genes, 
                 method = "knn", 
                 trControl = fitControl, 
                 tuneLength = 8,
                 metric = "ROC")

predict(knnFit1, newdata = testing_top_genes)


##### for KNN
test_set <- data.frame(obs=factor(testing_top_genes$group), N=predict(knnFit1, newdata = testing_top_genes, type = "prob")[,1], Y=predict(knnFit1, newdata = testing_top_genes, type = "prob")[,2], pred=predict(knnFit1, testing_top_genes))

confusionMatrix(data = test_set$pred, reference = test_set$obs)
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")



##############################################################################


#Calculate AUC for ground truth vs prediction

install.packages('cvAUC')
library(cvAUC)
library(pROC)

ground_truths <- as.vector(test_set$obs)
predictions <- as.vector(test_set$pred)

ground_truths <- gsub("Y", 1, ground_truths)
ground_truths <- gsub("N", 0, ground_truths)

predictions <- gsub("Y", 1, predictions)
predictions <- gsub("N", 0, predictions)

roc <- auc(as.numeric(ground_truths), as.numeric(predictions))
ci.auc(roc)
 

##############################################################################



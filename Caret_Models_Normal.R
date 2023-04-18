##Script for testing models using Caret (on normal samples with and without met)


###############################################################################

######## Set up for Caret model training 

library(caret)

#Load in output from limma analysis (limma_DE_analysis.R)
Limma_DE_genes_Norm <- read.delim("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Limma_DE_genes_Norm_met_vs_Norm_no_met.tsv", quote="")

#Subset all DE genes with p value <0.05 and abs(log fold) > 1
DE_genes <- subset(Limma_DE_genes_Norm, Limma_DE_genes_Norm$adj.P.Val < 0.05 & abs(Limma_DE_genes_Norm$logFC) > 1)
DE_genes <- as.data.frame(rownames(DE_genes))

#Load in normalized expression matrix for all samples, and subset normal samples only
exp_data.normalised_no_IS_merged_clean.all_532 <- read.delim("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/exp_data.normalised_no_IS_merged_clean.all_532.txt", row.names=1, quote="")
Norm_with_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Norm_with_met_IDs.txt", quote="\"", comment.char="")
Norm_no_met_IDs <- read.table("~/Desktop/bci/RNAseq_Sanofi/Transcriptome_reanalysis/Norm_no_met_IDs.txt", quote="\"", comment.char="")

Norm_IDs <- rbind(Norm_no_met_IDs, Norm_with_met_IDs)

Norm_exp_data<-exp_data.normalised_no_IS_merged_clean.all_532[which(rownames(exp_data.normalised_no_IS_merged_clean.all_532) %in% DE_genes$`rownames(DE_genes)`),]
Norm_exp_data<-t(Norm_exp_data[,which(colnames(Norm_exp_data) %in% Norm_IDs$V1)])

Norm_exp_data <- Norm_exp_data[order(match(rownames(Norm_exp_data), (as.vector(Norm_IDs$V1)))), , drop = FALSE]


#Labelling for with and without met
group_norm <- as.data.frame(rep(c("N", "Y"), c(132, 71)))
colnames(group_norm) <- "group"
#Add this data to normal samples data frame
Norm_exp_data <- cbind(Norm_exp_data, group_norm)

####
#Avoid rerunning below line 
#This will randomly create new training/testing sets
inTraining <- createDataPartition(Norm_exp_data$group, p = .75, list = FALSE)
####

#Extract testing and training data sets for Caret
training <- Norm_exp_data[inTraining,]
testing  <- Norm_exp_data[-inTraining,]


x <- training[,c(1:198)]
x <- as.data.frame(x)
y <- as.vector(training$group)
y <- factor(y, levels=c('N', 'Y'))


###############################################################################


###############################################################################

#K-nearest neighbours

set.seed(500)

ctrl <- rfeControl(functions = caretFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

knnProfile <- rfe(x, y,
                  sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),
                  rfeControl = ctrl,
                  ## pass options to train()
                  method = "knn")

plot(knnProfile)


##### Use top genes in the training and testing data to perform KNN 

## Tested a variety of numbers of top genes, but not shown here to reduce script length

top_genes <- knnProfile$fit$ptype[1:70]
top_genes <- t(top_genes)

training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[199,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[71] <- "group"
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
testing_t_sel <- rbind(testing_t_sel, testing_t[199,])
rownames(testing_t_sel)
rownames(testing_t_sel)[71] <- "group"
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

set.seed(500)
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




###############################################################################

#SVM

ctrl <- rfeControl(functions = caretFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

svmProfile <- rfe(x, y,
                  sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),
                  rfeControl = ctrl,
                  ## pass options to train()
                  method = "svmRadial")

plot(svmProfile)


##### Use top genes in the training and testing data to perform SVM

## Tested a variety of numbers of top genes, but not shown here to reduce script length

top_genes <- svmProfile$fit$ptype[1:22]
top_genes <- t(top_genes)

training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[199,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[23] <- "group"
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
testing_t_sel <- rbind(testing_t_sel, testing_t[199,])
rownames(testing_t_sel)
rownames(testing_t_sel)[23] <- "group"
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


###############################################################################


## Random Forest


set.seed(825)
rfProfile <- rfe(x, y, sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),rfeControl = rfeControl(functions = rfFuncs))


##### Use top genes in the training and testing data to perform random forest 

## Tested a variety of numbers of top genes, but not shown here to reduce script length

top_genes <- rfProfile$fit$forest$ncat[1:40]
top_genes$V1 <- rownames(top_genes)


training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[199,])
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
testing_t_sel <- rbind(testing_t_sel, testing_t[199,])
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
set.seed(825)
rfFit1 <- train(group ~ ., data = training_top_genes, 
                method = "rf", 
                trControl = fitControl, 
                tuneLength = 8,
                metric = "ROC")

predict(rfFit1, newdata = testing_top_genes)

#### for random forest
test_set <- data.frame(obs=factor(testing_top_genes$group), N=predict(rfFit1, newdata = testing_top_genes, type = "prob")[,1], 
                       Y=predict(rfFit1, newdata = testing_top_genes, type = "prob")[,2], pred=predict(rfFit1, testing_top_genes))


confusionMatrix(data = test_set$pred, reference = test_set$obs)
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")


###############################################################################


## Lda


set.seed(825)
ldaProfile <- rfe(x, y,
                  sizes = c(5, 10, 12, 15, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100),
                  rfeControl = rfeControl(functions = ldaFuncs, method = "repeatedcv", repeats = 5, verbose = FALSE))

##### Use top genes in the training and testing data to perform lda 

## Tested a variety of numbers of top genes, but not shown here to reduce script length

top_genes <- rownames(ldaProfile$fit$scaling)[1:15]
write.table(top_genes, file="top_15_genes_NORM_ldaProfile.txt", quote=F, sep="\t")
top_genes <- read.table(file="top_15_genes_NORM_ldaProfile.txt", sep="\t")
top_genes$V1 <- top_genes$x

training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[199,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[16] <- "group"
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
testing_t_sel <- rbind(testing_t_sel, testing_t[199,])
rownames(testing_t_sel)
rownames(testing_t_sel)[16] <- "group"
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

test_set
confusionMatrix(data = test_set$pred, reference = test_set$obs)
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")


###############################################################################


## PAM

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

top_genes <- pamProfile$fit$ptype[1:80]
top_genes <- t(top_genes)

training_t <- t(training)
training_t_sel <- training_t[match(top_genes$V1, rownames(training_t)),]
head(training_t_sel)
training_t_sel <- rbind(training_t_sel, training_t[199,])
head(training_t_sel)
tail(training_t_sel)
rownames(training_t_sel)
rownames(training_t_sel)[81] <- "group"
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
testing_t_sel <- rbind(testing_t_sel, testing_t[199,])
rownames(testing_t_sel)
rownames(testing_t_sel)[81] <- "group"
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


###############################################################################


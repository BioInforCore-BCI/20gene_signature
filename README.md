# 20gene_signature


Scripts to accompany the following paper: 

>Transcriptomic analysis of cutaneous squamous cell carcinoma reveals a multi-gene prognostic signature associated with metastasis (Wang et al., 2023)

***

1. [DE_analysis.R](#de_analysisr)
1. [Caret_Models_Tumour.R](#caret_models_tumourr)
1. [Caret_Models_Normal.R](#caret_models_normalr)
1. [Linear_predictor_analysis.R](#linear_predictor_analysisr)
1. [Heatmap_top20_script.R](#heatmap_top20_scriptr)

***

### DE_analysis.R

- Comparing samples: tumour with met vs tumour without met
- Includes ClusterProfiler over representation analysis 
- Code for supplemental fig 1
- Output from this script used to produce figure 1 in main paper
- Comparing samples: normal with met vs normal without met
- Ouput from this script used to produce supplemental table II

***

### Caret_models_tumour.R

- Scripts for testing predictive models using the Caret package
- Only for looking at tumour with met vs tumour without met
- Includes SVM, Random forest, Lda, PAM, KNN
- Used to identify the top 20 gene signature

***

### Caret_models_normal.R

- Scripts for testing predictive models using the Caret package
- Only for looking at normal with met vs normal without met
- Includes SVM, Random forest, Lda, PAM, KNN

***

### Linear_predictor_analysis.R

- Script for linear predictor analysis of top 20 gene signature
- Code for supplmental fig 3 - boxplots of linear predictors for training vs testing sets
- Code for ROC curves - figure 3 in main paper
- Code for supplemental fig 4 - boxplots of linear predictors across all four clinical groups 

***

### Heatmap_top20_script.R 

- Script to generate heatmap for supplemental fig 2
- Based on the top 20 gene signature for tumour samples (met vs without met)




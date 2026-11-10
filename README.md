---
title: "Research Project: Machine Learning-Based Detection of Coronary Artery Disease (CAD) Using Gene Expression Data"
subtitle: "Dataset: GSE20680 | Platform: Affymetrix | Binary Classification: Control vs CAD"
author: "Prateek & Likhit"
date: "`r Sys.Date()`"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

# Introduction

This report presents a machine learning pipeline for classifying individuals as either
Control (healthy) or CAD (Coronary Artery Disease) using whole-blood gene expression data
from the publicly available GSE20680 dataset.

The analytical workflow consists of two main stages:

1. Bioinformatics preprocessing
2. Machine learning analysis 


The original dataset contains three clinical groups: Control (0), Intermediate (1), and
Case/CAD (2). Samples belonging to the Intermediate group are excluded to create a clearly
defined binary classification problem. The Intermediate group represents an uncertain,
transitional clinical phenotype; including it would introduce label ambiguity and increase
noise, so excluding it allows the models to learn more distinct expression patterns.

---

# Install Packages

```{r install-packages, eval=FALSE}


install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "vsn", "hgu133plus2.db",
                       "org.Hs.eg.db", "MKomics"))

install.packages(c("matrixStats", "ggplot2",
                   "mlr3", "mlr3learners", "mlr3verse",
                   "ranger", "e1071", "xgboost", "precrec"))
```

# Load Packages

```{r load-packages}
library(GEOquery)      
library(limma)         
library(vsn)           
library(MKomics)       
library(matrixStats)   
library(ggplot2)       

library(mlr3)          
library(mlr3learners)  
library(mlr3verse)     
```

---

# Load the Dataset

```{r load-data}
library(GEOquery)

gse <- getGEO(filename = "C:/Users/prate/Downloads/GSE20680_series_matrix.txt.gz",
              getGPL = FALSE)

expr_data  <- exprs(gse)     
pheno_data <- pData(gse)

dim(expr_data)
table(pheno_data$characteristics_ch1)
```

The dataset contains roughly 45,000 gene probes and 195 samples: 87 CAD cases, 52 healthy
controls, and 56 intermediate samples.

# Assign Group Labels (Control vs CAD)

This is the most important step: everything downstream depends on `group` being correct and
perfectly aligned with the columns of `expr_data`.

```{r assign-groups}
state <- trimws(sub(".*:", "", as.character(pheno_data$characteristics_ch1)))

group <- ifelse(state == "Control (0)", "Control",
         ifelse(state == "Case (2)",    "CAD", NA))

group <- factor(group, levels = c("Control", "CAD"))

length(group) == ncol(expr_data)

table(group, useNA = "ifany")
```

Expected: 52 Control, 87 CAD, and 56 Intermediate (NA).

# Removing Intermediate Samples

```{r drop-intermediate}
keep <- !is.na(group)

expr_data <- expr_data[, keep]        # 195 -> 139 columns
group     <- droplevels(group[keep])  

ncol(expr_data) == length(group)      
table(group)
```

This leaves 139 samples: 52 Controls and 87 CAD cases. The TRUE result confirms samples and
labels are aligned.

---

# Quality Control (Before Normalization)

## Boxplot of Raw Intensities

```{r qc-boxplot-before}
par(mar = c(5, 4, 2, 1))
boxplot(log2(expr_data + 1),
        las  = 2,
        cex.axis = 0.5,
        outline  = FALSE,
        main = "Raw log2 intensity distributions",
        ylab = "log2 intensity")
```

Each box is one sample. Differences in box position/spread indicate technical variation before
normalization.

## Density Plot of Raw Intensities

```{r qc-density-before}
E_raw <- log2(expr_data + 1)
plot(density(E_raw[, 1]),
     xlab = "log2 expression", ylab = "density",
     main = "Raw log2 intensity densities", type = "n")
apply(E_raw, 2, function(x) lines(density(x)))
```

Each curve is one sample. Curves that do not overlap suggest technical differences before
normalization.

## Correlation Plot BEFORE Normalization (MKomics)

```{r corplot-before}

cor_before <- cor(expr_data, method = "spearman")
dim(cor_before)

corPlot(cor_before,
        minCor  = min(cor_before),
        labels  = FALSE,
        title   = "Sample correlation BEFORE normalization (Spearman)")
```

Each cell shows the correlation between two samples. Before normalization, off-diagonal
correlations are more variable, reflecting technical differences between arrays.

## PCA Before Normalization

```{r pca-before}
expr_raw    <- log2(expr_data + 1)
vars_raw    <- rowVars(as.matrix(expr_raw))                     
top_var_raw <- names(sort(vars_raw, decreasing = TRUE))[1:2000] 

pca_before <- prcomp(t(expr_raw[top_var_raw, ]), scale. = TRUE)

# Proportion of variance explained by each PC.
var_before <- summary(pca_before)$importance[2, ]   # row 2 = proportion of variance
pc1_before <- round(var_before[1] * 100, 1)
pc2_before <- round(var_before[2] * 100, 1)
cat("Variance explained BEFORE normalization:\n")
cat("  PC1:", pc1_before, "%\n")
cat("  PC2:", pc2_before, "%\n")

pca_before_df <- data.frame(PC1 = pca_before$x[, 1],
                            PC2 = pca_before$x[, 2],
                            Group = group)

ggplot(pca_before_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(Control = "#56B4E9", CAD = "#E69F00")) +
  theme_bw() +
  labs(title = "PCA before normalization",
       x = paste0("PC1 (", pc1_before, "%)"),
       y = paste0("PC2 (", pc2_before, "%)"))
```

Each point is one patient. Separation between Control and CAD suggests differing expression
patterns; overlap indicates similar profiles.

---

# Normalization (vsn)

```{r normalization-vsn}
vsn_fit   <- vsn2(expr_data)         
expr_norm <- as.matrix(exprs(vsn_fit)) 
head(expr_norm[, 1:4])
```

The vsn method removes technical differences between samples while stabilizing the
mean-variance relationship, producing values on a generalized-log scale that are well suited
for downstream statistics and machine learning.


```{r normalization-Quantile}
expr_log        <- log2(expr_data + 1)
expr_quantile   <- normalizeBetweenArrays(expr_log, method = "quantile")
```


## meanSdPlot: Comparing vsn vs Quantile Normalization

The `meanSdPlot` shows each gene's standard deviation against its mean rank. If a
normalization has stabilized the variance, the red running-median line should be
approximately **flat** -- meaning variability no longer depends on expression level.

```{r meansdplot-vsn}

meanSdPlot(expr_norm)
```

```{r meansdplot-quantile}

meanSdPlot(expr_quantile)
```

 Comparing the two meanSdPlots, the
quantile normalization produced a flatter and lower red line across the full range of
expression, whereas vsn showed a pronounced hump at intermediate expression ranks (the
red line rising and then falling), indicating that variance was less well stabilized. The
quantile method also kept the data on the familiar log2 scale, while vsn produced large
negative generalized-log values driven by many near-background probes. For these reasons,
quantile normalization was selected as the preferred method for this dataset, and it is
used for all downstream steps (`expr_norm` below refers to the quantile-normalized data).

```{r choose-normalization}

expr_norm <- expr_quantile
```

A flatter red line in the quantile meanSdPlot confirms that variance is well stabilized
across the range of expression for the data used downstream.

## Boxplot After Normalization

```{r qc-boxplot-after}
par(mar = c(7, 4, 2, 1))
boxplot(expr_norm , 
        las = 2,
        cex.axis = 0.5,
        main = "Intensity after quantile normalization",
        ylab = "log2 intensity",
        outline = FALSE)
```

The boxes now sit at similar levels, showing technical differences have been reduced.

## Correlation Plot AFTER Normalization (MKomics)

```{r corplot-after}
cor_after <- cor(expr_norm, method = "spearman")

corPlot(cor_after, minCor = min(cor_after),, labels = FALSE,
        title = "Sample correlation AFTER normalization (Spearman)")
```

Compared with the "before" heatmap, correlations are more uniform, indicating that samples are
now more directly comparable.

## PCA After Normalization

```{r pca-after}
vars_norm    <- rowVars(as.matrix(expr_norm))
top_var_norm <- names(sort(vars_norm, decreasing = TRUE))[1:2000]
pca_after    <- prcomp(t(expr_norm[top_var_norm, ]), scale. = TRUE)

# Proportion of variance explained by each PC.
var_after <- summary(pca_after)$importance[2, ]
pc1_after <- round(var_after[1] * 100, 1)
pc2_after <- round(var_after[2] * 100, 1)
cat("Variance explained AFTER normalization:\n")
cat("  PC1:", pc1_after, "%\n")
cat("  PC2:", pc2_after, "%\n")

pca_after_df <- data.frame(PC1 = pca_after$x[, 1],
                           PC2 = pca_after$x[, 2],
                           Group = group)

ggplot(pca_after_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(Control = "#56B4E9", CAD = "#E69F00")) +
  theme_bw() +
  labs(title = "PCA after normalization",
       x = paste0("PC1 (", pc1_after, "%)"),
       y = paste0("PC2 (", pc2_after, "%)"))
```

Each point is one patient after normalization. The plot also helps spot any remaining unusual
samples or technical variation.

---

# Preprocessing: Filter Low-Variance Genes (top 25%)

```{r filter-genes}
expr_norm_matrix <- as.matrix(expr_norm)

vars <- rowVars(expr_norm_matrix)

expr_filtered <- expr_norm_matrix[vars > quantile(vars, 0.75), ]

cat("Genes before filtering:", nrow(expr_norm_matrix), "\n")
cat("Genes after filtering :", nrow(expr_filtered),    "\n")

ncol(expr_filtered) == length(group)
```

Filtering at the 75th percentile keeps only the top 25% most variable genes. This focuses the
analysis on the most informative features, speeds up model training, and can reduce
overfitting. The TRUE result confirms samples and labels remain aligned.

---

# Differential Expression Analysis (DEG)

DEG identifies the genes that differ most between Control and CAD; the top genes become the
features for the machine learning models (feature selection).

```{r deg-design}
Design <- model.matrix(~ 0 + group)
colnames(Design) <- levels(group)  
```

```{r deg-fit}
fit <- lmFit(expr_filtered, Design)

cont.matrix <- makeContrasts(CADvsControl = CAD - Control, levels = Design)
fit2 <- contrasts.fit(fit, cont.matrix)

res.limma <- eBayes(fit2)
```

```{r deg-table}
res.tab <- topTable(res.limma, number = Inf, adjust.method = "BH")
head(res.tab)
```

The table summarizes each gene's expression difference and statistical significance. Genes at
the top show the strongest evidence of differences between CAD and Control.

## Volcano Plot

```{r volcano}
volcanoplot(res.limma, main = "Volcano Plot: CAD vs Control")
abline(v = c(-1, 1), col = "red", lwd = 2)   # two-fold change lines
```

Genes on the right are higher in CAD, on the left higher in Control; genes higher up are more
significant. Most genes cluster around a log2 fold change of 0, indicating generally small
expression differences.

---

# Select Top Genes for Machine Learning

```{r select-top-genes}
top_genes <- rownames(res.tab)[1:100]

ml_data <- as.data.frame(t(expr_filtered[top_genes, ]))

ml_data$label <- group

colnames(ml_data) <- make.names(colnames(ml_data), unique = TRUE)

dim(ml_data)
table(ml_data$label)
```

The final dataset contains 139 patients and 101 columns (100 gene features + 1 label), with 52
Controls and 87 CAD cases.

---

# Machine Learning with mlr3

We now switch to the mlr3 framework. mlr3 organizes machine learning into reusable objects:
a Task (the data + what to predict), Learners (the models), a Resampling strategy
(how we validate), and Measures (how we score). This is cleaner and more reproducible than
calling each model function separately.

## Create the Task

```{r mlr3-task}
set.seed(42)

task <- as_task_classif(ml_data, target = "label", positive = "CAD")
task
```

The Task holds the 100 gene features and the CAD/Control target, with CAD set as the positive
class (the class we most want to detect).

## Class Balancing on the TRAINING data only


# The data are imbalanced (52 Control vs 87 CAD). We balance by OVERSAMPLING the
smaller class (Control) up to the size of the larger class (CAD), following the
professor's approach:  sample(seq_len(52), 87, replace = TRUE)

# Balancing is done INSIDE the cross-validation pipeline (via a pipeline operator),
NOT on the whole dataset beforehand. If we oversampled first and then split,
duplicated Control patients could land in both training and test folds -- that leaks information and gives falsely high scores. Doing it inside the pipeline
means only each fold's TRAINING part is balanced, while the held-out part staysuntouched and honest.




```{r mlr3-oversample-pipeline}
library(mlr3pipelines)

po_balance <- po("classbalancing",
                 id        = "oversample",
                 adjust    = "minor",       
                 reference = "major",       
                 shuffle   = FALSE,
                 ratio     = 87 / 52)        # 52 -> ~87 (balanced)
```

## Define the Learners (models)

```{r mlr3-learners}
lrn_rf  <- lrn("classif.ranger",  predict_type = "prob", num.trees = 500)

lrn_svm <- lrn("classif.svm",     predict_type = "prob",
               type = "C-classification", kernel = "radial")

lrn_log <- lrn("classif.log_reg", predict_type = "prob")

glrn_rf  <- as_learner(po_balance %>>% lrn_rf)
glrn_svm <- as_learner(po_balance %>>% lrn_svm)
glrn_log <- as_learner(po_balance %>>% lrn_log)

glrn_rf$id  <- "RandomForest"
glrn_svm$id <- "SVM"
glrn_log$id <- "LogReg"
```

##  XGBoost learner

```{r mlr3-xgboost}
lrn_xgb  <- lrn("classif.xgboost", predict_type = "prob", nrounds = 100)
glrn_xgb <- as_learner(po_balance %>>% lrn_xgb)
glrn_xgb$id <- "XGBoost"
```

## Resampling: 5-fold Cross-Validation, Repeated 20 Times

```{r mlr3-resampling}
set.seed(42)
resampling <- rsmp("repeated_cv", folds = 5, repeats = 20)
resampling
```

Repeated cross-validation (5 folds x 20 repeats = 100 evaluations) provides a much more reliable
estimate of performance than a single split, directly addressing the small-sample concern.

## Benchmark All Models

```{r mlr3-benchmark}
set.seed(42)
design <- benchmark_grid(
  tasks       = task,
  learners    = list(glrn_rf, glrn_svm, glrn_log, glrn_xgb),
  resamplings = resampling
)

bmr <- benchmark(design)
```

## Compare Performance (Accuracy and AUC)

```{r mlr3-results}
measures <- list(
  msr("classif.acc", id = "Accuracy"),
  msr("classif.auc", id = "AUC")
)

results <- bmr$aggregate(measures)
results[, .(learner_id, Accuracy, AUC)]
```

The table shows each model's mean accuracy and mean AUC, averaged over all 100 cross-validation
rounds. Because these are averaged over many splits, they are far more trustworthy than a single
split, and the model with the highest mean AUC is the most reliable classifier.

## Visualize the AUC Comparison

```{r mlr3-boxplot}
library(mlr3viz)
autoplot(bmr, measure = msr("classif.auc")) +
  labs(title = "AUC across repeated cross-validation", y = "AUC")
```

Each box summarizes one model's AUC over all 100 cross-validation rounds. A high box that is also
narrow indicates a model that is both accurate and stable.

---

---

# Export Output File

The key results are written to disk so they can be shared and inspected outside R, as requested.

```{r export-results}
# 1. Differential expression results (all genes, ranked).
write.csv(res.tab, "DEG_results_CADvsControl.csv", row.names = TRUE)

# 2. Model performance comparison (accuracy + AUC per model).
perf <- as.data.frame(results[, .(learner_id, Accuracy, AUC)])
write.csv(perf, "ML_model_performance.csv", row.names = FALSE)

# 3. The processed ML dataset (100 gene features + label) used for modelling.
write.csv(ml_data, "ML_input_top100genes.csv", row.names = TRUE)

cat("Output files written to:", getwd(), "\n")
list.files(pattern = "\\.csv$")
```

These CSV files (differential expression results, model performance, and the processed
modelling dataset) constitute the output of the analysis and can be opened in any spreadsheet
program.

# Conclusion

This study developed a machine learning pipeline to distinguish CAD patients from healthy
Controls using the GSE20680 dataset. Intermediate samples were removed to create a clear binary
comparison.

The data were quality-checked (boxplots, density plots, Spearman correlation plots, and PCA).
Two normalization methods (quantile and vsn) were compared using meanSdPlot; quantile
normalization gave the flatter mean-standard-deviation relationship and was selected for
downstream analysis. The data were then filtered to the top 25% most variable genes, and
differential expression analysis selected the top 100 genes as features.

Classification was performed within the mlr3 framework, using class balancing (oversampling the
minority Control class inside the training folds) and 5-fold cross-validation repeated 20 times.
Random Forest, SVM, Logistic Regression, and XGBoost were compared on averaged accuracy and AUC,
giving a robust estimate of performance. Across the 100 resampling iterations, SVM achieved the
best performance (mean AUC approximately 0.90, accuracy approximately 0.81), followed by Random
Forest (AUC approximately 0.83) and XGBoost (AUC approximately 0.82). Logistic regression
performed close to chance (AUC approximately 0.58), consistent with complete separation in
high-dimensional data.

# Limitations

1. The dataset contains only 139 samples, which limits the power of the analysis even with
   repeated cross-validation.
2. The classes were imbalanced (52 Control vs 87 CAD). This was addressed by oversampling the
   minority class within the training folds, but oversampling duplicates existing samples rather
   than adding genuinely new information.
3. The top 100 genes were selected on the full dataset before cross-validation. A fully rigorous
   design would embed feature selection inside each resampling fold to avoid any selection bias.
4. No independent external dataset (e.g. GSE20681) was used for validation, so generalization to
   other cohorts is untested.
5. The important features are probe IDs and require gene annotation and biological validation
   before being described as confirmed CAD biomarkers.
6. Clinical factors (age, sex, medication, smoking) were not included and may also influence
   gene expression and CAD status.

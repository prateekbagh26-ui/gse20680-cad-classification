# Machine learning detection of coronary artery disease (CAD) from gene expression

Binary classification of Control vs CAD using whole-blood gene expression from the public GEO dataset
**GSE20680** (Affymetrix). Intermediate samples are excluded to give a clearly defined two-class problem.

Analysis by Prateek & Likhit.

## What the pipeline does

1. Quality control before and after normalisation (boxplots, density plots, Spearman correlation, PCA).
2. Normalisation: vsn and quantile are compared with `meanSdPlot`; quantile is used downstream.
3. Keeps the top 25% most variable genes.
4. Differential expression with `limma` (CAD vs Control); top 100 genes become the features.
5. Machine learning in `mlr3`: Random Forest, SVM, logistic regression and XGBoost, compared by
   accuracy and AUC over 5-fold cross-validation repeated 20 times. Class balancing happens inside each
   training fold.

## Repository layout

```
analysis/FINAL_Research_Project_v3.Rmd   main analysis, with explanations (knit this)
analysis/run_analysis.R                  the same code as a plain R script (no knitting needed)
data/                                    put the GEO download here (not tracked by git)
README.md
```

## Getting the data

The dataset is not stored in this repository. Download the series matrix file from GEO:

1. Go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20680
2. Under "Download family", download **GSE20680_series_matrix.txt.gz**.
3. Save it in the `data/` folder of this project.

## Running it

1. Open `analysis/FINAL_Research_Project_v3.Rmd` in RStudio.
2. Install the packages. The first chunk (`install-packages`) is set to `eval = FALSE`, so run it
   manually once. Two packages are used but not in that chunk, so install them too:
   ```r
   install.packages("BiocManager")
   BiocManager::install(c("GEOquery", "limma", "vsn", "hgu133plus2.db", "org.Hs.eg.db", "MKomics"))
   install.packages(c("matrixStats", "ggplot2", "mlr3", "mlr3learners", "mlr3verse",
                      "ranger", "e1071", "xgboost", "precrec",
                      "mlr3pipelines", "mlr3viz", "data.table"))
   ```
3. Click **Knit**. The data path is relative to the `.Rmd` file, so keep the folder structure above.

Or run the plain script without knitting, from the `analysis/` folder:

```bash
Rscript run_analysis.R
```

`run_analysis.R` contains the same code as the `.Rmd`, without the explanations. It is generated from
the `.Rmd` code chunks. If you change the `.Rmd`, regenerate the script or make the same change in both.

The script writes three CSV files into the working directory: DEG results, model performance, and the
processed machine learning dataset.

## Results (as reported in the report)

Mean over 100 cross-validation rounds, reported in the Conclusion section of the report:

| Model | Mean AUC | Mean accuracy |
|---|---|---|
| SVM | ~0.90 | ~0.81 |
| Random Forest | ~0.83 | - |
| XGBoost | ~0.82 | - |
| Logistic regression | ~0.58 | - |

These values come from the original run of the analysis. Re-running the code should reproduce them, but
it has not been re-run as part of this repository setup. Check the knitted output against them.

## Limitations

See the Limitations section at the end of the report. The main ones:

- Only 139 samples.
- Oversampling duplicates existing Control samples rather than adding new information.
- The top 100 genes were selected on the full dataset before cross-validation, which can inflate scores.
  Feature selection should be moved inside each resampling fold in a more rigorous design.
- No external validation cohort was used.
- Features are probe IDs and need gene annotation and biological validation before they are described as
  biomarkers.

## Changes from the original script

- The original `getGEO()` call pointed to a file on the author's own computer. It now reads
  `../data/GSE20680_series_matrix.txt.gz`, relative to the `.Rmd` file. No other code was changed.

## Credits

Dataset: GEO accession GSE20680 (Affymetrix). Analysis: Prateek & Likhit.

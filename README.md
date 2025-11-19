# _Comprehensive_DGE_Analysis_and_Visualization_using_DESeq2
_Comprehensive_DGE_Analysis_and_Visualization_using_DESeq2
# 🏃 Mouse Pancreatic Islet DGE Analysis: Exercise vs. Sedentary (GSE227516)

This R script automates a comprehensive bioinformatics pipeline for **Differential Gene Expression (DGE)** analysis of RNA-Seq count data from the **GSE227516** study. The analysis specifically compares gene expression in **mouse pancreatic islets** exposed to **Exercise** conditions versus **Sedentary** controls.

The pipeline utilizes the gold-standard **`DESeq2`** package, performs necessary Quality Control (QC) checks, and generates a suite of diagnostic and results plots for a complete and reproducible workflow.

## 🚀 Key Features

* **DESeq2 Implementation:** Follows the standard workflow for DGE analysis of count data using the `DESeq2` package.
* **QC & Data Transformation:** Includes pre-filtering of low-count genes and applies **Regularized Log (rlog)** and **Variance Stabilizing Transformation (VST)** for robust QC plots.
* **Comprehensive Visualization:** Generates a total of **8 diagnostic and results plots**, including PCA, Heatmaps, MA Plots, and Volcano Plots.
* **Top Gene Analysis:** Identifies, filters, and generates heatmaps for the **Top 10 Up- and Down-Regulated Genes** (LFC > 1, padj < 0.05).
* **Reproducibility:** Saves final significant gene lists to a CSV and includes `sessionInfo()` for environment tracking.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **DGE** | `DESeq2` (Negative Binomial GLM) | Identifies genes significantly altered by exercise relative to sedentary controls. |
| **QC & Clustering** | PCA, Sample Distance Heatmap | Assesses sample similarity, identifies potential batch effects, and confirms group separation. |
| **Normalization** | **Median of Ratios** (within `DESeq2`) | Corrects for sequencing depth and RNA composition differences. |
| **Significance Cutoffs** | $\text{padj} < 0.05$ and $|\text{log2FC}| > 1$ | Filters for high-confidence Differentially Expressed Genes (DEGs). |

---

## 🛠️ Prerequisites and Setup

### 📦 Input Files

The script assumes that the necessary count data and metadata files have been downloaded and are available in the specified directory:
1.  **`GSE227516_counts.csv`**: The raw gene count matrix.
2.  **`sample_information.csv`**: The metadata file containing the `condition` variable (e.g., Exercise/Sedentary).

***Note:*** *You must update the file paths in steps 11 and 12 of the R script to point to the correct location of your downloaded files.*

### ⚙️ Execution

1.  **Download** the R script and the two required CSV data files.
2.  **Ensure R Packages are Installed:** The script includes checks and installations for: `DESeq2`, `RUVSeq`, `pheatmap`, `RColorBrewer`, `ggplot2`, and `ggrepel`.
3.  **Optional:** Modify the `OUTPUT_PATH` variable (Step 1) to your desired saving location.
4.  **Execute** the script in your R environment:
    ```R
    source("GSE227516_Comprehensive_DGE_Analysis_and_Visualization_using_DESeq2.R")
    ```

---

## 📁 Output Files (8 Plots + 1 CSV)

All output files are saved to the specified `OUTPUT_PATH` (default: `D:/DOWNLOADS/`).

### Statistical Results

| Filename | Type | Description |
| :--- | :--- | :--- |
| `GSE227516_significant_DE_genes.csv` | CSV | Final list of genes where $\text{padj} < 0.05$ and $|\text{log2FC}| > 1$. |

### Visualization and QC Plots

| Filename | Analysis Stage | Description |
| :--- | :--- | :--- |
| `GSE227516_QC_Scatter_Plots.png` | QC | Comparison of $\log_2(x+1)$, rlog, and VST transformations. |
| `GSE227516_QC_Histograms.png` | QC | Histograms of raw counts, rlog, and VST counts to assess transformation effectiveness. |
| `GSE227516_QC_Sample_Distance_Heatmap.png` | QC | **Sample-to-Sample Distance Heatmap** (rlog transformed data). |
| `GSE227516_QC_PCA_Plot.png` | QC | **Principal Component Analysis (PCA)** plot for visual clustering of samples. |
| `GSE227516_QC_Dispersion_Plot.png` | QC | **Dispersion Plot** showing estimated and fitted dispersion functions. |
| `GSE227516_Results_MA_Plot.png` | Results | **MA Plot** (Mean-Difference Plot) showing log2FC vs. Average Expression. |
| `GSE227516_Results_Volcano_Plot.png` | Results | **Volcano Plot** highlighting significant up- and down-regulated genes ($\text{padj} < 0.05$, $|\text{log2FC}| > 1$). |
| `GSE227516_Results_Heatmap_Top10_Up.png` | Results | **Heatmap** of the Top 10 Up-Regulated DEGs. |
| `GSE227516_Results_Heatmap_Top10_Down.png` | Results | **Heatmap** of the Top 10 Down-Regulated DEGs. |

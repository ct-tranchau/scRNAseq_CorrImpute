# Correlation Methods for Single-Cell RNA-Seq: A Study of Promoter-Reporter and Native Gene Expression

Single-cell RNA sequencing (scRNA-seq) enables the exploration of transcriptomic heterogeneity at the individual cell level, facilitating the identification of distinct cell types, rare populations, and dynamic gene expression patterns over time. This technique has advanced plant biology by enabling cross-species comparisons and uncovering gene regulatory networks. While gene co-expression networks are essential for understanding these networks, commonly used correlation methods have primarily been applied to microarray and bulk RNA-seq datasets, with limited validation in single-cell contexts. This study uses an Arabidopsis scRNA-seq dataset to evaluate the performance of various correlation methods in comparing promoter-reporter genes with their native counterparts, addressing challenges such as zero-inflation and demonstrating the benefits of pseudo-bulk analysis and imputation techniques for improving correlation accuracy.

## Pipeline Overview
The pipeline consists of five main steps:

0. **Rank gene correlation** (`00_Rank_genes_correlation.R`)
1. **Compute gene-gene correlation using different methods** (`10_Correlation_comparisons.R`)
2. **Perform Pseudo bulk correlation analysis across all clusters** (`20_PseudoBulk_correlation.R`)
3. **Impute missing values** \
   3.1 *Variational Autoencoder (VAE)* (`30_VAE_imputation.py`) \
   3.2 *Generative adversarial network (GAN)* (`30_GAN_imputation.py`) \
   3.3 *Autoencoder (AE)* (`30_AE_imputation.py`) 
4. **TF-Target Gene Correlations** (`40_TFtargetgene_correlation.py`)
5. **Visualization of TF-Target Gene Correlations** (`50_Visualization_TF_original_vs_imputed_data.R`)

---

## Installation
### Requirements
- R (for correlation analysis)
- Python 3.7+ (for deep learning-based imputation)
- Dependencies:
  - R packages: `tidyselect_1.2.1`,`dplyr_1.1.4`,  `Seurat_5.2.1`, `WGCNA_1.73`, `infotheo_1.2.0.1`, `CSCORE_0.0.0.9000`, `optparse_1.7.5`, `ggplot2_3.5.1`, `reshape2_1.4.4`, `rstatix_0.7.2`, `ggpubr_0.6.0`
  - Python packages: `torch_2.3.1`, `numpy_1.26.4`, `pandas_2.0.3`, `argparse_1.1`, `scipy_1.13.1`, `matplotlib_3.8.4`

### Setup
#### Install R dependencies
```r
install.packages(c("tidyselect", "dplyr", "Seurat", "WGCNA", "infotheo", "optparse", "ggplot2", "reshape2", "rstatix", "ggpubr"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("ChangSuBiostats/CS-CORE")
```

#### Install Python dependencies
```bash
conda create -n CorrImpute python=3.7
conda activate CorrImpute
pip install torch==2.3.1 numpy==1.26.4 pandas==2.0.3 argparse==1.1 scipy==1.13.1 matplotlib==3.8.4
```

## Usage

### **Step 0: Rank gene correlation (`00_Rank_genes_correlation.R`)**
This script calculates and ranks the correlation of all genes with a specified gene of interest in a Seurat object. It extracts gene expression data, computes correlation coefficients using the specified method (Pearson, Spearman, or Kendall), and outputs a ranked list of genes from highest to lowest correlation.

#### **Run the script**
```bash
Rscript 00_Rank_genes_correlation.R --input_path "path/to/dataset.rds" --gene "gene_name" --output_path "path/to/correlation_rank.csv"
```
#### Arguments
- `--input_path` → Path to the input data (**Seurat RDS file**)
- `--data_type` *(optional)* → Type of data to extract from the Seurat object (counts, data, or scale.data) (**default: `"counts"`**)
- `--correlation_method` *(optional)* → Correlation method to use (pearson, spearman, or kendall) (**default: `"pearson"`**)
- `--gene` → The gene of interest to compute correlation against.
- `--output_path` → Path to save the output correlation results.

#### Example with Custom Parameters
```bash
Rscript 00_Rank_genes_correlation.R --input_path "Seurat_Obj_WER_wlabels.rds" --data_type "counts" --correlation_method "pearson" --gene "gene:AT5G14750" --output_path "correlation_rank.csv"
```


### **Step 1: Compute Gene-Gene Correlation (`10_Correlation_comparisons.R`)**
This script calculates gene-gene correlations using various methods, including Pearson, Spearman, Kendall, Mutual Information, Euclidean Distance, Manhattan Distance, Chebyshev Distance, and CS-CORE. It supports multiple data types (raw read counts, normalized data, and scaled data) and can be applied to different subsets of the dataset: all cells in the sample, cells expressing both genes, or cells expressing either gene 1, gene 2, or both.

#### **Run the script**
```bash
Rscript 10_Correlation_comparisons.R --input_path "path/to/dataset.rds or /path/to/dataset.h5" --gene1 "gene name" --gene2 "gene name"
```

#### **Arguments**
- `--input_path` → Path to the input data (**Seurat RDS file** or **HDF5 file**)
- `--gene1` → The name of Gene 1
- `--gene2` → The name of Gene 2

#### **Output**
- `correlation_results.csv` → File containing the computed correlation values.

#### **Example with Custom Parameters**
```bash
Rscript 10_Correlation_comparisons.R --input_path "Seurat_Obj_WER_wlabels.rds" --gene1 "GFP" --gene2 "gene:AT5G14750"
```

### **Step 2: Perform Pseudobulk Correlation Analysis (`20_PseudoBulk_correlation.R`)**
This script aggregates single-cell data into pseudo-bulk profiles by computing the average gene expression across all cells in each cluster. It measures the correlation between the average expression levels of the selected genes across clusters.

#### **Run the script**
```bash
Rscript 20_PseudoBulk_correlation.R --input_path "path/to/dataset.rds" --gene1 "gene name" --gene2 "gene name"
```
#### **Arguments**
- `--input_path` → Path to the input Seurat object data (**Seurat RDS file**)
- `--gene1` → The name of Gene 1
- `--gene2` → The name of Gene 2

#### **Output**
- `average_expression.csv` → Contains the average expression values and confidence intervals for all cells in each cluster.
- `correlation_plot.pdf` → Visualization of the correlation between average gene expression across clusters.
- `average_gene_expression_across_clusters.pdf` → Plot displaying the average gene expression across different clusters.

#### **Example with Custom Parameters**
```bash
Rscript 20_PseudoBulk_correlation.R --input_path "Seurat_Obj_WER_wlabels.rds" --gene1 "GFP" --gene2 "gene:AT5G14750"
```
### **Step 3: Impute missing values of single cell RNA-seq data**

### *3.1 Impute Data Using Variational Autoencoder (`30_VAE_imputation.py`)*
This script uses a **Variational Autoencoder (VAE)** to impute missing values in single-cell RNA-seq data.

#### **Run the script**
```bash
python 30_VAE_imputation.py --input_path "path/to/readcount.csv" --output_path "path/to/output.csv" --loss_plot_path "path/to/loss.pdf"
```
#### **Arguments**
- `--input_path` → Path to the expression matrix containing read count data (**CSV format**)
- `--output_path` → Path to save the imputed dataset
- `--loss_plot_path` →  Path to save the training loss graph
- `--epochs` *(optional)* → Number of training epochs (**default: 50**)
- `--batch_size` *(optional)* → Batch size for training (**default: 64**)
- `--lr` *(optional)* → Learning rate (**default: 0.0001**)

#### **Example with Custom Parameters**
```bash
python 30_VAE_imputation.py --input_path "WER_count.csv" --output_path "WER_VAE_impute.csv" --loss_plot_path "WER_VAE_loss.pdf" --epochs 50 --batch_size 64 --lr 0.0001
```

### *3.2: Impute Data Using Generative Adversarial Networks (`30_GAN_imputation.py`)*
This script uses a **Generative Adversarial Networks (GAN)** to impute missing values in single-cell RNA-seq data.

#### **Run the script**
```bash
python 30_GAN_imputation.py --input_path "path/to/readcount.csv" --output_path "path/to/output.csv" --loss_plot_path "path/to/loss.pdf"
```
#### **Arguments**
- `--input_path` → Path to the expression matrix containing read count data (**CSV format**)
- `--output_path` → Path to save the imputed dataset
- `--loss_plot_path` →  Path to save the training loss graph
- `--epochs` *(optional)* → Number of training epochs (**default: 50**)
- `--batch_size` *(optional)* → Batch size for training (**default: 64**)
- `--lr_g` *(optional)* → Learning rate of Generator (**default: 0.0001**)
- `--lr_d` *(optional)* → Learning rate of Discriminator (**default: 0.00001**)

#### **Example with Custom Parameters**
```bash
python 30_GAN_imputation.py --input_path "WER_count.csv" --output_path "WER_GAN_impute.csv" --loss_plot_path "WER_GAN_loss.pdf" --epochs 50 --batch_size 64 --lr_g 0.0001 --lr_d 0.00001
```

### *3.3: Impute Data Using Autoencoder (`30_AE_imputation.py`)*
This script uses a **Autoencoder (AE)** to impute missing values in single-cell RNA-seq data.

#### **Run the script**
```bash
python 30_AE_imputation.py --input_path "path/to/readcount.csv" --output_path "path/to/output.csv" --loss_plot_path "path/to/loss.pdf"
```
#### **Arguments**
- `--input_path` → Path to the expression matrix containing read count data (**CSV format**)
- `--output_path` → Path to save the imputed dataset
- `--loss_plot_path` →  Path to save the training loss graph
- `--epochs` *(optional)* → Number of training epochs (**default: 50**)
- `--batch_size` *(optional)* → Batch size for training (**default: 64**)
- `--lr` *(optional)* → Learning rate (**default: 0.0001**)

#### **Example with Custom Parameters**
```bash
python 30_AE_imputation.py --input_path "WER_count.csv" --output_path "WER_AE_impute.csv" --loss_plot_path "WER_AE_loss.pdf" --epochs 50 --batch_size 64 --lr 0.0001
```

### **Step 4: TF-Target Gene Correlations (`40_TFtargetgene_correlation.py`)**
This script computes pairwise correlations among known target genes for each transcription factor (TF) in both the original and imputed datasets.

#### **Run the script**
```bash
python 40_TFtargetgene_correlation.py \
    --original_data "path/to/original_data_readcount.csv" \
    --imputed_data "path/to/imputed_data.csv" \
    --tf_target_data "path/to/tf_target_data.csv" \
    --threshold 10 \
    --correlation_method pearson \
    --output_file "path/to/output.csv"
```

#### **Arguments**
- `--original_data` → Path to the original read count data (**CSV format**)
- `--imputed_data` → Path to the imputed dataset generated in Step 3 or its output (**CSV format**)
- `--tf_target_data` → Path to the transcription factor (TF) target gene file (**CSV format**)
- `--threshold` → Maximum number of target genes regulated by a single transcription factor (TF)
- `--correlation_method` *(optional)* → Correlation method to use. Choices: `pearson`, `spearman`, `kendall`
- `--output_file` →  Path to save the computed correlation results (**CSV format**)

#### **Example with Custom Parameters**
```bash
python 40_TFtargetgene_correlation.py \
--original_data "WER_count.csv" \
--imputed_data "WER_VAE_impute.csv" \
--tf_target_data "dap_dhs_filtered_root.csv" \
--threshold 3000 \
--correlation_method pearson \
--output_file TF_correlation_expression_comparison.csv
```

### **Step 5: Visualization of TF-Target Gene Correlations (`50_Visualization_TF_original_vs_imputed_data.R`)**
This script compares the correlation of target genes for each transcription factor (TF) between the original and imputed data or generates plots based on the output from Step 4.

#### **Run the script**
```bash
Rscript 50_Visualization_TF_original_vs_imputed_data.R --input_path "path/to/TF_correlation.csv" 
```
#### **Arguments**
- `--input_path` → Path to the transcription factor correlation (**CSV format**)

#### **Output**
The output will be saved in the "Visualization" folder.
- `lineplot.pdf` → . The line plot displays the average correlation of target genes for each TF in the original and imputed data.
- `Expression_correlation.pdf` → The scatter plots with fitted lines illustrate the relationship between gene expression and correlation for each TF.
- `Boxplot_expression_Ori_VS_Imputed.pdf` → The boxplot compares the average target gene expression for each TF.
- `Boxplot_correlation_Ori_VS_Imputed.pdf` → The boxplot compares the average target gene correlation for each TF.

#### **Example with Custom Parameters**
```bash
Rscript 50_Visualization_TF_original_vs_imputed_data.R --input_path "TF_correlation_expression_comparison.csv" 
```


## **Output Files**

| Step  | Output File                | Description                          |
|-------|----------------------------|--------------------------------------|
| **Step0** | `correlation_rank.csv` | Rank genes by correlation     |
| **Step1** | `correlation_results.csv` | Gene-gene correlation matrix        |
| **Step2** | `average_expression.csv` | PseudoBulk correlation results      |
| **Step3** | `VAE_impute.csv` `GAN_impute.csv` `AE_impute.csv` | VAE-imputed, GAN-imputed, AE-imputed expression matrix       |
| **Step4** | `TF_correlation_expression_comparison.csv` | Average correlation of all target genes in each TF    |
| **Step5** | `Lineplot.pdf` `Expression_correlation.pdf` `Boxplot_expression_Ori_VS_Imputed.pdf` `Boxplot_correlation_Ori_VS_Imputed.pdf`| Compare the correlation of target genes for each TF in the original and imputed data  |

## **License**
This project is licensed under the **MIT License**.


## **Contact**
For questions or issues, please contact:

📧 tnchau@vt.edu

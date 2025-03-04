# Correlation Methods for Single-Cell RNA-Seq: A Study of Promoter-Reporter and Native Gene Expression

Single-cell RNA sequencing (scRNA-seq) enables the exploration of transcriptomic heterogeneity at the individual cell level, facilitating the identification of distinct cell types, rare populations, and dynamic gene expression patterns over time. This technique has advanced plant biology by enabling cross-species comparisons and uncovering gene regulatory networks. While gene co-expression networks are essential for understanding these networks, commonly used correlation methods have primarily been applied to microarray and bulk RNA-seq datasets, with limited validation in single-cell contexts. This study uses an Arabidopsis scRNA-seq dataset to evaluate the performance of various correlation methods in comparing promoter-reporter genes with their native counterparts, addressing challenges such as zero-inflation and demonstrating the benefits of pseudo-bulk analysis and imputation techniques for improving correlation accuracy.

## Pipeline Overview
The pipeline consists of four main steps:

1. **Compute gene-gene correlation using different methods** (`10_correlation_method.R`)
2. **Perform Pseudo bulk correlation analysis across all clusters** (`20_PseudoBulk_correlation.R`)
3. **Impute missing values using Variational Autoencoder (VAE)** (`30_VAE_imputation.py`)
4. **Impute missing values using Autoencoder (AE)** (`40_AE_imputation.py`)

---

## Installation
### Requirements
- R (for correlation analysis)
- Python 3.8+ (for deep learning-based imputation)
- Dependencies:
  - R packages: `tidyverse`,`dplyr`,  `Seurat`, `WGCNA`, `infotheo`, `CSCORE`, `hdf5r`, `optparse`
  - Python packages: `torch`, `numpy`, `pandas`, `argparse`

### Setup
#### Install R dependencies
```r
install.packages(c("tidyverse", "dplyr", "Seurat", "WGCNA", "infotheo", "CSCORE", "hdf5r", "optparse"))
```

#### Install Python dependencies
```bash
pip install torch numpy pandas argparse
```

## Usage

### **Step 1: Compute Gene-Gene Correlation (`10_correlation_method.R`)**
This script calculates gene-gene correlations using various methods, including Pearson, Spearman, Kendall, Mutual Information, Euclidean Distance, Manhattan Distance, Chebyshev Distance, and CS-CORE. It supports multiple data types (raw read counts, normalized data, and scaled data) and can be applied to different subsets of the dataset: all cells in the sample, cells expressing both genes, or cells expressing either gene 1, gene 2, or both.

### **Run the script**
```bash
Rscript 10_Correlation_comparisons.R --dataset_path "path/to/dataset.rds or /path/to/dataset.h5" --gene1 "gene name" --gene2 "gene name"
```

### **Arguments**
- `--dataset_path` → Path to the input data (**Seurat RDS file** or **HDF5 file**)
- `--gene1` → The name of Gene 1
- `--gene2` → The name of Gene 2

### **Output**
- `correlation_results.csv` → File containing the computed correlation values.


### **Step 2: Perform Pseudobulk Correlation Analysis (`20_PseudoBulk_correlation.R`)**
This script aggregates single-cell data into pseudo-bulk profiles by computing the average gene expression across all cells in each cluster. It measures the correlation between the average expression levels of the selected genes across clusters.

### **Run the script**
```bash
Rscript 20_PseudoBulk_correlation.R --SeuratObj_path "path/to/dataset.rds" --gene1 "gene name" --gene2 "gene name"
```
### **Arguments**
- `--SeuratObj_path` → Path to the input Seurat object data (**Seurat RDS file**)
- `--gene1` → The name of Gene 1
- `--gene2` → The name of Gene 2

### **Output**
- `average_expression.csv` → Contains the average expression values and confidence intervals for all cells in each cluster.
- `correlation_plot.pdf` → Visualization of the correlation between average gene expression across clusters.
- `average_gene_expression_across_clusters.pdf` → Plot displaying the average gene expression across different clusters.


### **Step 3: Impute Data Using Variational Autoencoder (`30_VAE_imputation.py`)**
This script uses a **Variational Autoencoder (VAE)** to impute missing values in single-cell RNA-seq data.

### **Run the script**
```bash
python 30_VAE_imputation.py --input_path "path/to/readcount.csv" --output_path "path/to/output.csv"
```
### **Arguments**
- `--input_path` → Path to the expression matrix containing read count data (**CSV format**)
- `--output_path` → Path to save the imputed dataset
- `--epochs` *(optional)* → Number of training epochs (**default: 50**)
- `--batch_size` *(optional)* → Batch size for training (**default: 64**)
- `--lr` *(optional)* → Learning rate (**default: 0.0001**)

### **Example with Custom Parameters**
```bash
python 30_VAE_imputation.py --input_path expression_matrix.csv --output_path VAE_impute.csv --epochs 50 --batch_size 64 --lr 0.0001
```


### **Step 4: Impute Data Using Autoencoder (`40_AE_imputation.py`)**
This script uses a **Autoencoder (AE)** to impute missing values in single-cell RNA-seq data.

### **Run the script**
```bash
python 40_AE_imputation.py --input_path "path/to/readcount.csv" --output_path "path/to/output.csv"
```
### **Arguments**
- `--input_path` → Path to the expression matrix containing read count data (**CSV format**)
- `--output_path` → Path to save the imputed dataset
- `--epochs` *(optional)* → Number of training epochs (**default: 50**)
- `--batch_size` *(optional)* → Batch size for training (**default: 64**)
- `--lr` *(optional)* → Learning rate (**default: 0.0001**)

### **Example with Custom Parameters**
```bash
python 40_AE_imputation.py --input_path expression_matrix.csv --output_path AE_impute.csv --epochs 50 --batch_size 64 --lr 0.0001
```


## **Output Files**

| Step  | Output File                | Description                          |
|-------|----------------------------|--------------------------------------|
| **Step 1** | `correlation_results.csv`  | Gene-gene correlation matrix        |
| **Step 2** | `average_expression.csv`  | PseudoBulk correlation results      |
| **Step 3** | `VAE_impute.csv`         | VAE-imputed expression matrix       |
| **Step 4** | `AE_impute.csv`          | AE-imputed expression matrix        |

## **License**
This project is licensed under the **MIT License**.


## **Contact**
For questions or issues, please contact:

📧 tnchau@vt.edu

#!/usr/bin/env Rscript

# Load necessary packages
library(optparse)
library(Seurat)
library(tidyverse)
library(hdf5r) # for loading the h5 file
library(CSCORE) # for cs-core
library(sctransform) # for sctransform
library(WGCNA) # for bicor
library(infotheo) # for Mutual information

###############################
# Define command-line arguments
option_list <- list(
  make_option(c("--dataset_path"), type = "character", default = NULL, help = "Path to dataset (Seurat object .rds)"),
  make_option(c("--gene1"), type = "character", default = NULL, help = "First gene name"),
  make_option(c("--gene2"), type = "character", default = NULL, help = "Second gene name")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$dataset_path) || is.null(opt$gene1) || is.null(opt$gene2)) {
  print_help(opt_parser)
  stop("Missing required arguments: --dataset_path, --gene1, --gene2")
}

# Assign arguments to variables
dataset_path <- opt$dataset_path
gene1 <- opt$gene1
gene2 <- opt$gene2
##################################

# Check file type and load dataset
cat("\nLoading dataset...\n")
if (grepl("\\.h5$", dataset_path)) {
  matrix <- Read10X_h5(dataset_path)
  Seu_obj <- CreateSeuratObject(matrix)
} else if (grepl("\\.rds$", dataset_path)) {
  data_file <- readRDS(dataset_path)
  Seu_obj <- CreateSeuratObject(data_file)
} else {
  stop("\nError: Dataset file must be .h5 or .rds\n")
}
cat("Dataset loaded successfully.\n")

### **Data Processing**
cat("\nPerforming quality control...\n")
Subset <- subset(Seu_obj, subset = nFeature_RNA > 500 & nCount_RNA > 500)

# Normalize and scale the data
cat("\nNormalizing and scaling data...\n")
Sample <- NormalizeData(Subset, verbose = FALSE)
Sample <- FindVariableFeatures(Sample, selection.method = "vst", verbose = FALSE)
Sample <- ScaleData(Sample, features = rownames(Sample))

### **Gene Presence Check**
all_genes <- rownames(Sample)
if (!(gene1 %in% all_genes)) stop(paste("\nError:", gene1, "not found in dataset.\n"))
if (!(gene2 %in% all_genes)) stop(paste("\nError:", gene2, "not found in dataset.\n"))

### **Correlation Calculation**
cat("\nRunning CS-CORE...\n")
CSCORE_result <- CSCORE(Sample, genes = c(VariableFeatures(Sample), gene1, gene2))
cs_core_est <- as.numeric(CSCORE_result$est[gene1, gene2])
cs_core_pvalue <- as.numeric(CSCORE_result$p_value[gene1, gene2])
# Apply SCTransform
#cat("\nRunning SCTransform...\n")
#sctransform <- SCTransform(Subset, residual.features = c(VariableFeatures(Sample), gene1, gene2))
#sctransform_correlation <- cor(as.numeric(sctransform@assays$SCT@scale.data[gene1, ]), as.numeric(sctransform@assays$SCT@scale.data[gene2, ]), method = "pearson")

### **Function to Compute Gene Relationships**
compute_gene_relationships <- function(dat, gene1, gene2) {
  correlations <- list(
    pearson  = cor.test(dat[gene1, ], dat[gene2, ], method = "pearson"),
    spearman = cor.test(dat[gene1, ], dat[gene2, ], method = "spearman"),
    kendall  = cor.test(dat[gene1, ], dat[gene2, ], method = "kendall")
  )
  
  # Extract correlation estimates and p-values
  cor_results <- map(correlations, ~ list(correlation = .x$estimate, p_value = .x$p.value))
  
  # Biweight midcorrelation (bicor)
  cor_results$bicor <- bicor(dat[gene1, ], dat[gene2, ])
  
  # Mutual Information
  mutual_info <- function(x, y) {
    x_discrete <- discretize(x, disc = "equalfreq", nbins = 5)
    y_discrete <- discretize(y, disc = "equalfreq", nbins = 5)
    mutinformation(x_discrete, y_discrete)
  }
  cor_results$mutual_information <- mutual_info(dat[gene1, ], dat[gene2, ])
  
  # Distance calculations
  diffs <- abs(dat[gene1, ] - dat[gene2, ])
  cor_results$distances <- list(
    euclidean  = sqrt(sum(diffs^2)),
    manhattan  = sum(diffs),
    chebyshev  = max(diffs)
  )
  
  return(cor_results)
}

### **Extract Data and Compute Relationships**
cat("\nExtracting data layers...\n")
data_layers <- list(
  count_data    = GetAssayData(Sample, layer = "counts"),
  normalize_data = GetAssayData(Sample, layer = "data"),
  scale_data    = GetAssayData(Sample, layer = "scale.data")
)

# Compute relationships for all cells across different data layers
cat("\nComputing correlations for all cells...\n")
results_all_cells <- lapply(data_layers, function(x) compute_gene_relationships(x, gene1, gene2))
cat(paste0("\nTotal number of cells in dataset: ", length(colnames(data_layers$count_data)), "\n"))

# Identify cells expressing gene1 and gene2
nonzero_gene1_cells <- colnames(data_layers$count_data)[data_layers$count_data[gene1, ] != 0]
cat(paste0("Number of cells expressing ", gene1, ": ", length(nonzero_gene1_cells), "\n"))
nonzero_gene2_cells <- colnames(data_layers$count_data)[data_layers$count_data[gene2, ] != 0]
cat(paste0("Number of cells expressing ", gene2, ": ", length(nonzero_gene2_cells), "\n"))

# Cells expressing both genes
cells_both <- intersect(nonzero_gene1_cells, nonzero_gene2_cells)
cat(paste0("Number of cells expressing both genes: ", length(cells_both), "\n"))
results_both_cells <- if (length(cells_both) > 0) lapply(data_layers, function(x) compute_gene_relationships(x[, cells_both], gene1, gene2)) else NULL

# Cells expressing either gene1, gene2, or both genes
cells_either <- union(nonzero_gene1_cells, nonzero_gene2_cells)
cat(paste0("Number of cells expressing either gene1, gene2, or both genes: ", length(cells_either), "\n"))
results_either_cells <- lapply(data_layers, function(x) compute_gene_relationships(x[, cells_either], gene1, gene2))

##########################################
### **Format Results as a Tidy Dataframe**
format_results <- function(results, category) {
  if (is.null(results)) return(NULL)
  
  # List to store formatted data
  formatted_list <- list()
  
  for (data_layer in names(results)) {
    res <- results[[data_layer]]
    
    # Ensure all required elements exist
    cor_metrics <- c("pearson", "spearman", "kendall")
    distances <- c("euclidean", "manhattan", "chebyshev")
    
    formatted_list[[data_layer]] <- data.frame(
      Category = category,
      Data_Type = data_layer,
      Metric = c(
        "Pearson", "Spearman", "Kendall", "Biweight midcorrelation",
        "Mutual Information", "Euclidean Distance", "Manhattan Distance", "Chebyshev Distance"
      ),
      Value = c(
        if ("pearson" %in% names(res)) res$pearson$correlation else NA,
        if ("spearman" %in% names(res)) res$spearman$correlation else NA,
        if ("kendall" %in% names(res)) res$kendall$correlation else NA,
        if ("bicor" %in% names(res)) res$bicor else NA,
        if ("mutual_information" %in% names(res)) res$mutual_information else NA,
        if ("distances" %in% names(res)) res$distances$euclidean else NA,
        if ("distances" %in% names(res)) res$distances$manhattan else NA,
        if ("distances" %in% names(res)) res$distances$chebyshev else NA
      ),
      P_Value = c(
        if ("pearson" %in% names(res)) res$pearson$p_value else NA,
        if ("spearman" %in% names(res)) res$spearman$p_value else NA,
        if ("kendall" %in% names(res)) res$kendall$p_value else NA,
        NA,  # Biweight correlation has no p-value
        NA,  # Mutual Information has no p-value
        NA,  # Distances have no p-values
        NA,  
        NA   
      )
    )
  }
  
  # Combine results for different data layers into one dataframe
  return(bind_rows(formatted_list))
}

# Format and combine all results
df_all_cells <- format_results(results_all_cells, "All Cells")
df_both_cells <- format_results(results_both_cells, "Cells expressing both genes")
df_either_cells <- format_results(results_either_cells, "Cells expressing both genes and either gene")

# Combine all into one dataframe
final_results_df <- bind_rows(df_all_cells, df_both_cells, df_either_cells)

# Add CS-CORE result
final_results_df <- final_results_df %>%
  add_row(Category = "CS-CORE", Data_Type = "CS-CORE", Metric = "CS-CORE_Value", Value = cs_core_est, P_Value = cs_core_pvalue)

# Save results to CSV
write.csv(final_results_df, "correlation_results.csv", row.names = FALSE)
cat("\nAnalysis complete. Results saved to 'correlation_results.csv'.\n")

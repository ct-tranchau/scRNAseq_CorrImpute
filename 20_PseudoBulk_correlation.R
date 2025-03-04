#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(dplyr)
library(Seurat)
library(tidyverse)

# Define command-line arguments
option_list <- list(
  make_option(c("--SeuratObj_path"), type = "character", default = NULL, help = "Path to Seurat object (.rds)"),
  make_option(c("--gene1"), type = "character", default = NULL, help = "First gene name"),
  make_option(c("--gene2"), type = "character", default = NULL, help = "Second gene name")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$SeuratObj_path) || is.null(opt$gene1) || is.null(opt$gene2)) {
  print_help(opt_parser)
  stop("Missing required arguments: --SeuratObj_path, --gene1, --gene2")
}

# Assign arguments to variables
SeuratObj_path <- opt$SeuratObj_path
gene1 <- opt$gene1
gene2 <- opt$gene2


# Load Seurat object
cat("Loading Seurat object...\n")
Obj <- readRDS(SeuratObj_path)

# Extract count data
counts <- GetAssayData(Obj, layer = "counts")

# Check if genes exist in the dataset
if (!(gene1 %in% rownames(counts))) stop(paste("Error:", gene1, "not found in dataset."))
if (!(gene2 %in% rownames(counts))) stop(paste("Error:", gene2, "not found in dataset."))

# Transpose and convert to a data frame
subset <- data.frame(t(counts[c(gene1, gene2), ]), check.names = FALSE)
subset$celltype <- Idents(Obj)

# Ensure gene1 and gene2 columns are numeric
subset[[gene1]] <- as.numeric(subset[[gene1]])
subset[[gene2]] <- as.numeric(subset[[gene2]])

# Compute mean and confidence intervals
average_expression <- subset %>%
  group_by(celltype) %>%
  summarise(
    gene1_mean = mean(.data[[gene1]], na.rm = TRUE),
    gene1_sd = sd(.data[[gene1]], na.rm = TRUE),
    gene1_CI_low = gene1_mean - 1.96 * gene1_sd / sqrt(n()),
    gene1_CI_high = gene1_mean + 1.96 * gene1_sd / sqrt(n()),
    
    gene2_mean = mean(.data[[gene2]], na.rm = TRUE),
    gene2_sd = sd(.data[[gene2]], na.rm = TRUE),
    gene2_CI_low = gene2_mean - 1.96 * gene2_sd / sqrt(n()),
    gene2_CI_high = gene2_mean + 1.96 * gene2_sd / sqrt(n())
  ) %>%
  select(-gene1_sd, -gene2_sd)  # Remove unnecessary columns

# Save the average expression data as a CSV file
write.csv(average_expression, "average_expression.csv", row.names = FALSE)
cat("Saved average expression data of each cluster to 'average_expression.csv'.\n")

dat <- data.frame(average_expression)
# Function to compute and print correlation results
compute_correlation <- function(data, method) {
  cor_test <- cor.test(data$gene1_mean, data$gene2_mean, method = method)
  cat(paste0(method, " correlation: ", cor_test$estimate, "\n"))
  cat(paste0("P-value: ", cor_test$p.value, "\n\n"))
}

# Compute Pearson, Spearman, and Kendall correlations
compute_correlation(dat, "pearson")
compute_correlation(dat, "spearman")
compute_correlation(dat, "kendall")

# Save Correlation Plot as PDF
pdf("correlation_plot.pdf", width = 6, height = 6)
ggplot(dat, aes(x = gene1_mean, y = gene2_mean)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Best fit line
  labs(
    title = paste("Correlation between", gene1, "and", gene2, "Mean Expression"),
    x = paste(gene1, "Mean Expression"),
    y = paste(gene2, "Mean Expression")
  ) +
  theme_minimal()
dev.off()
cat("Saved correlation plot to 'correlation_plot.pdf'.\n")

# Save Line Plot with Error Bars as PDF
pdf("average_gene_expression_across_clusters.pdf", width = 8, height = 6)
ggplot(average_expression, aes(x = celltype)) +
  # Plot Gene 1 Mean with Error Bars
  geom_line(aes(y = gene1_mean, group = 1, color = gene1), size = 1) +
  geom_point(aes(y = gene1_mean, color = gene1), size = 3) +
  geom_errorbar(aes(ymin = gene1_CI_low, ymax = gene1_CI_high, color = gene1), width = 0.2) +
  
  # Plot Gene 2 Mean with Error Bars
  geom_line(aes(y = gene2_mean, group = 1, color = gene2), size = 1) +
  geom_point(aes(y = gene2_mean, color = gene2), size = 3) +
  geom_errorbar(aes(ymin = gene2_CI_low, ymax = gene2_CI_high, color = gene2), width = 0.2) +
  
  # Formatting
  labs(title = "Average gene expression across clusters",
       x = "Cell Type",
       y = "Expression Level") +
  scale_color_manual(values = c("blue", "red")) +  # Custom colors for Gene 1 & Gene 2
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
cat("Saved line plot to 'average_gene_expression_across_clusters.pdf'.\n")

cat("Analysis complete.\n")

library(Seurat)
library(optparse)

# Define command-line arguments
option_list = list(
  make_option(c("--input_path"), type="character", default=NULL, help="Path to the Seurat RDS object", metavar="character"),
  make_option(c("--data_type"), type="character", default="counts", help="Type of data to extract from Seurat object", metavar="character"),
  make_option(c("--correlation_method"), type="character", default="pearson", help="Correlation method (pearson, spearman, kendall)", metavar="character"),
  make_option(c("--gene"), type="character", default=NULL, help="Gene ID to compute correlation against", metavar="character"),
  make_option(c("--output_path"), type="character", default=NULL, help="Output file path for saving results", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Function to compute rank correlation
rank_correlation = function(Obj_path, data_type, correlation_method, gene, output_path){
  Obj = readRDS(Obj_path)
  
  valid_slots <- c("counts", "data", "scale.data")
  if (!(data_type %in% valid_slots)) {
    stop(paste("Error: Invalid data type. Choose from:", paste(valid_slots, collapse = ", ")))
  }
  
  data <- GetAssayData(Obj, layer = data_type)
  
  # Check if the gene exists
  if (!(gene %in% rownames(data))) {
    stop(paste("Error: Gene", gene, "not found in dataset"))
  }
  genes_to_check <- setdiff(rownames(data), gene)
  
  correlations <- apply(data[genes_to_check, ], 1, function(x) cor(as.numeric(x), as.numeric(data[gene, ]), method = correlation_method, use = "complete.obs"))
  
  df = data.frame(
    genes = genes_to_check,
    correlation = round(correlations, 5)
  )
  
  sorted_df <- df[order(-df$correlation), ]
  write.csv(sorted_df, output_path, row.names = FALSE)
}

# Run the function with command-line arguments
rank_correlation(
  Obj_path = opt$input_path,
  data_type = opt$data_type,
  correlation_method = opt$correlation_method,
  gene = opt$gene,
  output_path = opt$output_path
)

import numpy as np
import pandas as pd
import argparse
from scipy.stats import pearsonr, spearmanr, kendalltau
import os

parser = argparse.ArgumentParser(description="Compute correlation between TF and target genes in original and imputed datasets.")
parser.add_argument("--original_data", required=True, help="Path to the original data CSV file.")
parser.add_argument("--imputed_data", required=True, help="Path to the imputed data CSV file.")
parser.add_argument("--tf_target_data", required=True, help="Path to the TF-target gene CSV file.")
parser.add_argument("--threshold", type=int, required=True, help="Threshold number of target genes.")
parser.add_argument("--correlation_method", choices=["pearson", "spearman", "kendall"], default="pearson", help="Correlation method to use (default: pearson).")
parser.add_argument("--output_file", required=True, help="Path to the output CSV file.")

args = parser.parse_args()

# Load TF-target data
gene_variances = pd.read_csv(args.tf_target_data)
value_targets = gene_variances["TF"].value_counts()
filtered_TF = value_targets[value_targets < args.threshold].sort_values(ascending=True).index

# Load original and imputed data
original_data = pd.read_csv(args.original_data, index_col=0)
original_data.index = original_data.index.str.replace('gene:', '', regex=True)

imputed_data = pd.read_csv(args.imputed_data, index_col=0)
imputed_data.index = imputed_data.index.str.replace('gene:', '', regex=True)

# Lists to store results
TF_list = []
TF_ori_cor = []
TF_imputed_cor = []
TF_ori_express_list = []
TF_imputed_express_list = []

for i in filtered_TF:
    # Get target genes for this TF
    target_genes = gene_variances[gene_variances['TF'] == i]['Target']
    valid_indices_original = target_genes[target_genes.isin(original_data.index)]
    valid_indices_imputed = target_genes[target_genes.isin(imputed_data.index)]
    
    # Process original data
    if not valid_indices_original.empty:
        filtered_data_original = original_data.loc[valid_indices_original]
        TF_ori_express = filtered_data_original.mean().mean()
        
        original_data_transposed = filtered_data_original.T
        correlation_matrix_original = original_data_transposed.corr(method=args.correlation_method)
        melted_df_original = correlation_matrix_original.reset_index().melt(id_vars='index', var_name='colname', value_name='value')
        melted_df_original = melted_df_original.rename(columns={'index': 'gene1', 'colname': 'gene2', 'value': 'cor'})
        
        melted_df_original.to_csv(f'TF_{i}.csv', index=False)
        TF_ori = melted_df_original['cor'].mean()
    else:
        print(f"No valid targets found for TF in original data: {i}")
        TF_ori = None
    
    # Process imputed data
    if not valid_indices_imputed.empty:
        filtered_data_imputed = imputed_data.loc[valid_indices_imputed]
        TF_imputed_express = filtered_data_imputed.mean().mean()
        
        imputed_data_transposed = filtered_data_imputed.T
        correlation_matrix_imputed = imputed_data_transposed.corr(method=args.correlation_method)
        melted_df_imputed = correlation_matrix_imputed.reset_index().melt(id_vars='index', var_name='colname', value_name='value')
        melted_df_imputed = melted_df_imputed.rename(columns={'index': 'gene1', 'colname': 'gene2', 'value': 'cor'})
        
        melted_df_imputed.to_csv(f'Imputed_TF_{i}.csv', index=False)
        TF_imputed = melted_df_imputed['cor'].mean()
    else:
        print(f"No valid targets found for TF in imputed data: {i}")
        TF_imputed = None
    
    TF_list.append(i)
    TF_ori_cor.append(TF_ori)
    TF_imputed_cor.append(TF_imputed)
    TF_ori_express_list.append(TF_ori_express)
    TF_imputed_express_list.append(TF_imputed_express)

# Create final DataFrame
result_df = pd.DataFrame({
    'TF': TF_list,
    'TF_ori_cor': TF_ori_cor,
    'TF_imputed_cor': TF_imputed_cor,
    'TF_ori_express': TF_ori_express_list,
    'TF_imputed_express': TF_imputed_express_list
})

# Save the result to a CSV file
result_df.to_csv(args.output_file, index=False)
print(f"Results saved to {args.output_file}")

library(ActivePathways)
library(ggplot2)

path="data_ER-HER2-_1.csv"

pvals_FCs <-read.csv(path, stringsAsFactors = FALSE)

pval_matrix <- data.frame(
  row.names = pvals_FCs$gene,
  me = pvals_FCs$me_pval)
#  rna = pvals_FCs$rna_pval, 
#  protein = pvals_FCs$protein_pval,
#  p_protein = pvals_FCs$p_protein_pval)

pval_matrix <- as.matrix(pval_matrix)

pval_matrix[is.na(pval_matrix)] <- 1

dir_matrix <- data.frame(
  row.names = pvals_FCs$gene,
  me = pvals_FCs$me_log2fc)
#  rna = pvals_FCs$rna_log2fc,
#  protein = pvals_FCs$protein_log2fc,
#  p_protein = pvals_FCs$p_protein_log2fc)

dir_matrix <- as.matrix(dir_matrix)

dir_matrix <- sign(dir_matrix)
dir_matrix[is.na(dir_matrix)] <- 0

example_genes <- c('ACTN4','PIK3R4','PPIL1','NELFE','LUZP1','ITGB2')
dir_matrix[example_genes,]

constraints_vector <- c(1)

directional_merged_pvals <- merge_p_values(pval_matrix, 
                                           method = "DPM", dir_matrix, constraints_vector)
directional_merged_pvals_mtx <- as.matrix(directional_merged_pvals)


path2="./data/h.all.v2024.1.Hs.symbols.gmt"

enriched_pathways_directional <- ActivePathways(
  pval_matrix, gmt = path2, cytoscape_file_tag = "./Directional_", merge_method = "DPM",
  scores_direction = dir_matrix, constraints_vector = constraints_vector,significant = 1)

result_file <- paste('./ActivePathways_results_me_1.csv', sep = '/')
export_as_CSV(enriched_pathways_directional, result_file)


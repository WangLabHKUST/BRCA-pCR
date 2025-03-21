library(ActivePathways)
library(ggplot2)

path="data_ER-HER2-_1.csv"

pvals_FCs <-read.csv(path, stringsAsFactors = FALSE)

pval_matrix <- data.frame(
  row.names = pvals_FCs$gene,
  me = pvals_FCs$me_pval,
  rna = pvals_FCs$rna_pval, 
  protein = pvals_FCs$protein_pval,
  p_protein = pvals_FCs$p_protein_pval)

pval_matrix <- as.matrix(pval_matrix)

pval_matrix[is.na(pval_matrix)] <- 1

dir_matrix <- data.frame(
  row.names = pvals_FCs$gene,
  me = pvals_FCs$me_log2fc,
  rna = pvals_FCs$rna_log2fc,
  protein = pvals_FCs$protein_log2fc,
  p_protein = pvals_FCs$p_protein_log2fc)

dir_matrix <- as.matrix(dir_matrix)

dir_matrix <- sign(dir_matrix)
dir_matrix[is.na(dir_matrix)] <- 0
dir_matrix[, "p_protein"] <- 0
example_genes <- c('ACTN4','PIK3R4','PPIL1','NELFE','LUZP1','ITGB2')
dir_matrix[example_genes,]

constraints_vector <- c(1,-1,-1,0)

directional_merged_pvals <- merge_p_values(pval_matrix, 
                                           method = "DPM", dir_matrix, constraints_vector)
merged_pvals <- merge_p_values(pval_matrix, method = "Brown")
directional_merged_pvals_mtx <- as.matrix(directional_merged_pvals)

lineplot_df <- data.frame(original = -log10(merged_pvals),
                          modified = -log10(directional_merged_pvals))

ggplot(lineplot_df) +
  geom_point(size = 2.4, shape = 19,
             aes(original, modified,
                 color = ifelse(original <= -log10(0.05),"gray",
                                ifelse(modified > -log10(0.05),"#1F449C","#F05039")))) +
  labs(title = "",
       x ="Merged -log10(P)",
       y = "Directional Merged -log10(P)") + 
  geom_hline(yintercept = 1.301, linetype = "dashed",
             col = 'black', size = 0.5) +
  geom_vline(xintercept = 1.301, linetype = "dashed",
             col = "black", size = 0.5) + 
  geom_abline(size = 0.5, slope = 1,intercept = 0) +
  scale_color_identity()

path2="./data/h.all.v2024.1.Hs.symbols.gmt"

enriched_pathways_directional <- ActivePathways(
  pval_matrix, gmt = path2, cytoscape_file_tag = "./Directional_", merge_method = "DPM",
  scores_direction = dir_matrix, constraints_vector = constraints_vector,significant = 1)

result_file <- paste('./ActivePathways_results_all_1.csv', sep = '/')
export_as_CSV (enriched_pathways_directional, result_file)
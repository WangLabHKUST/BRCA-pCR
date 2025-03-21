setwd("/Users/zmoad/OneDrive - HKUST Connect/HKUST/BreastGGH/12.Submission/03.ScienceAdvance/Code/subcode")

# cd /scratch/PI/jgwang/zmoad/Projects/BreastGGH/01.Data/07.WGBS/04.wgstoolsClean
library("data.table")
library("matrixStats")
library("ggplot2")

block_matrix = as.data.frame(fread("zbeta_to_table.txt"))
dim(block_matrix)
rownames(block_matrix) = paste(block_matrix$chr,block_matrix$start,block_matrix$end,sep = "_")
block_matrix = block_matrix[,-c(1:5)]
Callrate = rowSums(!is.na(block_matrix))/ncol(block_matrix)
block_matrix.pass_callrate = block_matrix[Callrate == 1,]
dim(block_matrix.pass_callrate)
passCallRate.o = cbind(rownames(block_matrix.pass_callrate),block_matrix.pass_callrate)
colnames(passCallRate.o)[1] = "block"
write.table(passCallRate.o,file = "table.2-0.block.by.sample.100percent.callrate.txt",col.names = T,row.names = F,sep = "\t",quote = F)



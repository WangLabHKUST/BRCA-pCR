setwd('/Users/zmoad/OneDrive - HKUST Connect/HKUST/BreastGGH/12.Submission/03.ScienceAdvance/Code/subcode')

# this code is used to calculate the correlation between segment methylation level and expression level.

tpm = read.csv("../data/table.1-2.TPM.gene.by.Tumor.csv",header = T,row.names = 1)
tpm[1:10,1:10]
tpm_genes = colnames(tpm)
tpm$SampleID = rownames(tpm)


RNAanno = fread('../data/table.3-2.wgbstools_block.anno.txt')
dim(RNAanno)


wgbs_res = as.data.frame(fread(paste("../data/Methylation_segment_beta_value.tsv",sep = "")))
rownames(wgbs_res) = wgbs_res$block
wgbs_res = wgbs_res[,-c(1)]
wgbs_res = as.data.frame(t(wgbs_res))
wgbs_res$SampleID = rownames(wgbs_res)
wgbs_res[1:10,1:10]
dim(wgbs_res)
  
# /Users/zmoad/OneDrive - HKUST Connect/HKUST/BreastGGH/01.Data/table.4-2.clinical.with.multi-omics-allmarker.txt
clinical = fread("../data/clinical.txt")
clinical2 = clinical[!clinical$SampleID %in% c("PS066_T","PS106_T"),]


dmr = RNAanno[RNAanno$key %in% colnames(wgbs_res),]
dmr = dmr[dmr$SYMBOL %in% tpm_genes,]
dim(dmr)# 964627     18
dmr$key


wgbs_res.sub = wgbs_res[,c("SampleID",dmr$key[1:5])]

# Start timer
start_time <- Sys.time()

# Preallocate a data frame for results
dmr_row = nrow(dmr)
results <- data.frame(SYMBOL = character(dmr_row),
                      key = character(dmr_row),
                      p_value = numeric(dmr_row),
                      estimate = numeric(dmr_row),
                      stringsAsFactors = FALSE)
dim(results)
# Loop through rows
for (i in 1:dmr_row) {  # nrow(dmr)
  key.tmp <- dmr$key[i]
  SYMBOL.tmp <- dmr$SYMBOL[i]
  
  # Use data.table for faster subsetting (if applicable)
  
  # Merge data frames
  rna_meth <- merge(tpm[, c("SampleID", SYMBOL.tmp)], wgbs_res[, c("SampleID", key.tmp)], by = "SampleID")
  colnames(rna_meth) <- c("SampleID", "rna", "meth")
  
  # Calculate Spearman correlation
  corTest_unlog <- cor.test(rna_meth$rna, rna_meth$meth, method = "spearman")
  
  # Store results
  results[i, ] <- c(SYMBOL.tmp, key.tmp, corTest_unlog$p.value, corTest_unlog$estimate)
  
  cat(i, "\n")
}

end_time <- Sys.time()
elapsed_time <- end_time - start_time

# Print elapsed time
cat("Elapsed time:", elapsed_time, "seconds\n")

overctor.df = as.data.frame(results)
colnames(overctor.df) = c("GeneName","key","RNA_p","RNA_r")
write.table(overctor.df,file = paste("table.2-1.association.meth_and_RNA.txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)






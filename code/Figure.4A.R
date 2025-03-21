library("data.table")
library("ggplot2")
library("ggsci")
library("Rtsne")
library("cowplot")
library("ggpubr")

clinical = fread("./data/clinical.txt")

# Generate some example data (replace this with your actual data)
set.seed(123)

gene_data = as.data.frame(fread("./data/table.1-2.TPM.gene.by.Tumor.csv"))
rownames(gene_data) = gene_data$SampleID
gene_data = gene_data[,-1]
gene_data[1:10,1:10]
dim(gene_data)

if(T){
  # Step 1: Filter columns where the proportion of zeros is less than 10%
  gene_data <- gene_data[, colMeans(gene_data == 0) < 0.1]
  dim(gene_data)
  
  # Step 2: Apply log2 transformation after adding 0.01
  gene_data <- log2(gene_data + 0.01) #行为样本列为基因
  dim(gene_data)
}


# Run t-SNE

set.seed(123)
perplexity = 5
tsne_result <- Rtsne(gene_data, dims = 2, perplexity = perplexity, verbose = TRUE, check_duplicates = FALSE)

tsne_data <- data.frame(tsne_result$Y)
colnames(tsne_data) <- c("tSNE1", "tSNE2")
tsne_data$SampleID = rownames(gene_data)

pca_data2 = merge(tsne_data,clinical,by = "SampleID")
pca_data2$pCR = pca_data2$pCR_consensus
pca_data2$pCR = factor(pca_data2$pCR,levels = c("pCR","RD"))
pca_data2$subtype = factor(pca_data2$subtype, levels = c("ER+HER2+","ER-HER2+","ER-HER2-"))

bySubtype = ggplot(pca_data2, aes(x = tSNE1, y = tSNE2, col = subtype , shape = pCR)) + 
  geom_point(size = 2.8,stroke = 0.9) + #shape = 21,
  xlab(paste0("tSNE1")) +
  ylab(paste0("tSNE2")) +
  ggtitle("Subtype") + theme_bw(base_size = 20) + scale_color_manual(values=c("#CC9BC6", "#70B087", "#777777")) + scale_fill_manual(values=c("#CC9BC6", "#70B087", "#777777")) + scale_shape_manual(values=c(20, 21))  #+ scale_color_npg()
bySubtype


pdf(paste("figure.4-1.tSNE.perplexity.",perplexity,".pdf",sep = ""),width = 7,height = 5)
plot(bySubtype)
dev.off()






library("data.table")
library("dplyr")
library("qusage")

bc_exp_raw = as.data.frame(read.csv("./data/protein.csv"))
rownames(bc_exp_raw) = paste(bc_exp_raw[,1],"_T",sep = "")
bc_exp_raw = as.data.frame(t(bc_exp_raw[,-1]))
dim(bc_exp_raw)

if(TRUE){
  for (i in 1:ncol(bc_exp_raw)){
    bc_exp_raw[,i] = 2^bc_exp_raw[,i]
  }
}

clinical = as.data.frame(fread("./data/clinical.txt"))
clinical = clinical[clinical$SampleID %in% colnames(bc_exp_raw),]
clinical$pCR = clinical$pCR_consensus

bc_exp_raw = bc_exp_raw[,clinical$SampleID]
RNA_samples = colnames(bc_exp_raw)
clinical = clinical[clinical$SampleID %in% RNA_samples,]

subtypes = list("ER-HER2+","ER+HER2+","ER-HER2-",c("ER-HER2+","ER+HER2+","ER-HER2-"))
namesall = list("ER-HER2+","ER+HER2+","ER-HER2-","all")

for (K in 1:length(subtypes)){
  subtype.tmp = subtypes[[K]]
  subtype.name = namesall[[K]]
  
  clinical.sub = clinical[clinical$subtype %in% subtype.tmp,]
  table(clinical.sub$pCR)
  conditions = as.factor(clinical.sub$pCR)
  
  count_norm = bc_exp_raw[,clinical.sub$SampleID]
  
  rowmean_V = rowMeans(count_norm)
  
  # Run the Wilcoxon rank-sum test for each gene
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  # Calculate the fold-change for each gene
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  
  fold_calculator = "byMedian"
  if(fold_calculator == "byMean"){
    foldChanges=rowMeans(dataCon1,na.rm = T) - rowMeans(dataCon2,na.rm = T)
    pCRmean = rowMeans(dataCon1,na.rm = T)
    RDmean = rowMeans(dataCon2,na.rm = T)
    outRst<-data.frame(GeneName =rownames(count_norm),log2foldChange=foldChanges, pValues=pvalues, FDR=fdr,pCRmean = pCRmean, RDmean = RDmean)
  }else{
    foldChanges=log2(rowMedians(as.matrix(dataCon1),na.rm = T)/rowMedians(as.matrix(dataCon2),na.rm = T))
    pCRmedian = rowMedians(as.matrix(dataCon1),na.rm = T)
    RDmedian = rowMedians(as.matrix(dataCon2),na.rm = T)
    outRst<-data.frame(GeneName =rownames(count_norm),log2foldChange=foldChanges, pValues=pvalues, FDR=fdr,pCRmedian = pCRmedian, RDmedian = RDmedian)
  }
  
  dim(outRst)
  outRst = outRst[order(outRst$pValues,decreasing = F),]
  write.table(outRst, file=paste("table.5A.DEPs.",subtype.name,".txt",sep = ""),sep="\t", quote=F,row.names = F,col.names = T)
}



library("dplyr")
library("ggplot2")
library("readr")
library("ggrepel")
library("Matrix.utils")

res_gs <- list()

updated_Subtype.list = c("ER-HER2+","ER+HER2+","ER-HER2-","all")#,"HER2+","All"
# Volcano plots and GSEA analysis
for (i in 1:length(updated_Subtype.list)) {
  subtype.tmp = updated_Subtype.list[i]
  sig_c1 <- as.data.frame(fread(paste("table.5A.DEPs.",subtype.tmp,".txt",sep = "")))
  colnames(sig_c1)[1] = "GENE"
  rownames(sig_c1) = sig_c1$GENE
  
  # Highlight breast cancer drivers
  BRCA_related = read.csv('./data/BRCA_related.genes.csv')
  BRCA_related <- unique(c(BRCA_related$core_genes,BRCA_related$extended_genes))
  
  sig_c1$color = ifelse(sig_c1$pValues > 0.05,"#808080",ifelse(sig_c1$log2foldChange > 0,"#FB8173","#80B1D3"))
  table(sig_c1$color)
  
  index_highlight = sig_c1$GENE %in% BRCA_related & 
    sig_c1$pValues < 0.05 & 
    sig_c1$log2foldChange !=0
  
  sig_c1$size = 1
  size_updated = -log10(sig_c1$pValues)
  sig_c1$size[index_highlight] = size_updated[index_highlight]
  sig_c1[is.na(index_highlight)]
  
  sig_c1.back = sig_c1[sig_c1$size==1,]
  sig_c1.highlight = sig_c1[sig_c1$size!=1,]
  
  pvalue.top10 = sig_c1.back[order(sig_c1.back$pValues),]
  pvalue.top10 = pvalue.top10[1:10,]
  sig_c1.back$bg_label = NA
  index_addTop_label = sig_c1.back$GENE %in% pvalue.top10$GENE
  sig_c1.back$bg_label[index_addTop_label] = sig_c1.back$GENE[index_addTop_label]
  
  sig_c1.highlight$label = sig_c1.highlight$GENE
  sig_c1.highlight = sig_c1.highlight[order(sig_c1.highlight$size,decreasing = T),]
  top_number = 20
  if(nrow(sig_c1.highlight)> top_number){
    sig_c1.highlight$label[top_number:nrow(sig_c1.highlight)] = NA
  }
  
  
  addupdat = rBind.fill(sig_c1.back,sig_c1.highlight)
  x_lim = max(abs(addupdat$log2foldChange))
  head(sig_c1.back)
  volcano_DEG = ggplot(sig_c1.back, aes(x = log2foldChange, y = -log10(pValues))) +
    geom_point(size = sig_c1.back$size, col = sig_c1.back$color) +
    geom_point(data = sig_c1.highlight, aes(size = sig_c1.highlight$size), shape = 21,fill = sig_c1.highlight$color,stroke = 0.7) + 
    theme_classic(base_size = 24) +
    xlim(-2.5, 2.5) +
    ylim(0, 4) +
    labs(title = paste(subtype.tmp,'Pro'),
         x = 'Log2 FoldChange', #
         y = '-Log10 (P-value)') +
    theme(legend.position = "none")+ 
    theme(axis.text.x=element_text(size=20,color = "black"),
          axis.text.y=element_text(size=20,color = "black"))+ 
    geom_text_repel(data = sig_c1.highlight,label = sig_c1.highlight$label,size = 7.5,max.overlaps = 9,box.padding = 0.4)
  
  
  pdf(file = paste0('Figure.5A.DEPs.volcano_',subtype.tmp, '.pdf'),width = 6,height= 6)
  plot(volcano_DEG)
  dev.off()
  
}







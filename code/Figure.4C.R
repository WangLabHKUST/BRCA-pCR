library("DESeq2")
library("dplyr")

clinical = as.data.frame(fread("./data/clinical.txt"))
clinical = clinical[,c("SampleID","pCR_consensus","subtype")]
colnames(clinical) = c("SampleID","pCR","Molecular_Subtype")
clinical$pCR = ifelse(clinical$pCR == "pCR",1,0)


countdat.raw = as.data.frame(fread("./data/RSEM_RNA.count.txt"))
countdat.raw = countdat.raw[!duplicated(countdat.raw$GeneSymbol),]
countdat.raw = countdat.raw[,-c(1,3)]
countdat = countdat.raw[,c(intersect(clinical$SampleID,colnames(countdat.raw)))]
rownames(countdat) = countdat.raw$GeneSymbol
countdat = as.data.frame(t(countdat))
for (i in 1:ncol(countdat)){
  countdat[,i] = as.integer(countdat[,i])
}
countdat$SampleID = rownames(countdat)
countdat = as.data.frame(merge(clinical,countdat,by = "SampleID"))
rownames(countdat) = countdat$SampleID
countdat = countdat[,-1]
countdat[1:10,1:10]


tpm.raw = as.data.frame(fread("./data/table.1-2.TPM.gene.by.Tumor.csv"))
tpm = merge(clinical,tpm.raw,by = "SampleID")
rownames(tpm) = tpm$SampleID
tpm = tpm[,-1]
tpm[1:10,1:10]
dim(tpm)


# Initialize an empty list to store results
deres <- list()

Molecular_Subtype = list("ER-HER2+","ER-HER2-","ER+HER2+",c("ER-HER2+","ER-HER2-","ER+HER2+"))
names_verctors = c("ER-HER2+","ER-HER2-","ER+HER2+","all")
# Loop through unique molecular subtypes
for (i in 1:length(Molecular_Subtype)) {
  Molecular_Subtype.tmp = Molecular_Subtype[[i]]
  names.tmp = names_verctors[[i]]
  # Subset and sort the DataFrame
  pd_df_t <- countdat %>% filter(Molecular_Subtype %in% Molecular_Subtype.tmp) %>% arrange(pCR)
  
  tmp_t_ctr <- tpm %>% filter(Molecular_Subtype %in% Molecular_Subtype.tmp & pCR == 1)
  tmp_t_case <- tpm %>% filter(Molecular_Subtype %in% Molecular_Subtype.tmp & pCR == 0)
  
  # Create the count matrix
  count_matrix <- as.data.frame(t(pd_df_t[ , !(names(pd_df_t) %in% c("pCR", "Molecular_Subtype"))]))
  
  # Create the design matrix
  design_matrix <- data.frame(pCR = as.factor(pd_df_t$pCR))
  
  # Print the design matrix
  print(design_matrix)
  
  # Create the DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                colData = design_matrix,
                                design = ~ pCR)
  dds <- DESeq(dds)
  
  # Get results and convert to a data frame
  pd_from_r_df <- as.data.frame(results(dds))
  
  
  # Perform rank-sum test
  ranksum_pdict <- list()
  geneNames = setdiff(colnames(tmp_t_ctr),c("pCR","Molecular_Subtype"))
  for (j in geneNames) {
    ranksum_result <- wilcox.test(tmp_t_ctr[[j]], tmp_t_case[[j]])#, exact = TRUE,alternative = "two.sided",correct = TRUE
    ranksum_pdict[[j]] <- ranksum_result$p.value
  }
  
  # Store results in the list
  pd_from_r_df <- pd_from_r_df[order(pd_from_r_df$padj), ]  # Sort by padj
  pd_from_r_df$prank <- as.numeric(as.character(ranksum_pdict[rownames(pd_from_r_df)]))
  deres[[i]] <- pd_from_r_df
  
  
  # Write to CSV
  write.csv(pd_from_r_df, file = paste0('table.4C-1.DEGs-newTPM_', names.tmp, '.csv'), row.names = TRUE)
  
}



#######  plot volcano
library("dplyr")
library("ggplot2")
library("readr")
library("ggrepel")
library("tidyr")
library("data.table")

# Initialize an empty list to store results
res_gs <- list()

subtypes = c("ER-HER2+", "ER+HER2+", "ER-HER2-","all")

# Volcano plots and GSEA analysis
for (i in subtypes) {
  clinical = as.data.frame(fread("./data/clinical.txt"))
  clinical = clinical[clinical$RNA!="",]
  if(i!="all"){
    clinical = clinical[clinical$subtype == i,]
  }
  pCR_consensus = table(clinical$pCR_consensus)
  title.tmp = paste(paste(names(pCR_consensus),pCR_consensus,sep = "="),collapse = ",")
  
  deres = as.data.frame(fread(paste0('table.4C-1.DEGs-newTPM_', i, '.csv')))
  sig_c1 <- deres %>%
    drop_na() #%>%
  sig_c1$GENE = sig_c1$V1
  
  
  # Highlight breast cancer drivers
  BRCA_related = read.csv('./data/BRCA_related.genes.csv')
  BRCA_related <- unique(c(BRCA_related$core_genes,BRCA_related$extended_genes))
  length(BRCA_related)
  
  # sig_c1$pvalue要求DEseq2的p value显著
  fc_cutoff = 0
  sig_c1$color = ifelse(sig_c1$pvalue > 0.05,"#808080",ifelse(sig_c1$log2FoldChange > fc_cutoff,"#FB8173",ifelse(sig_c1$log2FoldChange < -fc_cutoff,"#80B1D3", "#808080")))
  table(sig_c1$color)
  
  # highlight的必须在网络定义的基因list里面，并且rank sum test的pvalue < 0.05
  index_highlight = sig_c1$GENE %in% BRCA_related & sig_c1$prank < 0.1 & sig_c1$pvalue < 0.05 & abs(sig_c1$log2FoldChange) > fc_cutoff
  sig_c1$size = 1
  # 如果属于这些需要highlight的基因里面，那么就更新圈圈的大小
  size_updated = -log10(sig_c1$pvalue) # abs(sig_c1$log2FoldChange) *  - 0.2
  sig_c1$size[index_highlight] = size_updated[index_highlight]
  
  
  sig_c1.back = sig_c1[sig_c1$size==1,]
  if(fc_cutoff > 0){
    sig_c1.highlight = sig_c1[sig_c1$size!=1 & sig_c1$color != "#808080",]# 
  }else{
    sig_c1.highlight = sig_c1[sig_c1$size!=1,]# & sig_c1$color != "#808080"
  }
  
  pvalue.top10 = sig_c1.back[order(sig_c1.back$pvalue),]
  pvalue.top10 = pvalue.top10[1:3,]
  sig_c1.back$bg_label = NA
  index_addTop_label = sig_c1.back$GENE %in% pvalue.top10$GENE
  sig_c1.back$bg_label[index_addTop_label] = sig_c1.back$GENE[index_addTop_label]
  
  
  xlim.v = max(abs(sig_c1.back$log2FoldChange))+0.1
  if(i == "all"){
    max.overlaps = 13
  }else{
    max.overlaps = 13
  }
  
  if(FALSE){
    sig_c1.highlight = sig_c1.highlight[order(sig_c1.highlight$size,decreasing = T),]
    sig_c1.highlight$label = sig_c1.highlight$GENE
    if(nrow(sig_c1.highlight) > 6){
      sig_c1.highlight$label[7:nrow(sig_c1.highlight)] = NA
    }
  }
  
  sig_c1.highlight = sig_c1.highlight[order(sig_c1.highlight$size,decreasing = T),]
  if(i == "ER-HER2-"){
    sig_c1.highlight$label = sig_c1.highlight$GENE
    sig_c1.highlight$label[!sig_c1.highlight$label %in% c("NR3C2","AR","MKI67","E2F3","ING1","CDK1")]= NA
    sig_c1.add = sig_c1[sig_c1$GENE %in% c("FDCSP","CCNE1","SERHL2"),]
    sig_c1.add$label = sig_c1.add$GENE
    sig_c1.highlight = rbind(sig_c1.highlight,sig_c1.add)
    sig_c1[sig_c1$GENE %in% c("CCNE1",'ING1',"SERHL2"),]
    
    
  }else if(i == "ER-HER2+"){
    sig_c1.highlight$label = sig_c1.highlight$GENE
    sig_c1.highlight$label[!sig_c1.highlight$label %in% c("GFRA3","FOXO3","DGKB","CBLC","IGF1","ERBB2")]= NA
    sig_c1[sig_c1$GENE %in% c("ERBB2",'IGF1',"CBLC"),]
    sig_c1.add = sig_c1[sig_c1$GENE %in% c("ERBB2"),]
    sig_c1.add$label = sig_c1.add$GENE
    sig_c1.highlight = rbind(sig_c1.highlight,sig_c1.add)
    
  }else if(i == 'ER+HER2+'){
    print("label all genes")
    sig_c1.highlight$label = sig_c1.highlight$GENE
    sig_c1.highlight$label[!sig_c1.highlight$label %in% c("STAT1","ABCC3","CHGA")]= NA
    sig_c1.add = sig_c1[sig_c1$GENE %in% c("LRP1B","ATG13","CASP1","FGFR1"),]
    sig_c1.add$label = sig_c1.add$GENE
    sig_c1.highlight = rbind(sig_c1.highlight,sig_c1.add)
    sig_c1[sig_c1$GENE %in% c("FGFR1",'ATG13',"CASP1","STAT1","LRP1B"),]
    
  }else if(i == 'all'){
    sig_c1.highlight = sig_c1.highlight[order(sig_c1.highlight$size,decreasing = T),]
    sig_c1.highlight$label = sig_c1.highlight$GENE
    sig_c1.highlight$label[!sig_c1.highlight$label %in% c("ERBB2","S100A8","S100A9","CASP3","NR3C2","TP63","DGKB")]= NA
    # if(nrow(sig_c1.highlight) > 6){
    #   sig_c1.highlight$label[7:nrow(sig_c1.highlight)] = NA
    # }
  }
  
  
  # Create the volcano plot
  volcano_DEG = ggplot(sig_c1.back, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(size = sig_c1.back$size, col = sig_c1.back$color) +
    geom_point(data = sig_c1.highlight, aes(size = sig_c1.highlight$size), shape = 21,fill = sig_c1.highlight$color,stroke = 0.6) +  # Highlighted points
    #scale_color_manual(values = c("grey", "blue", "red")) +
    theme_classic(base_size = 18) +
    xlim(-xlim.v, xlim.v) +
    ylim(-0.1, 11) +
    labs(title = paste(i,title.tmp),
         x = 'Log2 Foldchange',
         y = '-Log10 (P-value)') +
    theme(legend.position = "none")+ 
    geom_text_repel(data = sig_c1.highlight,label = sig_c1.highlight$label,size = 6,max.overlaps = max.overlaps,box.padding = 0.7,point.padding = 0.5,fontface = "plain")+
    # geom_text_repel(data = sig_c1.back,label = sig_c1.back$bg_label,max.overlaps = 20,size = 6,col = "black")+ 
    theme(axis.text.x=element_text(size=17,color = "black"), axis.text.y=element_text(size=17,color = "black"), axis.title.y=element_text(size=17,color = "black"), axis.title.x=element_text(size=17,color = "black")) + scale_colour_manual(values = c("#20BEE6","#FF7878")) 
  
  
  # Uncomment to save the plot
  pdf(file = paste0('figure.4C.DEG.volcano_', i, '.pdf'),width = 5,height= 5)#.oldLabel
  plot(volcano_DEG)
  dev.off()
  
  
}










library("DESeq2")
library("dplyr")
library("ggplot2")
library("readr")
library("ggrepel")
library("data.table")
library("stringr")
library("Matrix.utils")


DMRanno = as.data.frame(fread(paste("./data/table.3-2.wgbstools_block.anno.txt",sep = "")))
DMRanno$SYMBOL[DMRanno$key == "chr9_21991219_21993173"] = "CDKN2A"
DMRanno[DMRanno$key == "chr9_21991219_21993173",]
DMRanno = DMRanno[!is.na(DMRanno$SYMBOL),]

callrate100 = fread("./data/Methylation_segment_beta_value.txt")
callrate100 = callrate100$block


tpm = read.csv("./data/table.1-2.TPM.gene.by.Tumor.csv",header = T,row.names = 1)
tpm_genes = colnames(tpm)


block_var = fread("./data/table.2-2.100percent.callrate.Variance.txt")
colnames(block_var)[1] = "key"
dim(block_var)

new_names = c("ER-HER2-","ER-HER2+","ER+HER2+","all")


# Volcano plots and GSEA analysis
for (i in 1:length(new_names)) {
  new_names.tmp = new_names[i]
  
  sig_c1 = fread(paste("table.3C.DMRs.",new_names.tmp,".txt",sep = ""))
  colnames(sig_c1)[1] = "key"
  
  
  sig_c1 = merge(sig_c1,block_var,by = "key")
  sig_c1 = merge(sig_c1,DMRanno,by = "key")
  
  # add segment-RNA association information
  assocRNA = fread(paste("./data/table.2-1.association.meth_and_RNA.txt",sep = ""))
  sig_c1 = merge(sig_c1,assocRNA,by = "key",all.x = T)
  
  # add BRCANET information
  brcaGene = fread("./data/BRCA_related.genes.csv",header = T)
  brcaGene = unique(c(brcaGene$core_genes,brcaGene$extended_genes))
  sig_c1$InNetwork = ifelse(sig_c1$SYMBOL %in% brcaGene,"InNetwork","NotInNetwork")
  sig_c1$InNetwork[is.na(sig_c1$InNetwork)] = "NotDMR"
  
  sig_c1$pvalue = sig_c1$pValues
  sig_c1$delta = sig_c1$log2foldChange
  
  # 
  delat_cutoff = 0
  sig_c1$color = ifelse(sig_c1$pvalue > 0.05,"#808080",ifelse(sig_c1$delta > delat_cutoff,"#FCA65D",ifelse(sig_c1$delta < -delat_cutoff,"#A2D9A4","#808080")))
  
  
  if(new_names.tmp == "ER-HER2-"){
    sig_c1$RNA_r[sig_c1$key == "chr9_21991219_21993173"] = -0.5463234 
    sig_c1$RNA_p[sig_c1$key == "chr9_21991219_21993173"] = 7.111e-05
  }
  
  sig_c1$RNA_p[is.na(sig_c1$RNA_p)] = 0.99
  sig_c1$RNA_r[is.na(sig_c1$RNA_r)] = 0.01
  
  ## 
  delat_cutoff_show = 0.02
  if(TRUE){
    index_highlight = sig_c1$InNetwork == "InNetwork" & sig_c1$pvalue < 0.05 & abs(sig_c1$delta) > delat_cutoff_show  & sig_c1$RNA_p < 0.05 & sig_c1$RNA_r < 0
    
    size_for_all = -log10(sig_c1$pvalue) * -log10(sig_c1$RNA_p) * 4
    sig_c1[is.na(index_highlight),]
    
    sig_c1$size = 1
    sig_c1$size[index_highlight] = size_for_all[index_highlight]
    sig_c1$isHighlight = ifelse(index_highlight,"highlight","background")
  }
  
  sig_c1$RNA_direction = ifelse(sig_c1$RNA_r > 0,"a.pos_rna","b.neg_rna")
  sig_c1$RNA_pthreshold = ifelse(sig_c1$RNA_p < 0.001,"b.sig_rna","a.non_rna")
  sig_c1$pCR_sigDiff = ifelse(sig_c1$pValues < 0.05,"b.sig_pCR","a.nonsig_pCR")
  sig_c1$area = c(-log10(sig_c1$pValues)) * abs(sig_c1$log2foldChange)*10
  sig_c1 = sig_c1[order(sig_c1$SYMBOL,
                        sig_c1$RNA_direction,
                        sig_c1$RNA_pthreshold,
                        sig_c1$pCR_sigDiff,
                        sig_c1$area,decreasing = T),]
  
  
  codingGene = read.table("./data/Coding.gene.list.txt")$V1
  sig_c1 = sig_c1[sig_c1$SYMBOL %in% codingGene,]
  dim(sig_c1)
  write.table(sig_c1,file = paste("table.3C.DMRs.",new_names.tmp,".txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)
  
  ## after selecting representitive for each gene, now keep only top 50% gene with high variance
  OneBedForOneGene = TRUE
  if(OneBedForOneGene){
    
    sig_c1.rmdup = sig_c1[!duplicated(sig_c1$SYMBOL),]
    sig_c1.rmdup[sig_c1.rmdup$SYMBOL %in% c("CDKN2A"),]
  }
  
  # do variance filter
  ApplyVarFilter = TRUE
  VarPercent = 50
  
  cat("start:",g,"\n")
  delegate = sig_c1.rmdup
  
  TopVar = g/100
  delegate = delegate[order(delegate$variance,decreasing = T),]
  TopCutoff = delegate$variance[nrow(delegate)*TopVar]
  CDKN2A.tmp = delegate[delegate$SYMBOL == "CDKN2A",]
  delegate = delegate[delegate$variance > TopCutoff,]
  
  Total_Feature = nrow(delegate)
  NoteF = paste("# features: ",Total_Feature,"; Variance > ",signif(TopCutoff,3)," (Top",TopVar*100,"%)",sep = "")
  
  
  delegate.highlight = delegate[delegate$InNetwork == "InNetwork" & delegate$pValues < 0.05,]
  
  delegate.highlight.ForLabel = delegate[delegate$InNetwork == "InNetwork" & 
                                         delegate$pValues < 0.05 & 
                                         delegate$RNA_direction == "b.neg_rna" &
                                         delegate$RNA_p < 0.05,]
  delegate.highlight.ForLabel$size = -log10(delegate.highlight.ForLabel$RNA_p)
  
  # Create the volcano plot
  volcano_DEG.raw = ggplot(delegate, aes(x = delta, y = -log10(pvalue))) + #p_2D_KS
    geom_point(col = delegate$color,size = 1) +#size = sig_c1.back$size, 
    geom_point(data = delegate.highlight.ForLabel, aes(size = delegate.highlight.ForLabel$size), shape = 21,fill = delegate.highlight.ForLabel$color,stroke = 1) +  # Highlighted points
    theme_classic(base_size = 40) + ggtitle(paste(new_names.tmp),subtitle = NoteF) + 
    ylim(0, 4.5) +
    labs(
      x = 'Median Difference',
      y = '-Log10 (P-value)') +
    theme(legend.position = "right")+ 
    theme(axis.text.x=element_text(size=35,color = "black"),
          axis.text.y=element_text(size=35,color = "black"))+ 
    scale_size_continuous(limits = c(-log10(0.1), 9),range = c(0,10),breaks = c(1.3,2,4,6,8,9))
  
  volcano_DEG.plain = volcano_DEG.raw
  
  AddBackgrounLabel = TRUE
  if(AddBackgrounLabel){
    volcano_DEG.withLabel = volcano_DEG.raw + geom_text_repel(data = delegate.highlight.ForLabel,label = delegate.highlight.ForLabel$SYMBOL,size = 13,max.overlaps = 25,box.padding = 0.8,fontface = "bold")
  }
  
  # save two plots, one with label and the otehr without labels
  pdf(file = paste0("Figure.3C.DMGs.",new_names.tmp,".OnlyRepresent.Label.pdf"),width = 18*4,height= 10.4*4)
  plot(volcano_DEG.withLabel)
  dev.off()
  
  pdf(file = paste0("Figure.3C.DMGs.",new_names.tmp,".OnlyRepresent.Plain.pdf"),width = 18,height= 10.4)
  plot(volcano_DEG.plain)
  dev.off()
  
  write.table(sig_c1.rmdup,file = paste("table.3C.DMGs.",new_names.tmp,".txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)
  
}






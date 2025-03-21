library("data.table")
library("dplyr")
library("qusage")
library("matrixStats")
library("stringr")
library("ggpubr")
library("ggrepel")

Protein = as.data.frame(read.csv("./data/ssrna-GGH_Protein-Zscore.csv"))
colnames(Protein) = paste("pro",colnames(Protein),sep = "_");colnames(Protein)[1] = "SampleID"
phosphoprotein = as.data.frame(read.csv("./data/ssrna-GGH_phosphoprotein-Zscore.csv"))
colnames(phosphoprotein) = paste("phospro",colnames(phosphoprotein),sep = "_");colnames(phosphoprotein)[1] = "SampleID"

bc_exp_raw = merge(Protein,phosphoprotein,by = "SampleID")

rownames(bc_exp_raw) = paste(bc_exp_raw[,1],"_T",sep = "")
bc_exp_raw = as.data.frame(t(bc_exp_raw[,-1]))
dim(bc_exp_raw)

clinical = as.data.frame(fread("./data/clinical.txt"))

clinical = clinical[clinical$SampleID %in% colnames(bc_exp_raw),]
clinical$pCR = clinical$pCR_consensus


bc_exp_raw = bc_exp_raw[,clinical$SampleID]
RNA_samples = colnames(bc_exp_raw)

clinical = clinical[clinical$SampleID %in% RNA_samples,]

subtypes = list("ER-HER2+","ER+HER2+","ER-HER2-",c("ER-HER2+","ER+HER2+","ER-HER2-"))
namesall = list("ER-HER2+","ER+HER2+","ER-HER2-","all")

outRst.all = c()
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
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]# EAS
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]# EUR
  
  fold_calculator = "byMean"
  if(fold_calculator == "byMean"){
    foldChanges=rowMeans(dataCon1,na.rm = T) - rowMeans(dataCon2,na.rm = T) # EAS -VS- EUR
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
  outRst$subtype = subtype.name
  outRst.all = rbind(outRst.all,outRst)
}


### plot
outRst.all$subtype
namesall = list("ER-HER2+","ER+HER2+","ER-HER2-","all")

for (i in 1:length(namesall)){
  subtype.tmp = namesall[[i]]
  
  de_path = outRst.all[outRst.all$subtype == subtype.tmp,]
  de_path$omics = str_split_fixed(de_path$GeneName,"_",2)[,1]
  
  de_path$pathway = gsub("HALLMARK_","",str_split_fixed(de_path$GeneName,"_",2)[,2])
  
  de_path$label = gsub("_"," ",de_path$pathway)#str_to_title()
  
  SigFeature = de_path[de_path$pValues < 0.05,]
  if(nrow(SigFeature) > 4){
    de_path$label[de_path$pValues > 0.05] = NA
  }else{
    de_path = de_path[order(de_path$pValues),]
    de_path$label[6:nrow(de_path)] = NA
  }
  
  
  
  de_path$color = ifelse(de_path$pValues > 0.05, "#7F7F7F",ifelse(de_path$log2foldChange > 0, "#FF776B","#6EB2D7"))
  de_path$isSig = ifelse(de_path$pValues > 0.05, "non-Sig",ifelse(de_path$log2foldChange > 0, "pCR","RD"))
  
  de_path$size = c(-log10(de_path$pValues)+0.3) * 2 # *10 - 0.4# * 5 # *3.5
  max(de_path$size)
  
  de_path.fake = de_path[1,]
  de_path.fake$log2foldChange = 0
  de_path.fake$pValues = 1
  de_path.fake$size = 3
  
  de_path.in = de_path
  head(de_path.in)
  xlim.tmp = max(abs(de_path.in$log2foldChange))
  # Create the volcano plot
  volcano_DEG = ggplot(de_path.in, aes(x = log2foldChange, y = -log10(pValues), shape = omics, fill = isSig, size = size)) +
    geom_point(stroke = 0.3,col = "white") + 
    theme_classic(base_size = 18) +
    labs(title = subtype.tmp,
         x = 'Log2 Foldchange',
         y = '-Log10 (P-value)') +
    theme(legend.position = "right")+ 
    # ylim(0,ylim.max)+
    xlim(-xlim.tmp,xlim.tmp)+
    geom_text_repel(data = de_path.in,label = stringr::str_wrap(de_path.in$label,width = 20),size = de_path.in$size*0.5,col = "black",max.overlaps = 10)+
    theme(axis.text.x=element_text(size=17,color = "black"), axis.text.y=element_text(size=17,color = "black"), axis.title.y=element_text(size=17,color = "black"), axis.title.x=element_text(size=17,color = "black")) + scale_fill_manual(values = c("gray","#FF7878","#20BEE6")) + scale_color_manual(values = c("gray","#FF7878","#20BEE6"))+ scale_shape_manual(values = c(21,22)) + scale_size_continuous(range = c(3,7)) # 
  volcano_DEG
  
  pdf(file = paste0('Figure.5C.pathway.volcano_', subtype.tmp, '.pdf'),width = 6.2,height= 4.5)
  plot(volcano_DEG)
  dev.off()
  
}




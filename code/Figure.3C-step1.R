library("data.table")
library("dplyr")

subtypes = list("ER-HER2+","ER+HER2+","ER-HER2-",c("ER-HER2+","ER+HER2+","ER-HER2-"))
namesall = list("ER-HER2+","ER+HER2+","ER-HER2-","all")

for (K in 1:length(subtypes)){
  subtype.tmp = subtypes[[K]]
  subtype.name = namesall[[K]]
  # the input is the segment methylaton matrix, which could be downloaded from our websit: http://www.wang-lab-hkust.com:8050/
  bc_exp_raw = as.data.frame(fread(paste("./data/Methylation_segment_beta_value.txt",sep = "")))
  
  rownames(bc_exp_raw) = bc_exp_raw$block
  bc_exp_raw = bc_exp_raw[,-1]
  bc_exp_raw[1:10,1:10]
  
  dim(bc_exp_raw)
  
  clinical = as.data.frame(fread("./data/clinical.txt"))
  clinical = clinical[clinical$Methylation!="",]
  clinical = clinical[!clinical$SampleID %in% c("PS066_T","PS106_T"),]
  clinical$pCR = clinical$pCR_consensus
  clinical = clinical[,c("SampleID","pCR","subtype")]
  clinical = clinical[clinical$SampleID %in% colnames(bc_exp_raw),]
  dim(clinical)
  
  
  clinical.sub = clinical[clinical$subtype %in% subtype.tmp,]
  table(clinical.sub$pCR)
  conditions = as.factor(clinical.sub$pCR)
  
  
  count_norm = bc_exp_raw[,clinical.sub$SampleID]
  
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
  
  
  fold_calculator = "byMean" # byMean
  if(fold_calculator == "byMean"){
    foldChanges=rowMeans(dataCon1,na.rm = T)-rowMeans(dataCon2,na.rm = T)
    con1 = rowMeans(dataCon1)
    con2 = rowMeans(dataCon2)
  }else{
    foldChanges=log2(rowMedians(as.matrix(dataCon1))/rowMedians(as.matrix(dataCon2)))
    con1 = rowMedians(dataCon1)
    con2 = rowMedians(dataCon2)
  }
  
  # Output results based on FDR threshold
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr, Con1 = rowMeans(dataCon1,na.rm = T),Con2 = rowMeans(dataCon2,na.rm = T))
  colnames(outRst)[c(4,5)] = conditionsLevel
  rownames(outRst)=rownames(count_norm)
  
  outRst = outRst[order(outRst$pValues),]
  fdrThres=0.05
  dim(outRst)
  write.table(outRst, file=paste("table.3C.DMRs.",subtype.name,".txt",sep = ""),sep="\t", quote=F,row.names = T,col.names = T)
  
}




# Load necessary libraries
library("ggplot2")
library("dplyr")
library("readr")
library("ggpubr")
library("EnvStats")

# Parameters
subtype <- 'NegNeg'
gene <- 'CDKN2A'
tss <- c(21974857-1000,21974857+1000, 21995324-1000, 21995324+1000)# NM_000077短，NM_001363763长

clinical = fread("clinical.txt")
clinical = clinical[clinical$subtype == "ER-HER2-",]
clinical = clinical[!clinical$SampleID %in% c("PS066_T","PS106_T"),]
clinical = clinical[clinical$Methylation != "",]

cpgmet = as.data.frame(fread("./data/cpgs.met-c10.CDKN2A.tsv"))
cpgmet = cpgmet[,c("chrom","start","end",clinical$SampleID)]
cpgmet$CallRate = rowSums(!is.na(cpgmet[,-c(1,2,3)]))/c(ncol(cpgmet)-3)
cpgmet = cpgmet[cpgmet$CallRate > 0.8,]

subsetDat = cpgmet[cpgmet$start > 21992500-1200 & cpgmet$start < 21992500+1000,]
rownames(subsetDat) = paste(subsetDat$chrom,subsetDat$start,subsetDat$end,sep = "_")
subsetDat = as.data.frame(t(subsetDat[,-c(1,2,3)]))
subsetDat$SampleID = rownames(subsetDat)
subsetDat = merge(clinical,subsetDat,by = "SampleID")


CDKN2A_sig_DMR = as.data.frame(fread("table.3C.DMRs.ER-HER2-.txt"))
CDKN2A_sig_DMR = CDKN2A_sig_DMR[CDKN2A_sig_DMR$SYMBOL == "CDKN2A" & CDKN2A_sig_DMR$pValues < 0.05,]
CDKN2A_sig_DMR = CDKN2A_sig_DMR[order(CDKN2A_sig_DMR$log2foldChange,decreasing = F),]
# CDKN2A_sig_DMR = CDKN2A_sig_DMR[1:5,"key"]

# DMR annotation
dmrs <- as.data.frame(fread("./data/table.3-2.wgbstools_block.anno.txt"))
dmrs = dmrs[dmrs$seqnames == "chr9",]
dmrs = dmrs[dmrs$start > 21963755 & dmrs$end < 21999264,]
dmrs = merge(dmrs,CDKN2A_sig_DMR,by = "key")
colnames(dmrs)[1:4] = c("key","chr","start","end")


clinical.pCR = clinical[clinical$pCR_consensus == "pCR",]
clinical.RD = clinical[clinical$pCR_consensus == "RD",]

ZoomIN_dmr = cpgmet$start > 21996049 & cpgmet$end < 21996912
cpgmet.pCR = cpgmet[,c("chrom","start","end",clinical.pCR$SampleID)]#ZoomIN_dmr
cpgmet.RD = cpgmet[,c("chrom","start","end",clinical.RD$SampleID)]#ZoomIN_dmr

cpgmet.pCR$CallRate = rowSums(!is.na(cpgmet.pCR[,-c(1,2,3)]))/ncol(cpgmet.pCR)
cpgmet.pCR$MeanBeta = rowMeans(cpgmet.pCR[,clinical.pCR$SampleID],na.rm = T)

pCR_id = clinical.pCR$SampleID
RD_id = clinical.RD$SampleID


# Create the plot
dmrs.highlight = dmrs#[abs(dmrs$mean_methylation_difference) > 0.1,]
dmrs.highlight$pCR_RD_meanDiff = ifelse(dmrs.highlight$log2foldChange > 0,"pCR_higher","RD_higher")


cdkn2a_line = ggplot() + 
  geom_point(data = cpgmet.RD, aes(x = start, y = rowMeans(cpgmet.RD[,RD_id],na.rm = T)), color = '#fc0000', size = 2)+
  geom_point(data = cpgmet.pCR, aes(x = start, y = rowMeans(cpgmet.pCR[,pCR_id],na.rm = T)), color = 'lightblue', size = 2)+
  geom_line(data = cpgmet.RD, aes(x = start, y = rowMeans(cpgmet.RD[,RD_id],na.rm = T)), color = '#ff7878', size = 1.5) + 
  geom_line(data = cpgmet.pCR, aes(x = start, y = rowMeans(cpgmet.pCR[,pCR_id],na.rm = T)), color = 'lightblue', size = 1.5) +
  geom_segment(data = dmrs.highlight, aes(x = start, xend = end, y = -0.1, yend = -0.1,
                                color = pCR_RD_meanDiff),size = 2) + scale_color_manual(values = c('lightblue','#ff7878')) + 
  theme_classic(base_size = 30) + ylab("Beta") + 
  geom_vline(xintercept = c(tss), color = 'lightgray', linetype = "dashed", size = 1) +
  geom_vline(xintercept = c(21996194), color = '#C07ABB', linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, color = '#000000') 
cdkn2a_line

pdf("Figure.3D-part1.CDKN2A.all.CpG.pdf",width = 22,height = 4.5)
plot(cdkn2a_line)
dev.off()





### CpG level ranksum test
cpgmet = as.data.frame(fread("./data/cpgs.met-c10.CDKN2A.tsv"))
cpgmet = cpgmet[,c("chrom","start","end",clinical$SampleID)]
cpgmet$CallRate = rowSums(!is.na(cpgmet[,-c(1,2,3)]))/c(ncol(cpgmet)-3)
cpgmet = cpgmet[cpgmet$CallRate > 0.8,]
rownames(cpgmet) = paste(cpgmet$chrom,cpgmet$start,cpgmet$end,sep = "_")
cpgmet = cpgmet[,clinical$SampleID]

bc_exp_raw = cpgmet[,clinical$SampleID]
RNA_samples = colnames(bc_exp_raw)

clinical.sub = clinical
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

fold_calculator = "byMean"
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


rankSum_CpG = outRst
rankSum_CpG$start = as.numeric(str_split_fixed(rownames(rankSum_CpG),"_",3)[,2])
rankSum_CpG = rankSum_CpG[order(rankSum_CpG$start),]
rankSum_CpG$color = "gray60"
rankSum_CpG.add = rankSum_CpG[c(1,nrow(rankSum_CpG)),]
rankSum_CpG.add$color = "black"

rankSum_CpG = rankSum_CpG[rankSum_CpG$pValues < 0.05,]
rankSum_CpG = rbind(rankSum_CpG,rankSum_CpG.add)
rankSum_CpG$yvalue = 1
rankSum_CpG$logP = -log10(rankSum_CpG$pValues)

pvalue_plot = ggplot(rankSum_CpG, aes(x=start, y=yvalue,size = logP)) +
  geom_point(col = rankSum_CpG$color) + theme_classic(base_size = 20) + theme(legend.position = "right") + scale_size_continuous(range = c(1, 10))#,size = rankSum_CpG$logP*4+2

pdf("Figure.3D-part2.CDKN2A.all.CpG.Pvalue.pdf",width = 22,height = 2)
plot(pvalue_plot)
dev.off()



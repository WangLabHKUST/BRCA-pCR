library("rstatix")
library("data.table")
library("stringr")
library("dplyr")
library("ggpubr")
library("EnvStats")

#################################################################################
## gene level log2ratio from the somatic segment file (annotated using Funcotator)
#################################################################################
codingGene = read.table("./data/Coding.gene.list.txt")
codingGene = codingGene$V1

### 这个是我call的结果
facets_seg = fread("./data/BRCA_hg38.funcotator.anno.seg")
facets_seg$seg = str_split_fixed(facets_seg$seg,"-VS-",2)[,1]

facets_seg.k<- facets_seg[facets_seg$gene %in% codingGene,]
facets_seg.k$cnlr.median = round(facets_seg.k$cnlr.median,3)
facets_seg.k$seg_length = as.numeric(facets_seg.k$end)-as.numeric(facets_seg.k$start)
facets_seg.k$length_X_mean = as.numeric(facets_seg.k$cnlr.median)*as.numeric(facets_seg.k$seg_length)


facets_seg.k.dat = reshape2::dcast(facets_seg.k,seg~gene,value.var = "length_X_mean",sum)
dim(facets_seg.k.dat)
facets_seg.k.dat[1:10,1:10]

facets_seg.k.dat_melt = reshape2::melt(facets_seg.k.dat)
facets_seg.k.dat_melt$key = paste(facets_seg.k.dat_melt$seg,facets_seg.k.dat_melt$variable,sep = "_")
colnames(facets_seg.k.dat_melt)[3] = "length_X_mean"

facets_seg.k.dat_length = reshape2::dcast(facets_seg.k,seg~gene,value.var = "seg_length",sum)
facets_seg.k.dat_length_melt = reshape2::melt(facets_seg.k.dat_length)
facets_seg.k.dat_length_melt$key = paste(facets_seg.k.dat_length_melt$seg,facets_seg.k.dat_length_melt$variable,sep = "_")
head(facets_seg.k.dat_length_melt)
colnames(facets_seg.k.dat_length_melt)[3] = "seg_length"


final_dat2 = merge(facets_seg.k.dat_melt,facets_seg.k.dat_length_melt,by="key")
head(final_dat2)
final_dat2$final_seg_mean = final_dat2$length_X_mean/final_dat2$seg_length
final_dat2 = final_dat2[final_dat2$seg_length!=0,]
dim(final_dat2)

log2ratio2 = as.data.frame(reshape2::dcast(final_dat2,variable.x~seg.x,value.var = "final_seg_mean"))
samplenum = ncol(log2ratio2)-1
callrate = rowSums(!is.na(log2ratio2[,-1]))/samplenum
log2ratio_dat = log2ratio2[callrate > 0.8,]
rownames(log2ratio_dat) = log2ratio_dat$variable.x
log2ratio_dat = as.data.frame(t(log2ratio_dat[,-1]))
log2ratio_dat[1:10,1:10]
dim(log2ratio_dat)



# focus on these cytoband as these were detected by GISTIC2
Final_gene = c()
group_list = c("ER-HER2-_0","ER-HER2-_1","ER-HER2+_0","ER-HER2+_1","ER+HER2+_0","ER+HER2+_1")
for (i in 1:length(group_list)){
  type_list = c("amp","del")
  for (j in 1:length(type_list)){
    group_list.tmp = group_list[i]
    type_list.tmp = type_list[j]
    
    amp = fread(paste("./data/GISTIC2_bySubtype_byPCR/",group_list.tmp,"/",type_list.tmp,"_genes.conf_90.txt",sep = ""),header = T)
    amp = as.data.frame(amp[-c(1,2,3),-1])
    amp = amp[,-ncol(amp)]
    amp$Gene = "Gene"
    amp_long = reshape::melt(amp,id = "Gene")
    amp_long = amp_long[amp_long$value!="",]
    amp_long$Group = group_list.tmp
    amp_long$Type = type_list.tmp
    
    Final_gene = rbind(Final_gene,amp_long)
  }
}
Final_gene.keep = Final_gene[,-1]
colnames(Final_gene.keep) = c("cytoband","gene","subtype_outcome","type")
unique(Final_gene.keep$cytoband)



mapping_dat = fread("./data/GISTIC2_bySubtype_byPCR/ER-HER2-_0/all_data_by_genes.txt")[,c(1,3)]
colnames(mapping_dat) = c("geneid","cytoband")
mapping_dat$geneid = str_split_fixed(mapping_dat$geneid,"[|]",2)[,1]
mapping_dat = mapping_dat[mapping_dat$cytoband %in% unique(Final_gene.keep$cytoband),]
unique(mapping_dat$cytoband)



keep_gene = intersect(unique(Final_gene.keep$gene),colnames(log2ratio_dat))
length(keep_gene)


log2ratio_dat.keep = log2ratio_dat[,c(keep_gene)]
log2ratio_dat.keep[1:10,1:10]
log2Data = as.data.frame(t(log2ratio_dat.keep))
dim(log2Data)
log2Data[1:10,1:10]


sum_dat = data.frame(cytoband = NA)
geneNum = data.frame(cytoband = NA)
for (i in 1:ncol(log2Data)){
  sid = colnames(log2Data)[i]
  log2Data.tmp = data.frame(geneid = rownames(log2Data),
                            log2ratio = log2Data[,i])
  intersect(mapping_dat$geneid,log2Data.tmp$geneid)
  
  log2Data.tmp.merge = merge(mapping_dat,log2Data.tmp,by = "geneid",all.x = T)#
  log2Data.tmp.merge = log2Data.tmp.merge[!is.na(log2Data.tmp.merge$log2ratio),]
  
  # log2Data.tmp.merge$log2ratio[is.na(log2Data.tmp.merge$log2ratio)] = 0
  
  log2Data.tmp.merge = log2Data.tmp.merge[order(log2Data.tmp.merge$cytoband),]
  length(unique(log2Data.tmp.merge$cytoband))
  
  
  # take the median log2ratio of genes locating in a cytoband to represent a cytoband
  median_by_group.tmp <- log2Data.tmp.merge %>%
    group_by(cytoband) %>%
    summarise(median_value = median(log2ratio)) %>% as.data.frame()
  colnames(median_by_group.tmp)[2] = sid
  sum_dat = merge(sum_dat,median_by_group.tmp,all = T)
  
  geneNuminCytoband.tmp = reshape2::dcast(log2Data.tmp.merge,cytoband~.,value.var = "geneid")
  colnames(geneNuminCytoband.tmp)[2] = sid
  geneNum = merge(geneNum,geneNuminCytoband.tmp,all = T)
}

sum_dat.o = sum_dat[!is.na(sum_dat$cytoband),]
dim(sum_dat.o)

### do rank sum test for each of the subtype
subtypes = c("ER-HER2+","ER-HER2-","ER+HER2+")
res_dat = c()
for (k in 1:length(subtypes)){
  subtype.tmp = subtypes[k]
  clinical = fread("./data/clinical.txt")
  clinical.sub = clinical[clinical$WES!="" & !is.na(clinical$Purity),]
  
  clinical.sub = clinical.sub[clinical.sub$subtype==subtype.tmp,]
  dim(clinical.sub)
  mapping = clinical.sub[,c("SampleID","pCR_consensus","Purity")]
  colnames(mapping)[2] = "pCR"
  
  sum_dat.o.tmp = as.data.frame(sum_dat.o)
  rownames(sum_dat.o.tmp) = sum_dat.o.tmp$cytoband
  sum_dat.o.tmp = as.data.frame(t(sum_dat.o.tmp[,-1]))
  sum_dat.o.tmp$SampleID = rownames(sum_dat.o.tmp)
  
  mergedat = merge(sum_dat.o.tmp,mapping,by = "SampleID")
  mergedat[1:10,1:10]
  dim(mergedat)
  longdat = reshape2::melt(mergedat)
  
  
  head(longdat)
  pair_test = longdat %>% group_by(variable) %>% pairwise_wilcox_test(value ~ pCR, p.adjust.method = "bonferroni")
  pair_test = as.data.frame(pair_test)
  pair_test$p.adj =p.adjust(pair_test$p)
  pair_test = pair_test[order(pair_test$p,decreasing = F),]
  pair_test$subtype = subtype.tmp
  
  res_dat = rbind(res_dat,pair_test)
}

res_dat = res_dat[order(res_dat$p),]
head(res_dat,10)



clinical = fread("./data/clinical.txt")
clinical.sub = clinical[clinical$WES!="" & !is.na(clinical$Purity),]
clinical.sub = clinical.sub[clinical.sub$subtype=="ER-HER2-",]
dim(clinical.sub)
mapping = clinical.sub[,c("SampleID","pCR_consensus","Purity","subtype")]
colnames(mapping)[2] = "pCR"

sum_dat.o.tmp = as.data.frame(sum_dat.o)
rownames(sum_dat.o.tmp) = sum_dat.o.tmp$cytoband
sum_dat.o.tmp = as.data.frame(t(sum_dat.o.tmp[,-1]))
sum_dat.o.tmp$SampleID = rownames(sum_dat.o.tmp)

mergedat = merge(sum_dat.o.tmp,mapping,by = "SampleID")


p <- ggboxplot(mergedat, x = "pCR", y = "6q27",
               color = "pCR", palette = "nejm",
               add = "jitter") + stat_compare_means(label.x.npc = 0,label.y.npc = 0.9,label = "p.format") + stat_n_text() + theme(legend.position="none") + theme_classic(base_size = 20)+ scale_color_manual(values = c("#FE818B","#A4E2EE"))
p

pdf("Figure.2E.6q27.pdf",width = 4,height = 4)
plot(p)
dev.off()




rm(list = ls())

library("pheatmap")
library("ggrepel")
library("ggplotify")
library("data.table")

################################################ 
################################################
clinical = as.data.frame(fread("./data/clinical.txt"))
clinical = clinical[clinical$WES!="",]
clinical = clinical[!is.na(clinical$Purity),]
clinical$pCR = clinical$pCR_consensus
dim(clinical)
Assignment.SBS = clinical[,grepl("SBS",colnames(clinical))]
Assignment.SBS$SampleID = clinical$SampleID
rownames(Assignment.SBS) = Assignment.SBS$SampleID

##################
######
library("RColorBrewer")
display.brewer.all(n=11,type="div"); title(main = "Divergent color palette")
display.brewer.all(n=9,type=c("seq")); title(main = "Sequential color palette")
YlGnBu_n3<- brewer.pal(5, "YlGnBu")
mycolors_positive <- colorRampPalette(YlGnBu_n3)(500)


clinical.d_keep = clinical
dim(clinical.d_keep)

dat2 = as.data.frame(t(Assignment.SBS[,grepl("f_",colnames(Assignment.SBS))]))
Keep_sig = rownames(dat2)

dat2_keep = dat2[c(Keep_sig),]
dat2_com = dat2[setdiff(rownames(dat2),Keep_sig),]
dat2_com = colSums(dat2_com)

dat2_forMelt = rbind(dat2_keep,dat2_com)
rownames(dat2_forMelt)[nrow(dat2_forMelt)] = "Others"
dat2 = dat2_forMelt[-nrow(dat2_forMelt),]

res2.anno<- clinical[,c("SampleID","PIK3CA","TP53","subtype","pCR")]
rownames(res2.anno)<- res2.anno$SampleID
res2.anno<- res2.anno[,-1]
dim(res2.anno)# 344  15


ann_colors<- list(PIK3CA=c(wildtype="#F6F6F6",mutant="#636363"),
                  TP53=c(wildtype="#F6F6F6",mutant="#636363"),
                  
                  pCR=c("pCR"="#ADD8E6","RD"="#FF7878"),
                  subtype=c("ER-HER2-"="#B6B6B6","ER-HER2+"="#AEEBC5","ER+HER2+"="#FED6FA"),
                  "T"=c(M="#59A7D1",F="#FFE2DF"))


dat2.tmp<- as.data.frame(apply(dat2, 2, as.numeric))
rownames(dat2.tmp)<- rownames(dat2)
dim(dat2.tmp)

#########################################################################################################
####### ABSOLUTE activity
#########################################################################################################
Assignment_abs = clinical[,grepl("SBS",colnames(clinical))]
Assignment_abs = Assignment_abs[,!grepl("f_SBS",colnames(Assignment_abs))]

rownames(Assignment_abs) = clinical$SampleID

Assignment_abs_com = Assignment_abs


#########################################################################################################
annotation_col=res2.anno[colnames(dat2.tmp),]
dim(annotation_col)

dim(res2.anno)# 344  15
dim(dat2.tmp)# 30 324

dat2.tmp.plot.all = dat2.tmp

rele_simple = dat2.tmp.plot.all
obj_combine <- pheatmap(
  rele_simple,
  silent=F,
  cluster_rows=T,
  cluster_cols=T,
  annotation_col=annotation_col,
  annotation_colors = ann_colors,
  show_rownames=T,
  show_colnames=TRUE,
  col=mycolors_positive,
  legend=F,
  main="mutational signatures contribution",
  number_color="black",
  border_color="gray50",
  treeheight_row=15,
  treeheight_col=80
)


Sample_Order = colnames(rele_simple)[obj_combine$tree_col$order]

###############
###
library("reshape2")
library("ggplot2")

sig_rank<- rowSums(dat2.tmp.plot.all)
sig_rank<- names(sig_rank[order(sig_rank,decreasing = T)])
dat2_forMelt<- dat2.tmp.plot.all


dat2_forMelt$Signature<- rownames(dat2_forMelt)
m_contribution<- melt(dat2_forMelt)
colnames(m_contribution) = c("Signature","Sample", "Contribution")

m_contribution$Sample<- as.character(m_contribution$Sample)
sample_list<- unique(as.character(m_contribution$Sample))
m_contribution$Sample<- factor(m_contribution$Sample,levels = colnames(dat2.tmp.plot.all))
m_contribution$Signature<- factor(m_contribution$Signature,levels = c(Keep_sig,"Others"))
head(m_contribution)


Tri_abs<- Assignment_abs_com

cutoof_for_scaling = 1000
Tri_abs$scaling<- rowSums(Tri_abs) > cutoof_for_scaling
Tri_abs$k <- cutoof_for_scaling/rowSums(Tri_abs)

larger400 =  rownames(Tri_abs)[rowSums(Tri_abs) > cutoof_for_scaling]


index<- Tri_abs$scaling==TRUE
to_adjust = which(grepl("SBS",colnames(Tri_abs)))
for (J in to_adjust){
  Tri_abs[,J][index] = Tri_abs[,J][index]*Tri_abs[,"k"][index]
}

scaling_s<- rownames(Tri_abs)[Tri_abs$scaling]
Tri_abs<- Tri_abs[,grepl("SBS",colnames(Tri_abs))]

Tri_abs$Sample<- rownames(Tri_abs)
head(Tri_abs)
Tri_abs_melt<- melt(Tri_abs)
colnames(Tri_abs_melt) = c("Sample","Signature","Contribution")
Tri_abs_melt$Sample<- factor(Tri_abs_melt$Sample,levels = colnames(dat2.tmp.plot.all))
head(Tri_abs_melt)
Tri_abs_melt$label<- NA
To_Add_label<- unique(c(scaling_s,larger400))


sim.dat.tmp = clinical[,c("SampleID","Total_SNV")]
colnames(sim.dat.tmp)<- c("Sample","freq")


Tri_abs_melt<- merge(Tri_abs_melt,sim.dat.tmp,by="Sample",all.x = T)
head(Tri_abs_melt,10)


index_2<- Tri_abs_melt$Sample %in% gsub("-","-",To_Add_label)
Tri_abs_melt$label[index_2]<- as.character(paste(Tri_abs_melt$Sample," (n=",Tri_abs_melt$freq,")",sep = "")[index_2])
Tri_abs_melt$label[Tri_abs_melt$Signature!="SBS1"]<- NA


set.seed(2021)

#############################################################################
## just for y xias lable
#############################################################################
Final_color = c("#3DA9F9","#007EBD","#f22020","gray80","green","skyblue","#946aa2","#5d4c86")

head(Tri_abs_melt,10)
Tri_abs_melt$Sample = factor(Tri_abs_melt$Sample,levels = Sample_Order)
plot1 = ggplot(Tri_abs_melt, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = factor(Sample),label = label)) + geom_bar(stat = "identity",width = 1) + labs(x = "", y = "Absolute contribution") + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) + theme()+theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())+ scale_fill_manual(values=Final_color)+ theme(legend.position="bottom")+  geom_text_repel(
  force_pull   = 2,
  nudge_y      = 300,
  direction    = "x",
  hjust         = 1,
  angle        = 0,
  arrow = arrow(length = unit(0.01, "npc")),
  box.padding = 1,
)+ theme(legend.position=c(-0.5,0.8))+guides(fill=guide_legend(ncol=1))


reverse_abs= F
if (reverse_abs == T){
  Tri_abs_melt$Signature = factor(Tri_abs_melt$Signature,levels = rev(c("SBS1","SBS5","SBS3", "SBS8", "SBS30", "Others")))
}


m_contribution_rev= m_contribution
reverse_abs = F
if (reverse_abs == T){
  m_contribution_rev$Signature = factor(m_contribution_rev$Signature,levels = rev(c("rele_SBS1","rele_SBS5","rele_SBS3", "rele_SBS8", "rele_SBS30", "Others" )))
}


Col_Toused_updated = Col_Toused[1:length(unique(m_contribution_rev$Signature))]
Col_Toused_updated = c("gray80","#DD4C8F","#AC4425","#6D95A8","#7BCCEE","#007EBD","#946aa2","#5d4c86")

m_contribution_rev$Sample = factor(m_contribution_rev$Sample,levels = Sample_Order)
plot2 = ggplot(m_contribution_rev, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = Sample)) + geom_bar(position = "fill", stat = "identity",width = 1) + labs(x = "", y = "Relative contribution") + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size=4))+ scale_fill_manual(values=c(Final_color))+ theme(legend.position=c(-0.5,0.6),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y = element_blank())+guides(fill=guide_legend(ncol=2))+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

plot1 = ggplot(Tri_abs_melt, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = factor(Sample),label = label)) + geom_bar(stat = "identity",width = 1) + labs(x = "", y = "Absolute contribution") + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) + theme()+theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())+ scale_fill_manual(values=Final_color)+ theme(legend.position="bottom")+  geom_text_repel(
  force_pull   = 2,
  nudge_y      = 300,
  direction    = "x",
  hjust         = 1,
  angle        = 0,
  arrow = arrow(length = unit(0.01, "npc")),
  box.padding = 1,
)+ theme(legend.position=c(-0.1,0.8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y = element_blank())+guides(fill=guide_legend(ncol=1))



as.ggplot(obj_combine$gtable$grobs[[7]])


Plot_samples = unique(m_contribution_rev$Sample)

pdf(paste("Figure.2A.pdf",sep = ""),width = 0.12*length(Plot_samples),height = 6)
cowplot::plot_grid(as.ggplot(obj_combine$gtable$grobs[[7]]),
                   plot1,
                   plot2,
                   ncol=1,align = "v",rel_heights = c(0.3,
                                                      #0.3,
                                                      # 0.02,
                                                      #0.05,
                                                      0.5,
                                                      0.5))#
dev.off()


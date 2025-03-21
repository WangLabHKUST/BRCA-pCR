rm(list = ls())

library("data.table")
library("pheatmap")

mycolors_positive <- c("white","darkred","#8FBC94","#4FB0C6","#D7463F","#f8ca00","#4f953b","#99CCCC","#4F86C6","#C65146","#EC6A5C","#e97f02","#6E7783","#77AAAD","#548687","#3a5134")

WES_clinical = fread("./data/clinical.txt")
WES_clinical = WES_clinical[WES_clinical$WES!="",]
WES_clinical = WES_clinical[!is.na(WES_clinical$Purity),]
dim(WES_clinical)


onco = read.table("./data/table.2B.input.txt",header = T,sep = "\t")

# Type_arrange
Type_arrange = c("Amp","Gain","Neu","Loss","Deletion")
mycolor1<- rev(c("#2896AF","#A7DDEA","#ECECEC","#EFB9C2","#CF4C35"))
fix_col_cnv = data.frame(Type= Type_arrange,col=rev(mycolor1[1:length(Type_arrange)]),number_in_grid=c(seq(length(Type_arrange)-1,0)))


vc_cols = RColorBrewer::brewer.pal(n = 9, name = 'Paired')
vc_cols[length(vc_cols)] = "white"
snv_type = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  "wt"
)
fix_col_snv = data.frame(Type = snv_type, col = vc_cols, number_in_grid = seq(5,4+length(snv_type)))

fix_col = rbind(fix_col_cnv,fix_col_snv)

Type_arrange_both = sort(unique(as.vector(as.matrix(onco))))


for (j in 1:length(Type_arrange_both)){
  tmp<- Type_arrange_both[j]
  onco[onco==tmp]<- fix_col[fix_col$Type==tmp,"number_in_grid"]
}


onco = as.data.frame(onco)

for (k in 1:ncol(onco)){
  onco[,k]<- as.numeric(as.character(onco[,k]))
}
onco = as.matrix(onco)


table(WES_clinical$subtype)
WES_clinical$pCR = WES_clinical$pCR_consensus
clinical_col<- as.data.frame(WES_clinical[,c("subtype","pCR")])
rownames(clinical_col) = WES_clinical$SampleID

ann_colors<- list(pCR=c("pCR"="#ADD8E6","RD"="#FF7878"),
                  subtype=c("ER-HER2-"="#B6B6B6","ER-HER2+"="#AEEBC5","ER+HER2+"="#FED6FA"),
                  funcoSEG_RB1=c(Amplication="#632626",Gain = "#F76E11",Neutral="#F6F6F6",Loss="#90E0EF",Deletion="#1C658C"))

clinical_col2 = clinical_col
clinical_col2$SampleID = rownames(clinical_col2)


table(clinical_col$pCR)
clinical_col.pCR = clinical_col[clinical_col$pCR == "pCR",]
onco.pCR = onco[,rownames(clinical_col.pCR)]
onco.pCR = as.data.frame(t(onco.pCR))
onco.pCR = onco.pCR[order(onco.pCR$TP53,onco.pCR$PIK3CA,onco.pCR$GATA3),]
onco.pCR = as.matrix(t(onco.pCR))

pheatmap(
  filename = paste("Figure.2B_pCR.pdf",sep = ""),
  width = 9,
  height = 6,
  annotation_colors = ann_colors,
  annotation_col = clinical_col.pCR,
  onco.pCR,
  gaps_row = 8,
  cluster_rows=F,
  cluster_cols=F,
  show_rownames=TRUE,
  show_colnames=T,
  color=fix_col$col,
  # breaks=seq(0,1,length.out=500),
  #cluster_rows=TRUE,
  #cluster_cols=TRUE,
  legend=TRUE,
  legend_breaks=fix_col$number_in_grid,
  legend_labels=fix_col$Type,
  fontsize=10,
  fontsize_col=8,
  fontsize_row=8,
  main=paste("Somatic mutation landscape n=",ncol(onco),sep = ""),
  # display_numbers=F,
  # fontsize_number=4,
  # number_color="black",
  cellwidth = 8,
  cellheight= 12,
  border_color="gray95",
  treeheight_row=10,
  treeheight_col=10
)

###### 
clinical_col.RD = clinical_col[clinical_col$pCR == "RD",]
onco.RD = onco[,rownames(clinical_col.RD)]
onco.RD = as.data.frame(t(onco.RD))
onco.RD = onco.RD[order(onco.RD$TP53,onco.RD$PIK3CA,onco.RD$GATA3),]
onco.RD = as.matrix(t(onco.RD))

pheatmap(
  filename = paste("Figure.2B_RD.pdf",sep = ""),
  width = 9,
  height = 6,
  annotation_colors = ann_colors,
  annotation_col = clinical_col.RD,
  onco.RD,
  gaps_row = 8,
  cluster_rows=F,
  cluster_cols=F,
  show_rownames=TRUE,
  show_colnames=T,
  color=fix_col$col,
  # breaks=seq(0,1,length.out=500),
  #cluster_rows=TRUE,
  #cluster_cols=TRUE,
  legend=TRUE,
  legend_breaks=fix_col$number_in_grid,
  legend_labels=fix_col$Type,
  fontsize=10,
  fontsize_col=8,
  fontsize_row=8,
  main=paste("Somatic mutation landscape n=",ncol(onco),sep = ""),
  # display_numbers=F,
  # fontsize_number=4,
  # number_color="black",
  cellwidth = 8,
  cellheight= 12,
  border_color="gray95",
  treeheight_row=10,
  treeheight_col=10
)







library("fgsea")
library("org.Hs.eg.db")
library("tibble")
library("dplyr")
library("gridExtra")
library("ggsci")
library("pheatmap")


subtypes = c("ER-HER2-","ER-HER2+","ER+HER2+")
pathways.hallmark <- gmtPathways("./data/h.all.v7.2.symbols.gmt")
mergdat = data.frame(pathway = names(pathways.hallmark))

for (i in 1:length(subtypes)){
  subtypes.tmp = subtypes[i]
  res = read.csv(paste("table.4C-1.DEGs-newTPM_",subtypes.tmp,".csv",sep = ""))
  rownames(res) = res$X
  res$row <- rownames(res)
  
  # creating  a named vector [ranked genes]
  head(res)
  res2 = res[order(res$stat),]
  res2$SYMBOL = rownames(res2)
  ranks <- res2$stat
  names(ranks) <- res2$SYMBOL
  
  #Running fgsea algorithm:
  fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark,
                              stats=ranks,
                              eps = 0.0)
  
  
  # Tidy the results:
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) # order by normalized enrichment score (NES)
  fgseaResTidy$leadingEdge_Collapse = NA
  for (i in 1:nrow(fgseaResTidy)){
    fgseaResTidy$leadingEdge_Collapse[i] = paste(fgseaResTidy$leadingEdge[[i]],collapse = ",")
  }
  fgseaResTidy= fgseaResTidy[,setdiff(colnames(fgseaResTidy),"leadingEdge")]
  
  gesres = as.data.frame(fgseaResTidy)[,-8]
  colnames(gesres) = paste(subtypes.tmp,colnames(gesres),sep = ".");colnames(gesres)[1] = "pathway"
  mergdat = merge(mergdat,gesres,by = "pathway")
}



library("stringr")


subtypes = c("ER-HER2-","ER-HER2+","ER+HER2+")
pvalue = paste(subtypes,c(".padj"),sep = "")
nes = paste(subtypes,c(".NES"),sep = "")


cleandat = mergdat[,c("pathway",pvalue,nes)]
cleandat_long = reshape2::melt(cleandat)
cleandat_long$subtype = str_split_fixed(cleandat_long$variable,"[.]",2)[,1]
cleandat_long$datatype = str_split_fixed(cleandat_long$variable,"[.]",2)[,2]


cleandat_long.size = cleandat_long[cleandat_long$datatype == "padj",]
cleandat_long.size$value = -log10(cleandat_long.size$value)
cleandat_long.size$Key = paste(cleandat_long.size$pathway,cleandat_long.size$subtype,sep = ".")

cleandat_long.col = cleandat_long[cleandat_long$datatype == "NES",]
cleandat_long.col$Key = paste(cleandat_long.col$pathway,cleandat_long.col$subtype,sep = ".")
cleandat_long.col = cleandat_long.col[,c("Key","value")]
colnames(cleandat_long.col)[2] = "NES"


df = merge(cleandat_long.size,cleandat_long.col,by = "Key")
df$pathway = gsub("HALLMARK_","",df$pathway)
CL$pathway = gsub("HALLMARK_","",CL$pathway )
df$pathway = factor(df$pathway,levels = as.character(CL$pathway))
df$TrimMinorLog10 = df$value
df$TrimMinorLog10[df$TrimMinorLog10 < -log10(0.05)] = 1
head(df)

df$shape = ifelse(df$value > -log10(0.05),"sig","not_sig")

head(df)
df.matrix = as.data.frame(reshape2::dcast(df,pathway~subtype,value.var = "TrimMinorLog10"))
rownames(df.matrix) = df.matrix$pathway
df.matrix = df.matrix[,-1]

df.matrix.pheat = pheatmap(df.matrix)
pathway_order = rownames(df.matrix)[df.matrix.pheat$tree_row$order]
subtype_order = colnames(df.matrix)[df.matrix.pheat$tree_col$order]

df$pathway = factor(df$pathway,levels = pathway_order)
df$subtype = factor(df$subtype,levels = rev(subtype_order))

plotPathway = ggplot(df, aes(x = pathway, y = subtype)) +
  geom_point(aes(color = NES, size = TrimMinorLog10, shape = shape), alpha = 0.7) +
  #scale_color_manual(values = c("#AA4371", "#E7B800", "#FC4E07")) +
  scale_size(range = c(1, 10)) + # Adjust the range of points size
  theme_bw(base_size = 12)+theme(legend.position = "right") + scale_fill_nejm()+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+scale_color_gradient(low="blue", high="red")+scale_shape_manual(values = c(15,16))
plotPathway


pdf("Figure.4D.Pathway.Across.subtype.pdf",width = 14,height = 3.8)
plot(plotPathway)
dev.off()




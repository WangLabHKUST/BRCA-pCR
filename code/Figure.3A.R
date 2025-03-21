library("cowplot")
library("ggpubr")

clinical = as.data.frame(fread("./data/clinical.txt"))
clinical = clinical[!clinical$SampleID %in% c("PS066_T","PS106_T"),]
clinical = clinical[clinical$Methylation!=0,]
clinical$pCR = factor(clinical$pCR_consensus,levels = c("RD","pCR"))

clinical$subtype = gsub("HER2","",clinical$subtype)
clinical$subtype = gsub("ER","",clinical$subtype)

sumTyps = c("ccre_prom","ccre_enhP","ccre_enhD","ccre_K4m3","ccre_CTCF","meth_wholegenome_CpG_mean")

FigList = list()
for (i in 1:length(sumTyps)){
  sumTyps.tmp = sumTyps[i]
  clinical$FeatureNow = clinical[,sumTyps.tmp]
  meanBetaplot = plot_in_details_tmp(clinical,  "subtype", "pCR","FeatureNow","npg")+ scale_color_manual(values = c("#FE818B","#A4E2EE")) + theme_classic(base_size = 15) + ylab("mean beta") + ggtitle(sumTyps.tmp)  + theme(legend.position = "none")+ ylim(0,1) #+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
  FigList[[sumTyps.tmp]] = meanBetaplot
}


pdf(paste("Figure.3A.pdf",sep = ""),width = 16,height = 3.8)
plot_grid(FigList[[1]],FigList[[2]],FigList[[3]],FigList[[4]],FigList[[5]],FigList[[6]],nrow = 1)#,rel_widths = c(1.7,2.5)
dev.off()



plot_in_details_tmp = function(input_dat, outer_pair, inner_pair, y_val,col_pal){
  inner_pair_name = inner_pair
  colnames(input_dat)[which(colnames(input_dat)==outer_pair)] = "outer_pair"
  colnames(input_dat)[which(colnames(input_dat)==inner_pair)] = "inner_pair"
  colnames(input_dat)[which(colnames(input_dat)==y_val)] = "y_val"
  
  input_dat$outer_pair = factor(input_dat$outer_pair)
  
  p <- ggboxplot(input_dat, x = "outer_pair", y = "y_val",
                 color = "inner_pair", palette = col_pal,add.params = list(size = 2,alpha = 1),
                 add = "") + xlab(outer_pair) + ylab(y_val)+ theme(legend.position="none")
  p
  
  tmp.dat = as.data.frame(table(input_dat$inner_pair,input_dat$outer_pair))
  tmp.dat$com_flag = paste(tmp.dat$Var2,tmp.dat$Var1,sep = "_")
  tmp.inner.notUSED = tmp.dat[tmp.dat$Freq==0,"com_flag"]
  
  input_dat$com_flag = paste(input_dat$outer_pair,input_dat$inner_pair,sep = "_")
  if(length(tmp.inner.notUSED)>0){
    input_dat.tmp = input_dat[!input_dat$com_flag %in% tmp.inner.notUSED,]
  }else{
    input_dat.tmp = input_dat
  }
  stat.test <- input_dat.tmp %>%
    group_by(outer_pair) %>% 
    wilcox_test(y_val ~ inner_pair) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  
  
  stat.test = stat.test[stat.test$p < 0.1,]
  
  
  stat.test <- stat.test %>%
    add_xy_position(x = "outer_pair", dodge = 0.8)
  
  stat.test2 <- input_dat %>%
    wilcox_test(y_val ~ outer_pair, p.adjust.method = "bonferroni")
  if("p.adj.signif" %in% colnames(stat.test2)){
    stat.test2 = stat.test2
  }
  
  stat.test2 <- stat.test2 %>% add_xy_position(x = "outer_pair")
  stat.test2$y.position = stat.test2$y.position*1.2 
  
  p.complex <- p + stat_pvalue_manual(stat.test,  label = "p",step.increase = 0.03, tip.length = 0)
  p.complex
  p.complex <- p.complex +
    stat_pvalue_manual(
      stat.test2,  label = "p", tip.length = 0.02,
      step.increase = 0.3
    ) + ylab(y_val) + guides(color = guide_legend(title = inner_pair_name))
  cat("plotting",inner_pair_name,"\n")
  p.complex
  
}




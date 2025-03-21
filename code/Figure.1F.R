plot_in_details_brca = function(input_dat, outer_pair, inner_pair, y_val,col_pal){
  inner_pair_name = inner_pair
  colnames(input_dat)[which(colnames(input_dat)==outer_pair)] = "outer_pair"
  colnames(input_dat)[which(colnames(input_dat)==inner_pair)] = "inner_pair"
  colnames(input_dat)[which(colnames(input_dat)==y_val)] = "y_val"
  
  input_dat$outer_pair = factor(input_dat$outer_pair)
  
  p <- ggboxplot(input_dat, x = "outer_pair", y = "y_val",
                 color = "inner_pair", palette = col_pal,add.params = list(size = 2,alpha = 0.8),
                 add = "jitter") + xlab(outer_pair) + ylab(y_val)+ theme(legend.position="right")
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
      step.increase = 0.09
    ) + ylab(y_val) + guides(color = guide_legend(title = inner_pair_name))
  cat("plotting",inner_pair_name,"\n")
  p.complex
}


clinical = fread("./data/clinical.txt")
clinical$Freq = 1
clinical$pCR = clinical$pCR_consensus
clinical$pCR = gsub("RD","0RD",clinical$pCR)
clinical$pCR = gsub("pCR","1pCR",clinical$pCR)

clinical$Ki.67

table(clinical$subtype)
clinical$subtype = factor(clinical$subtype,levels = c("ER-HER2-","ER+HER2+","ER-HER2+"))


Ki.67 = plot_in_details_brca(clinical,  "subtype","pCR", "Ki.67","nejm") + scale_color_manual(values = c("#FF6C72","#A1D9E8")) + theme_classic2(base_size = 18) + theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))


pdf("Figure.1F.pdf",width = 5,height = 4.5)
plot(Ki.67)
dev.off()



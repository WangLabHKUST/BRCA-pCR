library("ggsci")

cal_matrix_brac = function(inpu.matrix,col1,col2,return_type,width_scale){
  colnames(inpu.matrix)[which(colnames(inpu.matrix)==col1)]<- "Factor1"
  colnames(inpu.matrix)[which(colnames(inpu.matrix)==col2)]<- "Factor2"
  Input.dat = plyr::count(inpu.matrix,c("Factor1","Factor2"))
  
  if(nrow(Input.dat)>2){
    fisher_input.dat<- as.data.frame(reshape2::dcast(Input.dat,Factor1~Factor2,value.var = "freq")[,-1])
    fisher_input.dat[is.na(fisher_input.dat)]<- 0
    
    if(nrow(fisher_input.dat) >= 2 & ncol(fisher_input.dat) > 2){
      m.fisher<- paste("chisq.test",signif(chisq.test(fisher_input.dat)$p.value,3))
    }else if(nrow(fisher_input.dat) > 2 & ncol(fisher_input.dat) >= 2){
      m.fisher<- paste("chisq.test",signif(chisq.test(fisher_input.dat)$p.value,3))
    }else if(nrow(fisher_input.dat) == 2 & ncol(fisher_input.dat) == 2){
      m.fisher<- paste("fisher",signif(fisher.test(fisher_input.dat)$p.value,3))
    }else{
      m.fisher = NA
    }
    
  }
  
  percentData <- Input.dat %>% group_by(Factor2)  %>%
    mutate(ratio=scales::percent(freq/sum(freq)))
  
  percentData$label<- paste("n=",percentData$freq,"\n(",percentData$ratio,")",sep = "")
  
  percentData.out<- as.data.frame(percentData)
  percentData.out$Key<- paste(percentData.out$Factor2,percentData.out$Factor1,sep = "_")
  percentData.out<- percentData.out[,c("Key","ratio","label")]
  Input.dat$Key<- paste(Input.dat$Factor2,Input.dat$Factor1,sep = "_")
  Input.dat<- merge(Input.dat,percentData.out,by="Key")
  if(nrow(Input.dat)>2){
    Input.dat$fisherp<- m.fisher
  }else{
    Input.dat$fisherp<- NA
  }
  
  stat.dat = Input.dat
  
  stat.dat<- stat.dat[!is.na(stat.dat$Factor2),]
  stat.dat$ratio<- as.numeric(gsub("%","",stat.dat$ratio))
  
  Factor2_Volumn = aggregate(freq ~ Factor2, stat.dat, sum)
  tmp_dat = stat.dat[,c("Factor2","freq")]
  Factor2_Volumn = merge(Factor2_Volumn,tmp_dat,by = "Factor2")
  Factor2_Num = Factor2_Volumn[!duplicated(Factor2_Volumn$Factor2),c(1,2)]
  stat.dat = merge(stat.dat,Factor2_Num,by = "Factor2")
  stat.dat$Factor2 = paste(stat.dat$Factor2,"\n(n=",stat.dat$freq.x,")",sep = "")
  Factor2_Volumn$freq.x = Factor2_Volumn$freq.x/sum(Factor2_Volumn$freq.x)
  
  hide_low_freq = T
  if(hide_low_freq){
    stat.dat$label[stat.dat$ratio < 5] = NA
  }
  
  if(nrow(stat.dat)<10){
    MNP_Race <- ggplot(stat.dat, aes(x = Factor2, y = ratio, fill = Factor1,width = 0.3)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = label),color="white",size=3,
                position = position_stack(vjust = 0.5)) +
      theme(legend.position = "right",
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank()) + scale_fill_npg() + theme_bw() +theme(panel.grid.major.x = element_blank() ,panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())+ggtitle(paste("p=",unique(stat.dat$fisherp,3)))+xlab(col2)+ylab(col1) + guides(fill=guide_legend(title=col1)) 
    
    MNP_Race
  }else{
    MNP_Race <- ggplot(stat.dat, aes(x = Factor2, y = ratio, fill = Factor1,width = 0.3)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = label),color="white",size=3,
                position = position_stack(vjust = 0.5)) +
      theme(legend.position = "right",
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank()) + theme_bw() +theme(panel.grid.major.x = element_blank() ,panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())+ggtitle(paste("p=",unique(stat.dat$fisherp,3)))+xlab(col2)+ylab(col1) + guides(fill=guide_legend(title=col1)) 
    
    MNP_Race
  }
  
  if(return_type == "matrix"){
    return(Input.dat)
  }else if(return_type == "figure"){
    return(MNP_Race)
  }
}


clinical = fread("./data/clinical.txt")
clinical$Freq = 1
clinical$pCR = ifelse(clinical$pCR_consensus == "pCR","1","0")

table(clinical$T)
clinical$T_B = ifelse(clinical$T %in% c("T3","T4"),"T3&T4","T1&T2")
table(clinical$N)
clinical$N_B = ifelse(clinical$N %in% c("N0"),"No","Yes")
# clinical$N_Num = as.numeric(gsub("N","",str_sub(clinical$N,1,2)))

clinical$T_B = factor(clinical$T_B,levels = c("T3&T4","T1&T2"))
T_B <- cal_matrix_brac(clinical,"pCR","T_B","figure",3) + theme_bw(base_size = 18) + scale_fill_manual(values = c("white","#ABABAB"))

N_B <- cal_matrix_brac(clinical,"pCR","N_B","figure",3) + theme_bw(base_size = 18) + scale_fill_manual(values = c("white","#ABABAB"))


pdf("Figure.1E.pdf",width = 7,height = 4)
plot_grid(T_B,N_B)
dev.off()



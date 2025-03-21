library("scales")
library("ggrepel")

pval = c()
OR = c()

clinical = as.data.frame(fread("./data/clinical.txt"))
clinical = clinical[clinical$Histologic_Grade!="",]

clinical$pCR = ifelse(clinical$pCR_consensus == "pCR",1,0)

clinical$T_B = ifelse(clinical$T %in% c("T3","T4"),1,0)
clinical$N_B = ifelse(clinical$N %in% c("N0"),0,1)
clinical$HER2_B = ifelse(clinical$HER2_Status == "Pos",1,0)
clinical$ER_B = ifelse(clinical$ER_Status == "Pos",1,0)
clinical$Histologic_Grade_B = clinical$Histologic_Grade == "III"
clinical$Chemo_B = 1*(clinical$Chemo=='TCb')
clinical$KI67 = clinical$Ki.67
clinical$Age = clinical$Age


table(clinical$HER2_Status,clinical$antiHER2)
clinical$antiHER2 = as.character(clinical$antiHER2)

df = clinical[,c("pCR","T_B","N_B","HER2_B","ER_B","Histologic_Grade_B","Chemo_B","Age","KI67","antiHER2")]
table(df$HER2_B,df$antiHER2)
rownames(df) = clinical$SampleID

#######################################
## figure 1C
#######################################
mpval = c()
mOR = c()

dim(df)

dim(na.omit(df))

myf = 'pCR ~ HER2_B+ER_B+T_B+N_B+KI67+Age+Chemo_B+Histologic_Grade_B+antiHER2'
mylogit2 <- stats::glm(myf, data = df, na.action=na.omit, family = 'binomial', maxit = 1000)
summary(mylogit2)
exp(cbind(OR = coef(mylogit2), confint(mylogit2, level=0.95)))

summary(mylogit2)$coefficients[,4]->mpval
exp(cbind(OR = coef(mylogit2), confint(mylogit2)))->mOR

mpval.df = as.data.frame(mpval)
mpval.df$variable = rownames(mpval.df)
mOR.df = as.data.frame(mOR)
mOR.df$variable = rownames(mOR.df)

mOR.p = merge(mOR.df,mpval.df,by = "variable")

mOR.p.in = mOR.p
mOR.p.in$logP = -log10(mOR.p.in$mpval)
mOR.p.in$color = "gray"
mOR.p.in$color[mOR.p.in$logP > 1 & mOR.p.in$OR < 1] = "#FF6C72"
mOR.p.in$color[mOR.p.in$logP > 1 & mOR.p.in$OR > 1] = "#A1D9E8"
mOR.p.in$label = c("Intercept","Age","antiHER2","TCb","ER+","HER2+","Grade III","KI67","Lymphnode metastasis","T3/T4")#"antiHER2",
mOR.p.in$size = (mOR.p.in$logP) * c(abs(-log2(mOR.p.in$OR)))
mOR.p.in = mOR.p.in[mOR.p.in$variable!="(Intercept)",]

oRplot = ggplot(mOR.p.in, aes(x=OR, y=logP,label = label)) + 
  geom_point(col = mOR.p.in$color,size = 4) + scale_x_continuous(trans='log10') +
  annotation_logticks(sides="b",outside = TRUE)+ theme_classic(base_size = 25) + 
  geom_text_repel(size = 6) + 
  geom_hline(yintercept=1, linetype="dashed", color = "gray",size = 1) + 
  geom_vline(xintercept=1, linetype="dashed", color = "gray",size = 1) + ylab("-Log10 (P value)") + xlab("Odds ratio")

pdf("Figure.1C.pdf",width = 4.4,height = 5)
plot(oRplot)
dev.off()


#######################################
## figure S1D
#######################################
##### univariable logistic regression
# Get the names of predictor variables
predictor_vars <- names(df)[-which(names(df) == "pCR")]

# Create an empty list to store the regression results
results_list <- c()

# Loop through each predictor variable
for (varTmp in predictor_vars) {
  formula <- as.formula(paste("pCR ~", varTmp))
  model <- glm(formula, data = df, family = binomial(link = "logit"))
  # results_list[[varTmp]] <- summary(model)
  
  tidy_model <- tidy(model) %>%
    mutate(predictor = varTmp) %>% 
    mutate(odds_ratio = exp(estimate)) %>%
    # Calculate confidence intervals for odds ratios
    mutate(conf_low = exp(confint.default(model)[, 1]),
           conf_high = exp(confint.default(model)[, 2])) %>% as.data.frame
    
  
  tidy_model$varnow = varTmp
  results_list = rbind(results_list,tidy_model)
}

results_list.o = results_list[results_list$term!="(Intercept)",]
results_list.o$label = c("T3/T4","Lymphnode metastasis","HER2+","ER+","Grade III","TCb","Age","KI67","antiHER2")
results_list.o$logP = -log10(results_list.o$p.value)
results_list.o$color = "gray"
results_list.o$color[results_list.o$logP > -log10(0.05) & results_list.o$odds_ratio < 1] = "#FF6C72"
results_list.o$color[results_list.o$logP > -log10(0.05) & results_list.o$odds_ratio > 1] = "#A1D9E8"

oRplot = ggplot(results_list.o, aes(x=odds_ratio, y=logP,label = label)) + 
  geom_point(col = results_list.o$color,size = 4) + scale_x_continuous(trans='log10') +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides="b",outside = TRUE)+ theme_classic(base_size = 25) + 
  geom_text_repel(size = 6) + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray",size = 1) + 
  geom_vline(xintercept=1, linetype="dashed", color = "gray",size = 1) + ylab("-Log10 (P value)") + xlab("Odds ratio")

pdf("Figure.S1D.univariable.pdf",width = 5,height = 5)
plot(oRplot)
dev.off()









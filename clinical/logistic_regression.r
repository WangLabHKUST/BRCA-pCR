pval = c()
OR = c()

df = read.csv('linear_nona.csv', row.names = 'Patient_ID')
mpval = c()
mOR = c()

myf = 'pCR ~ HER2+ER+T+N+KI67+Age+Chemo+Histologic_Grade+antiHER2'
mylogit2 <- stats::glm(myf, data = df, na.action=na.omit, family = 'binomial', maxit = 1000)
summary(mylogit2)
exp(cbind(OR = coef(mylogit2), confint(mylogit2, level=0.95)))

summary(mylogit2)$coefficients[,4]->mpval
exp(cbind(OR = coef(mylogit2), confint(mylogit2)))[,1]->mOR
write.csv(mpval, 'multi.pval.csv')
write.csv(mOR, 'multi.OR.csv')

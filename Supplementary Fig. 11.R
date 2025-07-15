#Supplementary Fig. 11a#

library(pheatmap)

exp = read.csv("expression.csv",header = T, row.names = 1)

annotation = read.csv("annotation.csv",header = T, row.names = 1)

bk <- c(seq(-3,-0.01,by=0.01),seq(0,3,by=0.01))

pheatmap(total, cluster_cols = FALSE, cluster_rows = FALSE,annotation_row = annotation_row,color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),breaks=bk,show_rownames = TRUE,show_colnames = FALSE

#Supplementary Fig. 11b#

ESCC = read.csv("COCA.csv",header = T, row.names = 1)

boxplot(TMB~MIMER, data = ESCC, col = c('#359fbc', '#e58099'), boxwex = 0.5)

#Supplementary Fig. 11c#

GSEA4.1.0

#Supplementary Fig. 11d#

library(survival)
library(survminer)

data = read.csv("PFS.csv",header = 1)

fit = surv_fit(Surv(time, status) ~ group, data = data)

ggsurvplot(fit, data = data, palette = c("#0054FF", "#FF0000"), conf.int = T, conf.int.style = "step", risk.table = T,pval = T)

#Supplementary Fig. 11e#

fit.cox = coxph(Surv(time,status)~Stage+Smoking+Grade+Gender+Drinking+Age+MIMER, data = ESCC)

x = summary(fit.cox)
pvalue = signif(as.matrix(x$coefficients)[,5],2)
HR = signif(as.amtrix(x@coefficients)[,2],2)
low = signif(x$conf.int[,3],2)
high = signif(x$conf.int[,4],2)

result_mul_cox = data.frame(p.value=pvalue, HR=HR, low=low, high=high, stringsAsFactors = F)

result_mul_cox$'HR(95%CI)' = paste0(HR, "(",low, "-", high, ")", sep = " ")

#Supplementary Fig. 11f#

data = read.csv("OS.csv",header = 1)

fit = surv_fit(Surv(time, status) ~ group, data = data)

ggsurvplot(fit, data = data, palette = c("#0054FF", "#FF0000"), conf.int = T, conf.int.style = "step", risk.table = T,pval = T)

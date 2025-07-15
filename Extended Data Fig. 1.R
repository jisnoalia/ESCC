#Extended Data Fig. 1 -----
#Extended Data Fig. 1 a & b & c-----
library(Seurat)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
library(pheatmap)
library(ggsci)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
############ color choice
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/escc_color.RData")
#endo
load("/workspace/wangxiliang/project/escc_baseline_zhanQiMin/sc_rna/scRNA.RData")
my<-subset(sc,L1_C=='Endothelium');Idents(my)<-my$L4_C

my<-subset(sc,L1_C=='Endothelium');Idents(my)<-my$L4_C
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('Endo_C4_CCL21','Endo_C3_RGCC','Endo_C2_FBLN5','Endo_C1_ACKR1'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="",cellwidth = 60,cellheight = 60)
roeENDO <- roe.sub

load("/workspace/wangxiliang/project/escc_baseline_zhanQiMin/sc_rna/scRNA.RData")
my<-subset(sc,L1_C=='Endothelium');Idents(my)<-my$L4_C
T_plotG <- c('ACKR1','SELP','SELE','CCL23','SEMA3G','FBLN5','GJA5','SERPINE2','PLVAP','RGCC','APLNR','IGFBP5',
             'CCL21','PROX1','PDPN','FLT4')
celltype <- c('Endo_C1_ACKR1','Endo_C2_FBLN5','Endo_C3_RGCC','Endo_C4_CCL21')
ctidx <- c(4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('Endo_C4_CCL21','Endo_C3_RGCC','Endo_C2_FBLN5','Endo_C1_ACKR1'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="",cellwidth = 60,cellheight = 60)
roeENDO <- roe.sub

Idents(my)<-my$tissue
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Tissue<-factor(t$Tissue,levels = c('nLN','pLN','Nor','Adj','Tu'))
ggplot(t, aes(Tissue,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values =c('#DEBA6A','#E63228','#7D5599','#3787BC'))
Idents(my)<-my$patient
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =c('#DEBA6A','#E63228','#7D5599','#3787BC'))

#Fibro -----
load("/workspace/wangxiliang/escc/sc_rna/scRNA.RData")
my<-subset(sc,L1_C=='Fibroblast');Idents(my)<-my$L4_C;my<-subset(my,ident='Smooth muscle cell',invert=T)
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('FB_C6_ACTA2','FB_C5_PDGFRB','FB_C4_APOE','FB_C3_COL1A1',
                   'FB_C2_IGF1','FB_C1_CFD'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")
roe_fb <- roe.sub

Idents(my)<-my$tissue
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Tissue<-factor(t$Tissue,levels = c('nLN','pLN','Nor','Adj','Tu'))
ggplot(t, aes(Tissue,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values =c('#D4C799','#4BAC3E','#BA6C35','#D5A65C','#fdbf73','#800000'))

Idents(my)<-my$patient
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =c('#D4C799','#4BAC3E','#BA6C35','#D5A65C','#fdbf73','#800000'))

# B
load("/workspace/wangxiliang/project/escc_baseline_zhanQiMin/sc_rna/scRNA.RData")
Idents(sc)<-sc$L1_C;my<-subset(sc,ident=c('B','Plasma'));Idents(my)<-my$L4_C
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('Plasma','B_C7_MKI67','B_C6_BCL6','B_C5_ISG15','B_C4_TCL1A',
                   'B_C3_GPR183','B_C2_DUSP4','B_C1_CCR7'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")
roeB <- roe.sub

Idents(my)<-my$tissue
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Tissue<-factor(t$Tissue,levels = c('PBMC','nLN','pLN','Nor','Adj','Tu'))
ggplot(t, aes(Tissue,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values =c('#ED5152','#F06C45','#6BAD9E','#579BC7','#FDB258','#AC8DC3','#E5DD99','#8E6C99'))

Idents(my)<-my$patient
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =c('#ED5152','#F06C45','#6BAD9E','#579BC7','#FDB258','#AC8DC3','#E5DD99','#8E6C99'))

# Myeloid
my<-subset(sc,ident=c('Myeloid'));Idents(my)<-my$L4_C
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('Neutrophil','Mono_C2_CD16','Mono_C1_CD14','Mac_C6_MKI67','Mac_C5_CXCL10','Mac_C4_LYVE1','Mac_C3_C1QC',
                   'Mac_C2_SPP1','Mac_C1_NLRP3','DC_C5_IL3RA','DC_C4_LAMP3','DC_C3_CD1A','DC_C2_CD1C','DC_C1_CLEC9A'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")
roeMYE <- roe.sub

Idents(my)<-my$tissue
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Tissue<-factor(t$Tissue,levels = c('PBMC','nLN','pLN','Nor','Adj','Tu'))
ggplot(t, aes(Tissue,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values =c('#EA4243','#C79B7D','#3989b8','#752773FF','#DAA520',
                                                                '#F88D8D','#c9c193','#9FD18F','#F06061','#48D1CC',
                                                                '#3d804e','#8A64AE','#b7a3c4','#FE9424'))

Idents(my)<-my$patient
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =c('#EA4243','#C79B7D','#3989b8','#752773FF','#DAA520',
                                                                 '#F88D8D','#c9c193','#9FD18F','#F06061','#48D1CC',
                                                                 '#3d804e','#8A64AE','#b7a3c4','#FE9424'))

# CD4
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD4'));Idents(my)<-my$L4_C
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('CD4_C9_IFNG','CD4_C8_ITGB1','CD4_C7_OX40','CD4_C6_CD25','CD4_C5_RTKN2','CD4_C4_GPR183','CD4_C3_ANXA1',
                   'CD4_C2_TCF7','CD4_C13_ISG15','CD4_C12_NKG7','CD4_C11_IL17A','CD4_C10_CXCR5','CD4_C1_CCR7'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")
roeCD4T <- roe.sub

Idents(my)<-my$tissue
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Tissue<-factor(t$Tissue,levels = c('PBMC','nLN','pLN','Nor','Adj','Tu'))
ggplot(t, aes(Tissue,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values =c('#6B3F99','#40E0D0','#999C63','#B19999','#8B4513',
                                                                '#F0E084','#267DB1','#F7861D','#E73335','#A4D880',
                                                                '#20B2AA','#CD853F','#69BA53'))

Idents(my)<-my$patient
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =c('#6B3F99','#40E0D0','#999C63','#B19999','#8B4513',
                                                                 '#F0E084','#267DB1','#F7861D','#E73335','#A4D880',
                                                                 '#20B2AA','#CD853F','#69BA53'))

# CD8
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD8'));Idents(my)<-my$L4_C
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('CD8_C9_KLRC2','CD8_C8_KIR','CD8_C7_CX3CR1','CD8_C6_CD39','CD8_C5_GZMK','CD8_C4_ANXA1','CD8_C3_TCF7','CD8_C2_IL7R',
                   'CD8_C12_CD161','CD8_C11_ISG15','CD8_C10_CD16','CD8_C1_CCR7'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")
roecd8 <- roe.sub

Idents(my)<-my$tissue
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Tissue<-factor(t$Tissue,levels = c('PBMC','nLN','pLN','Nor','Adj','Tu'))
ggplot(t, aes(Tissue,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values =c('#4169E1','#9F8399','#66A5CD','#7CB999','#77C25F','#A6CEE3',
                                                                '#F37F4F','#B79BCA','#5AA1A3','#C2B099','#CC934F','#E31F1E'))

Idents(my)<-my$patient
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =c('#4169E1','#9F8399','#66A5CD','#7CB999','#77C25F','#A6CEE3',
                                                                 '#F37F4F','#B79BCA','#5AA1A3','#C2B099','#CC934F','#E31F1E'))
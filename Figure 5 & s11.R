#Figure 5-----
#Figure 5 b c h -----
load('/workspace/duxiao/ESCC/spatial/spatial_pathology/sp_correct.RData')
gene1<-c('CXCL13','ACP5','LAG3','PHLDA1','HAVCR2','RGS2','PLPP1','RHOB','SNX9','CCL5','CD8A','CD3D') # PDCD1+CD8
gene2<-c('FOXP3','TIGIT','BATF','TNFRSF18','TNFRSF4','TNFRSF9','IL32','CD4','IL10','IL2RA') #Treg
gene3<-c('SPP1','APOC1','MMP12','MMP9','FBP1','APOE','CTSB','CD68','CCL3','TYROBP') # SPP1+MAC
gene4<-c('PLVAP','COL4A1','COL4A2','HSPG2','VWF','IGFBP7','PECAM1','SERPINE1','SPARC','INSR') # PLVAP+Endo
gene5<-c('COL1A1','COL3A1','COL1A2','SPARC','FN1','POSTN','CST1','MMP11','CTHRC1','COL6A3') # CAF
Idents(sp)<-sp$development;p<-subset(sp,ident=c('SD&CA','ICA','MCA'))
my_levels<-c('SD&CA','ICA','MCA')
Idents(p)<-factor(Idents(p),levels = my_levels)
p<-AddModuleScore(p,features = list(gene3),name = 'SPP1.TAM');p<-AddModuleScore(p,features = list(gene1),name = 'PDCD1.CD8')
p<-AddModuleScore(p,features = list(gene2),name = 'OX40.Treg');p<-AddModuleScore(p,features = list(gene4),name = 'PLVAP.Endo');p<-AddModuleScore(p,features = list(gene5),name = 'POSTN.CAF')
fea<-c('SPP1.TAM1','PDCD1.CD81','OX40.Treg1','PLVAP.Endo1','POSTN.CAF1')
myd=as.data.frame(scale(p@meta.data[,fea]))
myd=cbind(myd,p@meta.data[,c('file','development','tissue','metastasis')])
colnames(myd)<-c('SPP1+TAM_score','PDCD1+CD8_score','OX40+Treg_score','PLVAP+Endo_score','POSTN+CAF_score','file','Development','tissue','metastasis')
s1<-myd[myd$Development=='SD&CA',];s1<-aggregate(s1[,1:5],by=list(Group=s1$file),mean);colnames(s1)[1]<-'file';s1$Development<-'SD&CA'
s2<-myd[myd$Development=='ICA',];s2<-aggregate(s2[,1:5],by=list(Group=s2$file),mean);colnames(s2)[1]<-'file';s2$Development<-'ICA'
s3<-myd[myd$Development=='MCA',];s3<-aggregate(s3[,1:5],by=list(Group=s3$file),mean);colnames(s3)[1]<-'file';s3$Development<-'MCA'
new<-rbind(s1,s2);new<-rbind(new,s3)
mycolor<-c("#66A61E","#E6AB02","#A6761D")

my_comparison<-list(c('SD&CA','MCA'),c('ICA','MCA'));ggboxplot(new,x="Development",y="PDCD1+CD8_score",color = "Development",palette = mycolor,add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(comparisons = my_comparison)
my_comparison<-list(c('SD&CA','MCA'),c('ICA','MCA'));ggboxplot(new,x="Development",y="OX40+Treg_score",color = "Development",palette = mycolor,add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(comparisons = my_comparison)
my_comparison<-list(c('SD&CA','MCA'),c('ICA','MCA'));ggboxplot(new,x="Development",y="SPP1+TAM_score",color = "Development",palette = mycolor,add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(comparisons = my_comparison)



#Figure 5 d e -----
library("Startrac")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
tcr_meata <- read.csv("tcr_meta_filter.csv",header = T)
cd8 <- tcr_use[tcr_use$stype == 'CD8',]
out1 <- Startrac.run(cd8, proj="ESCC",verbose=F)
level_cd8 <- c("CD8_C1_CCR7","CD8_C2_IL7R","CD8_C3_TCF7","CD8_C4_ANXA1" ,"CD8_C5_GZMK" , "CD8_C6_PDCD1"  ,
              "CD8_C7_CX3CR1","CD8_C8_KIR", "CD8_C9_KLRC2", "CD8_C10_CD16" ,"CD8_C11_ISG15","CD8_C12_CD161"
)
out2 = as.data.table(out1@cluster.sig.data)[aid!=out1@proj,][order(majorCluster),]
out3 =out2[out2$index == "expa",]
out3$majorCluster <- factor(out3$majorCluster,
                            levels = level_cd8 )
library(ggpubr)
my_comparisons = list(
  c("CD8_C4_ANXA1","CD8_C6_PDCD1")
)
ggboxplot(out3,
          x="majorCluster",y="value",palette = col_cd8$color,
          color = "majorCluster", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = my_comparisons,paired = T,method ="t.test")+
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'none',
  )
ggboxplot(as.data.table(obj@pIndex.sig.migr)[aid!=obj@proj,][order(majorCluster),],
          x="majorCluster",y="value",palette = col_cd8$color,
          color = "majorCluster", add = "point", outlier.colour=NULL) +
  #                  geom_signif(comparisons = list(
  #                     c("CD4_C1_CCR7","CD4_C7_OX40")
  #                  )
  #                  ,test= wilcox.test,
  #             step_increase = 0.1,
  # map_signif_level = function(p) sprintf("p = %.2g", p)
  #             )+
  stat_compare_means(comparisons = my_comparisons,
                     paired = T,method ="t.test")+
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'none')  

#CD4T
cd4 <- tcr_use[tcr_use$stype == 'CD4',]
out1 <- Startrac.run(cd4, proj="ESCC",verbose=F)
level_cd4 <- c("CD4_C1_CCR7" ,"CD4_C2_TCF7" ,"CD4_C3_ANXA1","CD4_C4_GPR183","CD4_C5_RTKN2" ,
               "CD4_C6_CD25","CD4_C7_OX40", "CD4_C8_ITGB1","CD4_C9_IFNG" , "CD4_C10_CXCR5" ,
               "CD4_C11_IL17A" ,"CD4_C12_NKG7" , "CD4_C13_ISG15"
)
out2 = as.data.table(out1@cluster.sig.data)[aid!=out1@proj,][order(majorCluster),]
out3 =out2[out2$index == "expa",]
out3$majorCluster <- factor(out3$majorCluster,
                            levels = level_cd8 )
library(ggpubr)
my_comparisons = list(
      c("CD4_C3_ANXA1","CD4_C7_OX40")
)
ggboxplot(out3,
          x="majorCluster",y="value",palette = col_cd8$color,
          color = "majorCluster", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = my_comparisons,paired = T,method ="t.test")+
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'none',
  )
ggboxplot(as.data.table(obj@pIndex.sig.migr)[aid!=obj@proj,][order(majorCluster),],
          x="majorCluster",y="value",palette = col_cd8$color,
          color = "majorCluster", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = my_comparisons,
                     paired = T,method ="t.test")+
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'none')  


#Figure 5 f-----
my = subset(sc, L3_C == 'CD4');Idents(my)<-my$L4_C
gene<-c('S1PR5','S1PR1','CCR8','CXCR3','CXCR6')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 6,scale.min = 0)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())

#Figure 5 g-----
my = subset(sc, L3_C == 'CD8');Idents(my)<-my$L4_C
gene<-c('S1PR5','S1PR1','CXCR6')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 6,scale.min = 0)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())

#Figure 5 i-----
my = subset(sc, L4_C == 'Mac_C2_SPP1' & metastasis == 'Y');Idents(my)<-my$tissue;my<-subset(my,ident=c('pLN','Tu'));my$tissue<-Idents(my)
gene<-c('MARCO','CCL2','SPP1','CCL7','MRC1','IL1B')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 8)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())


#Figure 5 j -----
load("ellphonedb.RData")
library(stringr)
library(tidyverse)
library(circlize)
library(reshape2)
com$L = str_split(com$interacting_pair,"_",simplify = T)[,1]
com$R = str_split(com$interacting_pair,"_",simplify = T)[,2]

com$sender = str_split(com$celltype_pairs,"\\.",simplify = T)[,1]
com$recptor = str_split(com$celltype_pairs,"\\.",simplify = T)[,2]

df = com %>% select(c(L,R,sender,mean))

df$L = paste0(df$L,"_",df$sender)

df$mean = as.numeric(df$mean)

rownames(df) <- 1:nrow(df)

fd <- dcast(df,L~R,value.var = "mean")

mat2 <- matrix(rnorm(25), nrow = 5)
rownames(mat2) <- paste0("A", 1:5)
colnames(mat2) <- paste0("C", 1:5)

mat3 <- matrix(rnorm(25), nrow = 5)
rownames(mat3) <- paste0("B", 1:5)
colnames(mat3) <- paste0("C", 1:5)

mat <- matrix(0, nrow = 10, ncol = 10)
rownames(mat) <- c(rownames(mat2), rownames(mat3))
colnames(mat) <- c(colnames(mat1), colnames(mat2))
mat[rownames(mat1), colnames(mat1)] <- mat1
mat[rownames(mat2), colnames(mat2)] <- mat2
mat[rownames(mat3), colnames(mat3)] <- mat3

nm <- unique(unlist(dimnames(mat)))
group <- structure(gsub("\\d", "", nm), names = nm)
grid.col <- structure(
  c(rep("#fb8072", 5), rep("#80b1d3", 5), rep("#fdb462", 5)),
  names = c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5))
)
group <- structure(gsub("\\d", "", nm), names = nm)
chordDiagram(
  mat, group = group, grid.col = grid.col, 
  annotationTrack = c("grid", "axis"),
  preAllocateTracks = list(
    track.height = mm_h(4),
    track.margin = c(mm_h(4), 0)
  )
)

circos.track(
  track.index = 2, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(
      mean(xlim), mean(ylim),
      sector.index, cex = 0.6,
      niceFacing = TRUE
    )
  }, 
  bg.border = NA
)

highlight.sector(
  rownames(mat1), track.index = 1, col = "#fb8072", 
  text = "A", cex = 0.8, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat1), track.index = 1, col = "#80b1d3", 
  text = "B", cex = 0.8, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat2), track.index = 1, col = "#fdb462", 
  text = "C", cex = 0.8, text.col = "white", 
  niceFacing = TRUE
)
circos.clear()



# Extended Data Fig. 11 ------
# Extended Data Fig. 11 a -----
my<-subset(sc,L2_C=='T')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='CD8_C6_PDCD1',]
colnames(sub)[5]<-"Ratio of CD8_C6_PDCD1 in T";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of CD8_C6_PDCD1 in T",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)


my<-subset(sc,L2_C=='T')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='CD4_C7_OX40',]
colnames(sub)[5]<-"Ratio of CD4_C7_OX40 in T";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of CD4_C7_OX40 in T",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

my<-subset(sc,L1_C=='Endothelium')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='Endo_C3_RGCC',]
colnames(sub)[5]<-"Ratio of Endo_C3_RGCC in Endothelium";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of Endo_C3_RGCC in Endothelium",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

my = subset(sc, L3_C == 'Macrophage' | L3_C == 'Monocyte')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='Mac_C2_SPP1',]
colnames(sub)[5]<-"Ratio of Mac_C2_SPP1 in Macrophage+Monocyte";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of Mac_C2_SPP1 in Macrophage+Monocyte",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

my<-subset(sc,L1_C=='Fibroblast')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='FB_C3_COL1A1',]
colnames(sub)[5]<-"Ratio of FB_C3_COL1A1 in Fibroblast";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of FB_C3_COL1A1 in Fibroblast",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)


# Extended Data Fig. 11 b -----
overlap = intersect(VDJ_lpf$clone.id,VDJ_tf$clone.id)
VDJ_lt$share = ifelse(VDJ_lt$clone.id %in% overlap,paste0("share",VDJ_lt$tissue),VDJ_lt$tissue )
tab.1 <- table(VDJ_lt$share,VDJ_lt$patient)
tab.1 <- tab.1[,c(3,4,6,7)]
tab.2 <- melt(tab.1)
colnames(tab.2) <- c("type","patient","number")
tab.2$tt <- tab.2$type
tab.2$tt <- gsub("share","",tab.2$tt)
tab.2 %>% group_by(type,patient)
tab.2<- ddply(tab.2,.(patient),transform,percent=number/sum(number)*100) 
tab.2<- ddply(tab.2,.(patient,tt),transform,percent2=number/sum(number)*100) 
tab.2$label = paste0(sprintf("%.1f", tab.2$percent), "%")
tab.2$label2 = paste0(sprintf("%.1f", tab.2$percent2), "%")
tab.2$sub_type = rep(pro,nrow(tab.2))
tab.2$percentage = tab.2$percent/100
tab.2$percentage2 = tab.2$percent2/100 
tab.2 %>%
    ggplot(aes(fill=type,x=patient,y=percentage,label = label))+
    geom_bar(position = position_stack(),stat = "identity")+
    geom_text(aes(label =label), 
              position = position_stack(vjust = 0.5), size = 6)+
    scale_y_continuous(labels = scales::percent)+
    # geom_text(aes(label=label),vjust=3,size=6,color="black")+
    theme_classic()+
    theme(legend.position = "top",
          legend.text = element_text(size=12),
          axis.text = element_text(size=12),
          axis.title = element_text(size=12,hjust = 0.5),
          plot.title  = element_text(hjust = 0.5,size = 16,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold"),
          axis.text.y = element_text(size = 20,face = "bold"),
          axis.title.x = element_text(size = 20,face = "bold"),
          axis.title.y = element_text(size = 20,face = "bold"),
          legend.title = element_text(size = 20,face = "bold")
    )+labs(title = pro)+scale_fill_d3()
  
  
# Extended Data Fig. 11 c -----
#Endothelium
my<-subset(sc,L1_C=='Endothelium' & metastasis == 'Y')
u<-subset(my,tissue=='Tu');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(u@meta.data$patient,u@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(u@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'Tu'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('pLN','Tu')) 
sub<-com[com$Subcluster=='Endo_C3_RGCC',]
colnames(sub)[5]<-"Ratio of Endo_C3_RGCC in Endothelium";sub$Tissue<-factor(sub$Tissue,levels = c('pLN','Tu'))
ggboxplot(sub,x="Tissue",y="Ratio of Endo_C3_RGCC in Endothelium",color = "Tissue",palette = c("#8DA0CB","#c87137"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

#Fibroblast
my = subset(sc, L1_C=='Fibroblast');my<-subset(my,metastasis == 'Y')
u<-subset(my,tissue=='Tu');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(u@meta.data$patient,u@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(u@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'Tu'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('pLN','Tu')) 
sub<-com[com$Subcluster=='FB_C3_COL1A1',]
colnames(sub)[5]<-"Ratio of FB_C3_COL1A1 in Fibroblast";sub$Tissue<-factor(sub$Tissue,levels = c('pLN','Tu'))
ggboxplot(sub,x="Tissue",y="Ratio of FB_C3_COL1A1 in Fibroblast",color = "Tissue",palette = c("#8DA0CB","#c87137"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

# Extended Data Fig. 11 d -----
my<-subset(sc,L1_C=='Endothelium')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='Endo_C3_RGCC',]
colnames(sub)[5]<-"Ratio of Endo_C3_RGCC in Endothelium";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of Endo_C3_RGCC in Endothelium",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

my<-subset(sc,L1_C=='Fibroblast')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='FB_C3_COL1A1',]
colnames(sub)[5]<-"Ratio of FB_C3_COL1A1 in Fibroblast";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of FB_C3_COL1A1 in Fibroblast",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

# Extended Data Fig. 11 e -----
my = subset(sc, L3_C == 'Cancer' & metastasis == 'Y');Idents(my)<-my$tissue
deg<-FindMarkers(my,assay='RNA',logfc.threshold = 0,ident.1 = 'pLN',ident.2 = 'Tu',only.pos = F,max.cells.per.ident = 1000)
logFC_cutoff <- 0.3
deg$change = as.factor(
  ifelse(deg$p_val < 0.05 & abs(deg$avg_log2FC) > logFC_cutoff,
         ifelse(deg$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
dat<-deg;dat$gene_name<-rownames(dat)
m2d = function(x){
  mean(abs(x))+2*sd(abs(x))
}  
gene<-c('CD24','RAC1','EGFR','S100A9','CXCL17','S100A8','SERPINE1','MDK','STMN1','TOP2A','TUBA1B','TUBB','TUBA1A','HLA-DRB1','HLA-DRA','HLA-A','CD74')
dat<-arrange(dat,desc(abs(avg_log2FC)))
dat<-dat[!(duplicated(dat$gene_name)),]
dat$gene <- as.factor(ifelse(dat$gene_name %in% gene, 'Y', 'N'))
dat$labels <- ''; for (i in gene) {dat[dat$gene_name == i, "labels"] <- i}
ggplot(data = dat, aes(x = avg_log2FC, y = -log10(p_val))) + 
  geom_point(alpha = 0.4, size = 3.5,aes(color = change)) + 
  geom_text_repel(aes(label=labels),size=5,max.overlaps = 10000,force=2,show.legend=F)+
  scale_color_manual(values = c("blue", "gray", "red")) + 
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 4, col = "black", lwd = 0.8) + 
  geom_hline(yintercept = -log10(0.05), lty = 4,col = "black", lwd = 0.8) + 
  theme_bw() + 
  labs(x = "log2FC", y = "-log10(P)") + 
  theme(plot.title = element_text(hjust = 0.5),panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),axis.text= element_text(colour = "black",size=14),axis.title = element_text(colour = "black",size=14),
        legend.text = element_text(colour = "black",size=14),legend.title = element_text(colour = "black",size=14))


# Extended Data Fig. 11 f -----
library(GseaVis)
library(clusterProfiler)
my = subset(sc, L3_C == 'Cancer' & metastasis == 'Y');Idents(my)<-my$tissue
deg<-FindMarkers(my,assay='RNA',ident.1 = 'pLN',ident.2='Tu',only.pos = T);deg$Gene<-rownames(deg)
rank<-deg[,c('Gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
human_H<-msigdbr(species = "human", category = "H")
human_C2<-msigdbr(species = "human", category = "C2")
human_C2<-human_C2[human_C2$gs_subcat%in%c('CP:BIOCARTA','CP:KEGG','CP:REACTOME'),]
human<-rbind(human_C2,human_H);gmt<-human %>% dplyr::select(gs_name,gene_symbol)
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000,
                pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('REACTOME_DEVELOPMENTAL_BIOLOGY','REACTOME_ADAPTIVE_IMMUNE_SYSTEM')
gseaNb(object = gseaRes,geneSetID = geneSetID,addPoint = F, addPval = T,
       pCol = 'black', pHjust = 0,subPlot = 2,pvalX = 1,pvalY = 1.1,curveCol = brewer.pal(n = 2, name = "Set1"))

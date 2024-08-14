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
dotplot <- function(mca,features,celltype,ctidx,color){
  exp <- mca@assays$RNA@data[features,] %>% as.data.frame() %>% t() %>% as.data.frame()
  exp$subC <- mca@meta.data[rownames(exp),'L1_C']
  
  trunc_z <- function(x) MinMax((x-mean(x))/sd(x),-2,2)
  
  exp %>% reshape2::melt() %>% group_by(subC,variable) %>% #summarise(N = n()) %>% as.data.frame()
    summarise(mean = mean(value),frac = sum(value > 0)/length(value)) %>% 
    group_by(variable) %>% mutate(mean_z = trunc_z(mean)) -> plot_df
  
  gType <- rep(celltype,ctidx)
  names(gType) <- features
  
  plot_df$gType <- gType[plot_df$variable]
  plot_df$gType <- factor(plot_df$gType,
                          levels = celltype)
  
  ggplot(plot_df,aes(x = variable,y = subC, fill = mean_z,size = frac)) +
    geom_point(shape=21,color='black') +
    scale_fill_gradientn(colors = color) +
    scale_y_discrete(limits=rev(levels(plot_df$subC))) +
    scale_size_continuous(range = c(1,8),name = 'Proportion') +
    facet_grid(~gType,space = 'free_x',scale = 'free_x') +
    theme_bw() + xlab('') + ylab('') +
    theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank())
}


# RO/E
library(pheatmap)
library(Seurat)
library(plyr)
library(ggpubr)
library(ggpmisc)

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}




######Figure 1 #####
#Fig1 b -----
#major
load("scRNA.RData")
FetchData(sc,c("UMAP_1",'UMAP_2','L1_C'))%>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = L1_C),size =.4/.pt) +
  scale_color_manual(values =pal_igv("default")(20)) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank()
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))
# endo -----
load("scRNA_endo.RData")
FetchData(sc,c("UMAP_1",'UMAP_2','subC')) %>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#DEBA6A','#E63228','#7D5599','#3787BC')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(t = 10,  # 顶部边缘距离
        # r = 40,  # 右边边缘距离
        # b = 40,  # 底部边缘距离
        # l = 10)
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))
# FB ----
load("scRNA_FB.RData")
FetchData(my,c("UMAP_1",'UMAP_2','subC'))  %>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#D4C799','#4BAC3E','#BA6C35','#D5A65C','#fdbf73','#800000','#95D075')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(t = 10,  # 顶部边缘距离
        # r = 40,  # 右边边缘距离
        # b = 40,  # 底部边缘距离
        # l = 10)
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))
# B -----
load("scRNA_B.RData")
FetchData(my,c("UMAP_1",'UMAP_2','subC'))  %>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#ED5152','#F06C45','#6BAD9E','#579BC7','#FDB258','#AC8DC3','#E5DD99','#8E6C99')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(t = 10,  # 顶部边缘距离
        # r = 40,  # 右边边缘距离
        # b = 40,  # 底部边缘距离
        # l = 10)
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))
# CD4T -----
load("scRNA_CD4T.RData")
FetchData(my,c("UMAP_1",'UMAP_2','subC'))  %>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#ED5152','#F06C45','#6BAD9E','#579BC7','#FDB258','#AC8DC3','#E5DD99','#8E6C99')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(t = 10,  # 顶部边缘距离
        # r = 40,  # 右边边缘距离
        # b = 40,  # 底部边缘距离
        # l = 10)
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))

# Myeloid ------
load("scRNA_myeloid.RData")
FetchData(my,c("UMAP_1",'UMAP_2','subC')) +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#EA4243','#C79B7D','#3989b8','#752773FF','#DAA520',
                               '#F88D8D','#c9c193','#9FD18F','#F06061','#48D1CC',
                               '#3d804e','#8A64AE','#b7a3c4','#FE9424')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(t = 10,  # 顶部边缘距离
        # r = 40,  # 右边边缘距离
        # b = 40,  # 底部边缘距离
        # l = 10)
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))
# CD4T -----
load("scRNA_CD4T.RData")
FetchData(my,c("UMAP_1",'UMAP_2','subC'))  %>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#6B3F99','#40E0D0','#999C63','#B19999','#8B4513',
                               '#F0E084','#267DB1','#F7861D','#E73335','#A4D880',
                               '#20B2AA','#CD853F','#69BA53')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(t = 10,  # 顶部边缘距离
        # r = 40,  # 右边边缘距离
        # b = 40,  # 底部边缘距离,header = T,row.names = 1),header = T,row.names 
        # l = 10)
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))
# CD8T -----
load("scRNA_CD8T.RData")
FetchData(my,c("UMAP_1",'UMAP_2','subC'))  %>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#4169E1','#9F8399','#66A5CD','#77C25F','#A6CEE3',
                               '#F37F4F','#B79BCA','#5AA1A3','#C2B099','#CC934F',
                               '#E31F1E','#7CB999')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(t = 10,  # 顶部边缘距离
        # r = 40,  # 右边边缘距离
        # b = 40,  # 底部边缘距离
        # l = 10)
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))

#Fig s1 a -----
load("/workspace/wangxiliang/esca/sc_rna/scRNA_t.RData")
my$L4_C<-Idents(my)
my$L4_C<-factor(my$L4_C,levels = c('T_C2_MKI67','T_C1_HSPA1A','CD8','CD4','NK_C2_XCL1','NK_C1_CD16'))
T_plotG <- c('SPON2','FCGR3A','CX3CR1','CD160','AREG','XCL1','FCER1G','XCL2',
             'MAL','CD4','ICOS','IL6ST','CD8A','CCL4L2','CD8B','GZMK','HSPA1B','HSPA1A','DNAJB1','JUN','STMN1',
             'MKI67','RRM2','TOP2A')
celltype <- c('NK_C1_CD16',  'NK_C2_XCL1','CD4','CD8', 'T_C1_HSPA1A',  'T_C2_MKI67')
ctidx <- c(4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Fig s1 b -----
FetchData(my,c("UMAP_1",'UMAP_2','subC'))  %>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#EC593B','#F7F499','#1E90FF','#20B2AA','#DF9A89','#B0DD8A')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(t = 10,  # 顶部边缘距离
        # r = 40,  # 右边边缘距离
        # b = 40,  # 底部边缘距离
        # l = 10)
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))
#Fig s1 c -----
T_plotG <- c('SPON2','FCGR3A','CX3CR1','CD160','AREG','XCL1','FCER1G','XCL2',
             'MAL','CD4','ICOS','IL6ST','CD8A','CCL4L2','CD8B','GZMK','HSPA1B','HSPA1A','DNAJB1','JUN','STMN1',
             'MKI67','RRM2','TOP2A')
plot_df <- read.csv("/home/users/zhangzhichao/workspace/project/ESCC/NKT_marker_gene_mean_z.csv",header = T)
plot_df$subC = factor(plot_df$subC,
                       levels = c("NK_C1_CD16","NK_C2_XCL1","CD4","CD8","T_C2_MKI67","T_C1_HSPA1A"))
plot_df$variable = factor(plot_df$variable ,
                          levels = T_plotG)

celltype <- c('NK_C1_CD16',  'NK_C2_XCL1','CD4','CD8', 'T_C1_HSPA1A',  'T_C2_MKI67')
ctidx <- c(4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(plot_df,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Fig s1 d -----
Idents(sc)<-sc$tissue
t<-as.matrix(table(Idents(sc),sc@meta.data$L1_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(sc)))
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
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values=pal_igv("default")(20))

Idents(sc)<-sc$patient
t<-as.matrix(table(Idents(sc),sc@meta.data$L1_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(sc)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =pal_igv("default")(20))


##Extended Data Fig. 2-----
#Extended Data Fig. 2 a 
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD4'));Idents(my)<-my$L4_C
T_plotG <- c('CD55','CCR7','ABLIM1','LEF1','KLF2','TOX2','PDCD1','TOX','IL6ST','CXCR5','IL22','IL17A','IL17F','CCR6','RORA',
             'CCL5','GZMK','NKG7','PRF1','CTSW','ISG15','IFI6','IFI44L','IFIT1','IRF7','TCF7','TXNIP','MYC','NOSIP','TLE5','ANXA1','FTH1','MYADM','IL7R','RGCC',
             'CD69','PTGER4','SLC2A3','BTG2','GPR183','SESN3','FCRL3','RTKN2','IKZF2','TGIF1','HLA-DRB1','FOXP3','HLA-DRA',
             'GBP5','IL2RA','TNFRSF4','TNFRSF18','LAIR2','CTLA4','TNFRSF9','ITGB1','AQP3','TRADD','CDC25B','SH3BP5','CXCL13',
             'IFNG','IL21','BHLHE40','GADD45G')
celltype <- c('CD4_C1_CCR7', 'CD4_C10_CXCR5', 'CD4_C11_IL17A','CD4_C12_NKG7', 'CD4_C13_ISG15','CD4_C2_TCF7',  
              'CD4_C3_ANXA1', 'CD4_C4_GPR183','CD4_C5_RTKN2','CD4_C6_CD25','CD4_C7_OX40','CD4_C8_ITGB1','CD4_C9_IFNG')
ctidx <- c(5,5,5,5,5,5,5,5,5,5,5,5,5)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Extended Data Fig. 2 b 
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD8'));Idents(my)<-my$L4_C
T_plotG <- c('CCR7','LEF1','SELL','RPS13','RPL32','CCL4','FCGR3A','CCL3','SPON2','BHLHE40','ISG15','IFI44L','IFIT1','IFIT3','OAS1','KLRB1','NCR3','ZBTB16',
             'SLC4A10','TRAV1-2','IL7R','RGCC','SLC2A3','ZNF331','PTGER4','LTB','TCF7','TXNIP','LDLRAP1','FLT3LG','LMNA','ANXA1','CRIP1','CAPG','PLP2',
             'GZMK','TNFSF9','SH2D1A','CRTAM','TRAT1','CXCL13','ENTPD1','TIGIT','HAVCR2','CTLA4','CX3CR1','ADGRG1','PLEK','S1PR5',
             'FCRL6','KIR3DL2','KIR2DL4','ZNF683','KIR2DL3','XCL1','TYROBP','KLRC2','KLRF1','KLRC3','CD160')
celltype <- c('CD8_C1_CCR7',  'CD8_C10_CD16', 'CD8_C11_ISG15',  'CD8_C12_CD161','CD8_C2_IL7R',   'CD8_C3_TCF7',  'CD8_C4_ANXA1',   'CD8_C5_GZMK',
              'CD8_C6_CD39','CD8_C7_CX3CR1','CD8_C8_KIR',  'CD8_C9_KLRC2')
ctidx <- c(5,5,5,5,5,5,5,5,5,5,5,5)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Extended Data Fig. 2 c
Idents(sc)<-sc$L1_C;my<-subset(sc,ident=c('B','Plasma'));Idents(my)<-my$L4_C
T_plotG <- c('NME1','EIF4A1','CCR7','EIF5A','CRIP1','DUSP4','ITGB1','ANXA2','HSPA1A','GPR183','FOSB','KLF6','TCL1A','FCER2',
             'IGHD','FCRL1','IFITM1','ISG15','IFI6','IFIT2','BCL6','BCL7A','CD38','SERPINA9','STMN1','HMGB2','UBE2C','MKI67',
             'JCHAIN','IGHG1','MZB1','SSR4')
celltype <- c('B_C1_CCR7','B_C2_DUSP4','B_C3_GPR183','B_C4_TCL1A','B_C5_ISG15', 'B_C6_BCL6','B_C7_MKI67','Plasma')
ctidx <- c(4,4,4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Extended Data Fig. 2 d
my<-subset(sc,L1_C=='Endothelium');Idents(my)<-my$L4_C
T_plotG <- c('ACKR1','SELP','SELE','CCL23','SEMA3G','FBLN5','GJA5','SERPINE2','PLVAP','RGCC','APLNR','IGFBP5',
             'CCL21','PROX1','PDPN','FLT4')
celltype <- c('Endo_C1_ACKR1','Endo_C2_FBLN5','Endo_C3_RGCC','Endo_C4_CCL21')
ctidx <- c(4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Extended Data Fig. 2 e
Idents(sc)<- sc$L2_C
my<-subset(sc,ident=c('Myeloid'));Idents(my)<-my$L4_C
T_plotG <- c('CLEC9A','XCR1','CPNE3','SNX3','CD1C','FCER1A','CLEC10A','HLA-DQA1','CD1A','CD1E','LTB','CD207','LAMP3','CCR7',
             'CCL19','EBI3','IL3RA','GZMB','JCHAIN','LILRA4','NLRP3','EREG','IL1B','AREG','SPP1','MMP12','MMP9','CSTB','TREM2',
             'APOC1','C1QC','LIPA','LYVE1','CCL18','FOLR2','CCL2','CXCL10','CXCL9','ISG15','GBP1','MKI67','STMN1','TUBA1B','HMGN2',
             'S100A9','S100A8','FCN1','VCAN','RHOC','LST1','MS4A7','LILRB2','CXCL8','G0S2','FCGR3B','CXCR2')
celltype <- c('DC_C1_CLEC9A','DC_C2_CD1C','DC_C3_CD1A', 'DC_C4_LAMP3','DC_C5_IL3RA','Mac_C1_NLRP3', 'Mac_C2_SPP1','Mac_C3_C1QC',
              'Mac_C4_LYVE1', 'Mac_C5_CXCL10','Mac_C6_MKI67', 'Mono_C1_CD14', 'Mono_C2_CD16', 'Neutrophil')
ctidx <- c(4,4,4,4,4,4,4,4,4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))
#Extended Data Fig. 2 f
my<-subset(sc,L1_C=='Fibroblast');Idents(my)<-my$L4_C;my<-subset(my,ident='Smooth muscle cell',invert=T)
T_plotG <- c("CFD","IGFBP6","GPX3","PI16","IGF1","DPT","CXCL12","COL4A4","MMP11","POSTN","COL1A1","CTHRC1","APOD","COL13A1",
             "CXCL14","CXCL1","HIGD1B","COX4I2","KCNJ8","PDGFRB","MUSTN1","ADIRF","ACTA2","ADAMTS9")
celltype <- c('FB_C1_CFD','FB_C2_IGF1','FB_C3_COL1A1','FB_C4_APOE','FB_C5_PDGFRB','FB_C6_ACTA2')
ctidx <- c(4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))




#Extended Data Fig. 3 -----
#B/Plasma
load("scRNA_b.RData")
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

#myeloid
Idents(sc)<- sc$L2_C
my<-subset(sc,ident=c('Myeloid'));Idents(my)<-my$L4_C
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

#Endothelium
my<-subset(sc,L1_C=='Endothelium');Idents(my)<-my$L4_C
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


#Fibroblast
my<-subset(sc,L1_C=='Fibroblast');Idents(my)<-my$L4_C;my<-subset(my,ident='Smooth muscle cell',invert=T)
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


#CD4T
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD4'));Idents(my)<-my$L4_C
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


#CD8T 
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD8'));Idents(my)<-my$L4_C
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


#NK/T
my<-subset(sc,L1_C=='NK/T');Idents(my)<-my$L3_C
t<-as.matrix(table(my$tissue,Idents(my)))
t<-as.data.frame(t)
total<-as.data.frame(table(my$tissue))
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
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values =c('#EC593B','#F7F499','#1E90FF','#20B2AA','#DF9A89','#B0DD8A'))

t<-as.matrix(table(my$patient,Idents(my)))
t<-as.data.frame(t)
total<-as.data.frame(table(my$patient))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Patient<-factor(t$Patient,levels = c('826','128','308','1104','1119','1204','1209','1229','118','1201','1019','1022','1123','1210'))
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =c('#EC593B','#F7F499','#1E90FF','#20B2AA','#DF9A89','#B0DD8A'))


#Epithelium/Cancer
my<-subset(sc,L1_C=='Epithelium/Cancer');Idents(my)<-my$L3_C
t<-as.matrix(table(my$tissue,Idents(my)))
t<-as.data.frame(t)
total<-as.data.frame(table(my$tissue))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t<-t[t$Tissue!='PBMC',]
t$Tissue<-factor(t$Tissue,levels = c('nLN','pLN','Nor','Adj','Tu'))
ggplot(t, aes(Tissue,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values =c('#1a759f','#a44a3f'))

t<-as.matrix(table(my$patient,Idents(my)))
t<-as.data.frame(t)
total<-as.data.frame(table(my$patient))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Patient<-factor(t$Patient,levels = c('826','128','308','1104','1119','1204','1209','1229','118','1201','1019','1022','1123','1210'))
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =c('#1a759f','#a44a3f'))





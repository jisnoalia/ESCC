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

#Fig 1d ISCHIA -----
cut = 15
sp <- Composition.cluster(sp, norm_weights, cut)
sp$cc_15 <- sp$CompositionCluster_CC
sp.umap <- Composition_cluster_umap(sp, norm_weights)
emb.umap <- sp.umap$umap.table
emb.umap$CompositionCluster_CC <- NULL
emb.umap$Slide <- NULL
emb.umap <- as.matrix(emb.umap)
colnames(emb.umap) <- c("UMAP1", "UMAP2")
sp[['umap.ischia15']] <- CreateDimReducObject(embeddings = emb.umap, key = 'umap.ischia15_', assay = 'rctd_full')
image_names <- unique(sp$file)
paletteMartin <- color$color
all_ccs <- unique(sp$cc_15)
color_mapping15 <- setNames(dittoColors()[1:cut][1:length(all_ccs)], sort(all_ccs))
df<-read.table("sp_table.txt", sep='\t',header = T)
sp$cellID = rownames(sp@meta.data)
mye <- sp[['umap.ischia15']]@cell.embeddings %>% as.data.frame()
mye$cellID = rownames(mye)
phe = sp@meta.data %>% select(cellID,cc_15,file)
mye = mye %>% left_join(phe,by = "cellID")
mye %>%
  ggplot() +
  ggrastr::geom_point_rast(aes(umapischia15_1,umapischia15_2,color = cc_15  ),size =.4/.pt) +
  scale_color_manual(values = color_mapping15) +
  theme_classic() +
  # NoLegend()+
  coord_fixed(ratio = 1) +
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank()
  )+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) +
  guides(colour = guide_legend(ncol = 2,
                               override.aes=list(shape=19, size=3, linetype=0)))
}
#Fig 1e -----
development = subset(sp,development=='unknown',invert=T)

Idents(development)<-development$development
t<-as.matrix(table(Idents(development),development@meta.data$cc_15))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(development)))
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
  ylab('Proportion')+scale_fill_manual(values =c('#ED5152','#F06C45','#6BAD9E','#579BC7','#FDB258','#AC8DC3','#E5DD99','#8E6C99'))


ME = subset(sp,ME=='unknown',invert=T)

Idents(ME)<-ME$ME
t<-as.matrix(table(Idents(ME),ME@meta.data$cc_15))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(ME)))
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
  ylab('Proportion')+scale_fill_manual(values =c('#ED5152','#F06C45','#6BAD9E','#579BC7','#FDB258','#AC8DC3','#E5DD99','#8E6C99'))

#Fig 1f -----
common_xlim <- c(50, 550)
common_xbreaks <- seq(50, 550, by = 100)
pdf("spatial_plots_K15.pdf", width = 7, height = 6)
for (image_name in image_names) {
  plot <- SpatialDimPlot(sp, group.by = "cc_15", images = df[df$file==image_name,"image"],stroke = NA,image.alpha = 0,pt.size.factor = 1.5) +
    scale_fill_manual(values =  color_mapping15) +
    theme_minimal() +
    ggtitle(image_name)+
    guides(fill = guide_legend(override.aes = list(size=4)))+
    coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
    scale_x_continuous(breaks = common_xbreaks) +
    scale_y_continuous(breaks = common_xbreaks) +
    theme(aspect.ratio = 1)
  print(plot)
}
dev.off()

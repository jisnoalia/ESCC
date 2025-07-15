#Figure 7 -----
#Figure a -----
load("atac_cd8.Rdata")
Idents(atac)<-atac$L4_C;Idents(atac)<-gsub('MAIT','CD8_C12_CD161',Idents(atac))
mylevel<-c('CD8_C1_CCR7',  'CD8_C10_CD16', 'CD8_C11_ISG15',  'CD8_C12_CD161','CD8_C2_IL7R',   'CD8_C3_TCF7',  'CD8_C4_ANXA1',   'CD8_C5_GZMK',
  'CD8_C6_CD39','CD8_C7_CX3CR1','CD8_C8_KIR',  'CD8_C9_KLRC2')
Idents(atac)<-factor(Idents(atac),levels = mylevel)
DimPlot(atac, pt.size=0.1, reduction="umap", label=F, cols = c('#4169E1','#9F8399','#66A5CD','#7CB999','#77C25F','#A6CEE3',
#Figure b -----
library(Seurat)
library(ggrepel)
library(Signac)
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/liuliqiu/ATAC/CD8_ATAC.RData")

DefaultAssay(atac3) <- 'ATAC'
da_peaks <- FindMarkers(
  object = atac3,
  ident.1 = "CD8_C6_CD39",
  logfc.threshold = 0,
  latent.vars = 'peak_region_fragments',
  only.pos = FALSE,
  max.cells.per.ident = 1000
)

open_cd8 <- rownames(da_peaks)
closest_genes_cd8 <- ClosestFeature(
  object = atac3,
  regions = open_cd8
)
da_peaks$query_region<-rownames(da_peaks);da_peaks<-merge(da_peaks,closest_genes_cd8[,c(2,7)]);
da_peaks<-arrange(da_peaks,desc(avg_log2FC))
mean(da_peaks$avg_log2FC) +sd(da_peaks$avg_log2FC)

logFC_cutoff <- 0.3
da_peaks$change = as.factor(
  ifelse(da_peaks$p_val < 0.05 & abs(da_peaks$avg_log2FC) > logFC_cutoff,
         ifelse(da_peaks$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
dat<-da_peaks
m2d = function(x){
  mean(abs(x))+2*sd(abs(x))
}  
gene<-c('GEM','LAYN','CXCL13','DUSP4','PDCD1','ENTPD1','TOX','S1PR5','LYAR','SORL1')
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
  labs(x = "log2FC", y = "-log10(P.value)") + 
  theme(plot.title = element_text(hjust = 0.5),panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),axis.text= element_text(colour = "black",size=14),
        axis.title = element_text(colour = "black",size=14),
        legend.text = element_text(colour = "black",size=14),legend.title = element_text(colour = "black",size=14))+xlim(-0.9,0.96)

dat_atac = dat

## RNA
load("/workspace/wangxiliang/project/escc_baseline_zhanQiMin/sc_rna/scRNA.RData")
#CD8 volcano
my<-subset(sc,L3_C=='CD8');Idents(my)<-my$L4_C
deg<-FindMarkers(my,assay='RNA',logfc.threshold = 0,ident.1 = 'CD8_C6_CD39',only.pos = F,max.cells.per.ident = 1000)
logFC_cutoff <- 0.3
deg$change = as.factor(
  ifelse(deg$p_val < 0.05 & abs(deg$avg_log2FC) > logFC_cutoff,
         ifelse(deg$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
dat<-deg;dat$gene_name<-rownames(dat)
m2d = function(x){
  mean(abs(x))+2*sd(abs(x))
}  
gene<-c('GEM','LAYN','CXCL13','DUSP4','PDCD1','ENTPD1','TOX','S1PR5','LYAR','SORL1')
dat<-arrange(dat,desc(abs(avg_log2FC)))
dat<-dat[!(duplicated(dat$gene_name)),]
dat$gene <- as.factor(ifelse(dat$gene_name %in% gene, 'Y', 'N'))
dat$labels <- ''; for (i in gene) {dat[dat$gene_name == i, "labels"] <- i}
ggplot(data = dat_rna, aes(x = avg_log2FC, y = -log10(p_val))) + 
  geom_point(alpha = 0.4, size = 3.5,aes(color = change)) + 
  geom_text_repel(aes(label=labels),size=5,max.overlaps = 10000,force=2,show.legend=F)+
  scale_color_manual(values = c("blue", "gray", "red")) + 
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 4, col = "black", lwd = 0.8) + 
  geom_hline(yintercept = -log10(0.05), lty = 4,col = "black", lwd = 0.8) + 
  theme_bw() + 
  labs(x = "log2FC", y = "-log10(P.value)") + 
  theme(plot.title = element_text(hjust = 0.5),panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),axis.text= element_text(colour = "black",size=14),axis.title = element_text(colour = "black",size=14),
        legend.text = element_text(colour = "black",size=14),legend.title = element_text(colour = "black",size=14))+xlim(-3,3)

dat_rna = dat
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/liuliqiu/RCTD/dat_all.Rdata")
all = intersect(dat_rna$gene_name,dat_atac$gene_name)

dat_rna_use = dat_rna %>% filter(gene_name %in% all)
dat_atac_use = dat_atac %>% filter(gene_name %in% all)
logFC_cutoff= 0.3
dat_atac_use$change = as.factor(
  ifelse(dat_atac_use$p_val < 0.05 & abs(dat_atac_use$avg_log2FC) > logFC_cutoff,
         ifelse(dat_atac_use$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
head(dat_rna_use)
head(dat_atac_use)
table(dat_atac_use$change)
dat_rna_use %>% left_join(dat_atac_use,by = "gene_name") -> dat_all
colnames(dat_all) = gsub("\\.x","_rna",colnames(dat_all))
colnames(dat_all) = gsub("\\.y","_atac",colnames(dat_all))
head(dat_all)
dat_all$change = paste0(dat_all$change_rna,"_",dat_all$change_atac)
table(dat_all$change )
dat_all %>%
  mutate(change = recode(change, 
                         "DOWN_DOWN" = "Both_DOWN", 
                         "DOWN_NOT" = "Only_DOWN_in_RNA",
                         "NOT_DOWN" = "Only_DOWN_in_ATAC",
                         "NOT_NOT" = "Both_NOT",
                         "NOT_UP" = "Only_UP_in_ATAC",
                         "UP_NOT" = "Only_UP_in_RNA",
                         "UP_UP" = "Both_UP")) ->dat_all_df
top5_rna = dat_all_df %>% arrange(desc(avg_log2FC_rna)) %>% head(5)
top5_rna
top5_atac = dat_all_df %>% arrange(desc(avg_log2FC_atac)) %>% head(5)
top5_atac
inter_5 = intersect(top5_rna$gene_name,top5_atac$gene_name)
head(dat_all_df)
dat_all_df$top5_log2FC = ''
dat_all_df$top5_log2FC = ifelse(dat_all_df$gene_name %in% unique(c(top5_rna$gene_name,top5_atac$gene_name)),"top5 log2FC",dat_all_df$top5_log2FC)
dat_all_df$top5_log2FC = ifelse(dat_all_df$gene_name==inter_5,"both in top5 log2FC",dat_all_df$top5_log2FC)
dat_all_df$top5_log2FC  = factor(dat_all_df$top5_log2FC ,
                                 levels = c("","top5 log2FC","both in top5 log2FC"))
dat_all_df$labels_atac = ifelse(dat_all_df$gene_name %in% "GALNT2","GALNT2",dat_all_df$labels_atac)
dat_all_df$labels_atac = ifelse(dat_all_df$gene_name %in% "MAST4","MAST4",dat_all_df$labels_atac)
dat_all_df$labels_atac = ifelse(dat_all_df$gene_name %in% "DAPK2","DAPK2",dat_all_df$labels_atac)
dat_all_df$labels_atac = ifelse(dat_all_df$gene_name %in% "INPP5F","INPP5F",dat_all_df$labels_atac)
dat_all_df$labels_atac = ifelse(dat_all_df$gene_name %in% "ITGAE","ITGAE",dat_all_df$labels_atac)
dat_all_df %>% 
  ggplot(.,aes(y= avg_log2FC_rna , x = avg_log2FC_atac ))+
  geom_point_rast(aes(color=change,size=top5_log2FC),alpha = 0.7)+
  geom_text_repel(aes(label=labels_atac,color=change,segment.color=change),
                  size=5,max.overlaps = 10000,force=2,
                  arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                  segment.color="black",
                  segment.size=0.5,segment.alpha=0.8,
                  nudge_y=0.4,
                  show.legend=F)+
  # geom_smooth(method = "lm", se = FALSE, color = "blue") +
  scale_color_manual(values = c("Both_DOWN"="#007abb", 
                                "Only_DOWN_in_RNA" = "#2ac2d1",
                                "Only_DOWN_in_ATAC" = "#58e07c",
                                "Both_NOT" = "grey",
                                "Only_UP_in_ATAC" = "#ffc149" ,
                                "Only_UP_in_RNA" = "#ff6d36" ,
                                "Both_UP" = "#fc0302"))+ # #fc0302  #f02656
  geom_vline(xintercept = c(-0.3, 0.3), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = c(-0.3, 0.3), lty = 4, col = "black", lwd = 0.8) +
  scale_size_manual(values = c(2,4, 8)) + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black",size=2),
        panel.grid.major = element_blank(),
        axis.text= element_text(colour = "black",size=14),
        axis.title = element_text(colour = "black",size=14),
        legend.text = element_text(colour = "black",size=14),
        legend.title = element_text(colour = "black",size=14),
        # legend.position = "top"
  )+
  guides(colour = guide_legend(ncol =1,
                               override.aes=list(shape=19, size=4, linetype=0)))+
  labs(x="Log2FC in scATAC",y="Log2FC in scRNA")
#Figure 3 c ------
load("/workspace/wangxiliang/esca/sc_rna/scRNA.RData")
#CD8 volcano
my<-subset(sc,L3_C=='CD8');Idents(my)<-my$L4_C
T_plotG <- c('TOX','PDCD1','LAYN','GEM','ENTPD1','DUSP4','CXCL13')
celltype <- c('CD8_C1_CCR7',  'CD8_C10_CD16', 'CD8_C11_ISG15',  'CD8_C12_CD161','CD8_C2_IL7R',   'CD8_C3_TCF7',  'CD8_C4_ANXA1',   'CD8_C5_GZMK',
              'CD8_C6_CD39','CD8_C7_CX3CR1','CD8_C8_KIR',  'CD8_C9_KLRC2')
ctidx <- c(1,1,1,1,1,1,1,1,1,1,1,1)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))
#Figure 3 d------
CoveragePlot(
  object = atac,
  region = 'GEM',features = 'GEM',
  expression.assay = "SCT",
  extend.upstream = 0,
  extend.downstream = 1000,
  ncol = 1
)&scale_fill_manual(values=mycolor)& 
  theme(axis.text= element_text(colour = "black",size=10))

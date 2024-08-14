#Figure 3 a------
load("atac_cd8.Rdata")
Idents(atac)<-atac$L4_C;Idents(atac)<-gsub('MAIT','CD8_C12_CD161',Idents(atac))
mylevel<-c('CD8_C1_CCR7',  'CD8_C10_CD16', 'CD8_C11_ISG15',  'CD8_C12_CD161','CD8_C2_IL7R',   'CD8_C3_TCF7',  'CD8_C4_ANXA1',   'CD8_C5_GZMK',
  'CD8_C6_CD39','CD8_C7_CX3CR1','CD8_C8_KIR',  'CD8_C9_KLRC2')
Idents(atac)<-factor(Idents(atac),levels = mylevel)
DimPlot(atac, pt.size=0.1, reduction="umap", label=F, cols = c('#4169E1','#9F8399','#66A5CD','#7CB999','#77C25F','#A6CEE3',
                                                               '#F37F4F','#B79BCA','#5AA1A3','#C2B099','#CC934F','#E31F1E'),raster = F)

#Figure 3 b------
DefaultAssay(atac) <- 'ATAC'
da_peaks <- FindMarkers(
  object = atac,
  ident.1 = "CD8_C6_CD39",
  logfc.threshold = 0,
  latent.vars = 'peak_region_fragments',
  only.pos = FALSE,
  max.cells.per.ident = 1000
)

open_cd8 <- rownames(da_peaks)
closest_genes_cd8 <- ClosestFeature(
  object = atac,
  regions = open_cd8
)
da_peaks$query_region<-rownames(da_peaks);da_peaks<-merge(da_peaks,closest_genes_cd8[,c(2,7)]);da_peaks<-arrange(da_peaks,desc(avg_log2FC))
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
        panel.grid.major = element_blank(),axis.text= element_text(colour = "black",size=14),axis.title = element_text(colour = "black",size=14),
        legend.text = element_text(colour = "black",size=14),legend.title = element_text(colour = "black",size=14))+xlim(-3,3)

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

#Extended Data Fig. 7 -----

#Extended Data Fig. 7 a 
co=subset(sc,subset =L4_C=='CD8_C6_CD39')
FeatureScatter(co, pt.size = 2, feature1 = 'GEM', feature2 = 'ENTPD1') +stat_cor(method = "pearson",size=6,vjust=0.5)+
  stat_regline_equation(aes(label=paste(..eq.label..)),size=6,vjust=1.6)+stat_smooth(method='lm', color="black", se=TRUE)+labs(title='')+guides(size=FALSE)+theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  scale_color_manual(values ='#F37F4F')
#Extended Data Fig. 7 b 
load('ESCC_RNA_FPKM_edit.RData')
head(meta)
library(ggrepel)
library(ggrastr)

exprSet %>% rownames_to_column("gene") %>% 
  filter(gene %in%c("GEM",'PDCD1','CD8A'))  %>%
  column_to_rownames("gene") %>% t()  %>%
  as.data.frame()->tmp

tmp$gem  = tmp$GEM/tmp$CD8A
tmp$pdcd1 = tmp$PDCD1/tmp$CD8A

ggplot(tmp, aes(x=log2(gem), y=log2(pdcd1))) +
  geom_point(color="#F37F4F",size=4)+ 
  stat_cor(method = "pearson",size=6,vjust=0.5,digits = 2)+
  stat_smooth(method='lm', 
              color="black", se=TRUE)+labs(title='')+guides(size=FALSE)+
  theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  theme_classic()+
  theme(axis.text = element_text(size=12,colour = "black"),
        axis.title = element_text(size=15))+
  labs(x="log2(GEM/CD8A FPKM)", y="log2(PDCD1/CD8A FPKM")

#Extended Data Fig. 7 c 
DefaultAssay(sp)<-'spatial'
p<-subset(sp,CD8A>0&ENTPD1>0&GEM>0&PDCD1>0)
DefaultAssay(p)<-'spatial';p$group<-'CD8+CD39+';Idents(p)<-p$group
FeatureScatter(p, pt.size = 2, feature1 = 'GEM', feature2 = 'PDCD1') +stat_cor(method = "pearson",size=6,vjust=0.5)+
  stat_regline_equation(aes(label=paste(..eq.label..)),size=6,vjust=1.6)+stat_smooth(method='lm', color="black", se=TRUE)+labs(title='')+guides(size=FALSE)+theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  scale_color_manual(values ='#F37F4F')




#Extended Data Fig. 8 b c ------------------ 
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsci)
library(patchwork)
bulk <- read.table("../GSE194116/GSE194116_Raw_gene_counts_matrix.txt",header = T)
bulk <- bulk[!duplicated(bulk$gene_name),]
rownames(bulk) <- bulk$gene_name
exprSet <- bulk[,c(3:14)]
colnames(exprSet) <- gsub(".count","",colnames(exprSet))
keep_feature <- rowSums (exprSet > 1) > 1 
ensembl_matrix <- exprSet[keep_feature, ] #筛选表达的基因
dat=log2(edgeR::cpm(ensembl_matrix)+1)
use_gene <- c("ANTXR1")
ensembl_matrix_use <- ensembl_matrix[rownames(ensembl_matrix) %in% use_gene,]
dat_use <- dat[rownames(dat) %in% use_gene,]
cpm_mat <- melt(dat_use)
cpm_mat$group <- gsub("\\.\\d*","",cpm_mat$Var2)
cpm_mat$sample <- gsub("\\w*\\.","",cpm_mat$sample)
cpm_mat$group <- gsub("ECA","Tumor",cpm_mat$group)
cpm_mat$group <- factor(cpm_mat$group,
                       levels = c("Tumor","Normal"))
ggplot(cpm_mat,aes(x=group,y=log2cpm,fill=group))+
    geom_boxplot()+
    geom_line(aes(group = sample), color = 'grey', lwd = 1) +
    geom_point(size = 3,color='black',shape=21) +

    stat_compare_means(aes(group = group),
                       paired = TRUE,
                       method = "t.test",
                       label = "p.format",
                       label.y =max(cpm_mat$log2cpm+0.1),
                       bracket.size = 1)+
    geom_boxplot()+scale_fill_manual(values = c("#ff6666","#7f7f7f"))+
    theme_bw()+
    theme(axis.text.x = element_text(face = "bold",angle = 45, hjust=1, vjust=1,size = 14), 
          legend.position= "none",
          strip.background = element_rect(color="black",fill =  NA, linetype = "solid",size=1),
          strip.text = element_text(face = "bold", color = "black",hjust = 0.5, size = 18,),
          plot.title=element_text(face="bold.italic",size="20", color="brown",hjust = 0.5),
          axis.text.y = element_text(face = "bold",size = 14) ,
          axis.title = element_text(size = 20,face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black", size=1, linetype="solid")
    )+labs(x = pro)+ylab(expression(Expression-log[2](CPM+1)))+ylim(min(group1$log2cpm-0.1),max(group1$log2cpm+0.1))

#Extended Data Fig. 8 e------------------ 

#Extended Data Fig. 8 f------------------ 
load('sp_correct.RData')
df<-read.table("sp_table.txt", sep='\t',header = T)
Idents(sp)<-sp$development;p<-subset(sp,development=='unknown',invert=T)
my_levels<-c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')
Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
p<-subset(sp,file=='IT2');Idents(p)<-p$development;p<-subset(p,ident='unknown',invert=T);p$development<-Idents(p)
SpatialPlot(p,alpha = c(1,0.1),ncol=1,images = df[df$file=='IT2',]$image,group.by = 'development')
SpatialFeaturePlot(p,features = c('COL1A1','ANTXR1'),alpha = c(1,1),ncol=2,images = df[df$file=='IT2',]$image)

#Extended Data Fig. 8 g------------------ 
DefaultAssay(sp)<-'spatial'
p<-subset(sp,COL1A1>0 & POSTN>0 & ANTXR1>0)
DefaultAssay(p)<-'spatial';p$group<-'POSTN+CAF';Idents(p)<-p$group
FeatureScatter(p, pt.size = 2, feature1 = 'COL1A1', feature2 = 'ANTXR1',raster = T) +
  stat_cor(method = "pearson",size=6,vjust=0.5)+
  stat_regline_equation(aes(label=paste(..eq.label..)),size=6,vjust=1.6)+
  stat_smooth(method='lm', color="black", se=)+labs(title='')+
  guides(size=FALSE)+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values ='#f37f4f')+NoLegend()

#Extended Data Fig. 8 h------------------ 
my_comparison<-list(c('Nor','Hyp'),c('Hyp','MiD'),c('MiD','MoD'),
                    c('MoD','SD&CA'),c('SD&CA','ICA'),c('ICA','MCA'))
VlnPlot(p,features = c('ANTXR1'),pt.size = 0,assay = 'spatial',cols = mycolor)+
  geom_signif(comparisons = my_comparison,y_position = 2.3,textsize=4,step_increase=0)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('ANTXR1 expression')


#Supplementary Fig. 10 -----
#Supplementary Fig. 10  a b c d -----

library(Seurat)
load("sp_correct.RData")
library(xlsx)
library(dplyr)
bcells <- xlsx::read.xlsx("Table S2.xlsx",sheetName = "Sheet4")
tcells <-  xlsx::read.xlsx("Table S2.xlsx",sheetName = "Sheet5")
cd4 <-  xlsx::read.xlsx("Table S2.xlsx",sheetName = "Sheet6")
cd8 <-  xlsx::read.xlsx("Table S2.xlsx",sheetName = "Sheet7")
mye <-  xlsx::read.xlsx("Table S2.xlsx",sheetName = "Sheet8")
endo <-  xlsx::read.xlsx("Table S2.xlsx",sheetName = "Sheet9")
fibro <-  xlsx::read.xlsx("Table S2.xlsx",sheetName = "Sheet10")
load("modules_meta_0702.RData")
modules_ct = split(modules$sub,modules$module)
modules_ct
all_deg <- rbind(bcells,cd4,cd8,mye,endo,fibro)
top5 = all_deg %>% 
  group_by(cluster) %>% 
  top_n(5,avg_log2FC)
m1  = top5 %>% 
  filter(cluster%in% modules_ct$G1)
m1
m2  = top5 %>% 
  filter(cluster%in% modules_ct$G2)
m3 = top5 %>% 
  filter(cluster%in% modules_ct$G3)

m4 = top5 %>% 
  filter(cluster%in% modules_ct$G4)
m5 = top5 %>% 
  filter(cluster%in% modules_ct$G5)
m6 = top5 %>% 
  filter(cluster%in% modules_ct$G6)
m7 = top5 %>% 
  filter(cluster%in% modules_ct$G7)
genes_to_check <- list(
  T_cells = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A'),
  NK = c('KLRB1','NCR1',"FGFBP2", "CX3CR1","NCAM1","NKG7","GNLY"),
  B_cell=c("CD79A", "SLAMF7", "BLNK", "FCRL5",'IGHG1', 'MZB1', 'SDC1'),
  Fibro=c('DCN','SFRP4','LUM','COL1A1','PCLAF'),
  Endo = c('PECAM1', 'VWF'),
  Myeloid = c(  'MMP19',"CD163", "CD14","FCGR3A","FCGR3B",
             "CD68", 'C1QA',  'C1QB',
             'TPSAB1' , 'TPSB2',
             "MNDA","CSF3R","S100A9",
             'XCR1','CLEC9A',"CD1C","LAMP3","LILRA4")
)
module_genge <- list(
  m2 = m2$gene,
  m3 = m3$gene,
  m4 = m4$gene,
  m5 = m5$gene,
  m6 = m6$gene,
  m7 = m7$gene
)
colnames(sp@meta.data)

gene<-c("CXCL13","ACP5","LAG3","PHLDA1","HAVCR2","RGS2","FOXP3","TIGIT","BATF","TNFRSF18","TNFRSF4","TNFRSF9","COL1A1","COL3A1","COL1A2","SPARC","FN1","POSTN","PLVAP","COL4A1","COL4A2","HSPG2","VWF","IGFBP7","SPP1","APOC1","MMP12","MMP9","FBP1","APOE")
p<-subset(sp,development=='unknown',invert=T)
Idents(p)<-p$development
p<-AddModuleScore(p,features = list(gene),name = 'MIMER')
length(colnames(p@meta.data))
colnames(p@meta.data)[84]<-'MIMER'
my_levels<-c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')
Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor','Hyp'),c('Hyp','MiD'),c('MiD','MoD'),c('MoD','SD&CA'),c('SD&CA','ICA'),c('ICA','MCA'))

library(ggpubr)
p<-AddModuleScore(p,features = module_genge)
length(colnames(p@meta.data))
colnames(p@meta.data)
colnames(p@meta.data)[c(85:90)] = names(module_genge)
VlnPlot(p,features = c('MIMER'),pt.size = 0,assay = 'spatial',cols = mycolor)+
  geom_signif(comparisons = my_comparison,y_position = 1.4,textsize=4,step_increase=0,hjust=0,angle=45)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('MIMER expression')+ylim(-1,2.5)


lapply( names(module_genge), function(pro){
 
  VlnPlot(p,features =pro,pt.size = 0,assay = 'spatial',cols = mycolor)+
    geom_signif(comparisons = my_comparison,
                # y_position = 1.4,
                textsize=4,step_increase=0.00,angle=45,hjust=0)+NoLegend()+
    theme(axis.title.x = element_blank())+ggtitle("")+ylab(paste0(pro,' expression'))+ylim(-1,2.5)


  })
p<-AddModuleScore(p,features = genes_to_check)
length(colnames(p@meta.data))
colnames(p@meta.data)
colnames(p@meta.data)[c(91:96)] = names(genes_to_check)
lapply( names(genes_to_check), function(pro){
  VlnPlot(p,features =pro,pt.size = 0,assay = 'spatial',cols = mycolor)+
    geom_signif(comparisons = my_comparison,
                # y_position = 1.4,
                textsize=4,step_increase=0.00,angle=45,hjust=0)+NoLegend()+
    theme(axis.title.x = element_blank())+ggtitle("")+ylab(paste0(pro,' expression'))+ylim(min(p@meta.data[,pro])-0.5,
                                                                                           max(p@meta.data[,pro])+1.5)
  
})

######### ME #################
p<-subset(sp,ME=='unknown',invert=T)
Idents(p)<-p$ME
p<-AddModuleScore(p,features = list(gene),name = 'MIMER')
length(colnames(p@meta.data))
colnames(p@meta.data)[84]<-'MIMER'
my_levels<-c('Nor-ME','Hyp-ME','MiD-ME','MoD-ME','SD&CA-ME','ICA-ME','MCA-ME')
Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor-ME','Hyp-ME'),c('Hyp-ME','MiD-ME'),c('MiD-ME','MoD-ME'),c('MoD-ME','SD&CA-ME'),c('SD&CA-ME','ICA-ME'),c('ICA-ME','MCA-ME'))
p$group = Idents(p)
library(ggpubr)
p<-AddModuleScore(p,features = module_genge)
length(colnames(p@meta.data))
colnames(p@meta.data)
colnames(p@meta.data)[c(85:91)] = c("group",names(module_genge))
VlnPlot(p,features = c('MIMER'),pt.size = 0,assay = 'spatial',cols = mycolor)+
  geom_signif(comparisons = my_comparison,y_position = 1.5,textsize=4,step_increase=0.00,hjust=0,angle=45)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('MIMER expression')+ylim(-1,2.5)


# seurat violin plot 
VlnPlot(p, features = "Tcells", pt.size = 0,assay = 'spatial',cols = mycolor)

#  violin plot without noise

vln_df = data.frame(PPBP =p$Tcells, cluster = p$group)
ggplot(vln_df, aes(x = cluster, y = PPBP, fill = cluster)) + geom_violin( adjust =1,trim=TRUE, scale = "width")+
  geom_signif(comparisons = my_comparison,y_position = 1.4,textsize=4,step_increase=0.00,)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('Tcells expression')+theme_classic()

lapply( names(module_genge), function(pro){

  VlnPlot(p,features =pro,pt.size = 0,assay = 'spatial',cols = mycolor)+
    geom_signif(comparisons = my_comparison,
                # y_position = 1.4,
                textsize=4,step_increase=0.00,angle=45,hjust=0)+NoLegend()+
    theme(axis.title.x = element_blank())+ggtitle("")+ylab(paste0(pro,' expression'))+ylim(-1,2.5)

  
})
p<-AddModuleScore(p,features = genes_to_check)
length(colnames(p@meta.data))
colnames(p@meta.data)
colnames(p@meta.data)[c(92:97)] = names(genes_to_check)
lapply( names(genes_to_check), function(pro){
  VlnPlot(p,features =pro,pt.size = 0,assay = 'spatial',cols = mycolor)+
    geom_signif(comparisons = my_comparison,
                # y_position = 1.4,
                textsize=4,step_increase=0.00,angle=45,hjust=0)+NoLegend()+
    theme(axis.title.x = element_blank())+ggtitle("")+ylab(paste0(pro,' expression'))+ylim(min(p@meta.data[,pro])-0.5,
                                                                                           max(p@meta.data[,pro])+1.5)
  
})

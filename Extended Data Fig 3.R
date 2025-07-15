#Extended Data Fig 3 -----
#Extended Data Fig 3b -----
library(Seurat)
path = '/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/IO60+6/'
samples = list.files("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/IO60+6/")

codex_list = lapply(samples[1:134], function(sample){
  # sample = samples[1]
  print(sample)
  if (file.exists(paste0(path,sample,"/new_codex_processor.csv"))) {

    codex.obj <- LoadAkoya(filename = paste0(path,sample,"/new_codex_processor.csv"),
                           type = "processor", fov = sample)
    colnames(codex.obj)
    codex.obj = RenameCells(codex.obj, new.names = paste0(sample,"_",colnames(codex.obj)))
    codex.obj@assays$Akoya@layers$counts[1:4,1:4]
    codex.obj <- NormalizeData(object = codex.obj, normalization.method = "CLR", margin = 2)
    codex.obj <- ScaleData(codex.obj)
    VariableFeatures(codex.obj) <- rownames(codex.obj)  # since the panel is small, treat all features as variable.
    codex.obj <- RunPCA(object = codex.obj, npcs = 20, verbose = FALSE)
    codex.obj <- RunUMAP(object = codex.obj, dims = 1:20, verbose = FALSE)
    codex.obj <- FindNeighbors(object = codex.obj, dims = 1:20, verbose = FALSE)
    codex.obj <- FindClusters(object = codex.obj, verbose = FALSE, resolution = 2, n.start = 1)
    
    }
  
  return(codex.obj)
})

names(codex_list) = samples[1:134]
save(codex_list,file = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhang zhi chao/codex/codex_seurat.Rdata")

load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/codex_seurat.Rdata")
all_features = rownames(codex_list[[1]])
codex = merge(codex_list[[1]],codex_list[-1],
              add.cell.ids= samples[1:134])
save(codex,file = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/codex_all.Rdata")


ESCC.anchors <- FindIntegrationAnchors(object.list = codex_list, anchor.features = all_features, reduction = "rpca")
ESCC.combined <- IntegrateData(anchorset = ESCC.anchors)
DefaultAssay(ESCC.combined) <- "integrated"

ESCC.combined <- ScaleData(ESCC.combined, verbose = FALSE)
ESCC.combined <- RunPCA(ESCC.combined, npcs = 30, verbose = FALSE)
ESCC.combined <- RunUMAP(ESCC.combined, reduction = "pca", dims = 1:30)
ESCC.combined <- FindNeighbors(ESCC.combined, reduction = "pca", dims = 1:30)
ESCC.combined <- FindClusters(ESCC.combined, resolution = 0.5)
seuratObj <- RunHarmony(ESCC.combined, "Studylevel2",dims.use = 1:20, verbose = TRUE)

sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",
                     dims = 1:20)
sce <- FindClusters(sce, resolution = c(0.1,0.2,0.3,0.4,0.8))
names(sce@reductions)
sce <- RunUMAP(sce,  dims = 1:20, 
               reduction = "harmony")
DimPlot(sce,reduction = "umap",label=T,group.by = "RNA_snn_res.0.4")  

DimPlot(sce,reduction = "umap",label=T,
        group.by =  'Studylevel2') 
save(sce,file = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/codex_all_harmony.Rdata")

table(sce$Akoya_snn_res.0.1)
table(sce$Akoya_snn_res.0.8)
table(sce$Studylevel2)
DimPlot(ESCC.combined,reduction = "umap",label=T,
        group.by =  'Studylevel2') 
library(openxlsx)
meta_info = read.csv("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/meta_info/codex_meta_info.csv")
head(meta_info)
sce$Sample = sce$Studylevel2
sce@meta.data %>% 
  mutate(cellID = rownames(sce@meta.data)) %>%
  left_join(meta_info,by = "Sample") %>%
  column_to_rownames("cellID") -> sce@meta.data

sce$Patient = str_split(sce$SampleID,"-",simplify = T)[,1]
DimPlot(sce,reduction = "umap",label=T,
        group.by =  'Tissue') 
DimPlot(sce,reduction = "umap",label=T,
        group.by =  'Patient') 

check_IO66 <- function(group.by = "Akoya_snn_res.0.1",seurat = sce,features  = IO66){
  IO66 = list(
    immune = c("CD68",'CD4','CD44','CD45RO','CD45','CD11c','CD8','HLA-A','HLA-DR','CD14','CD20','CD3e','CD56'),
    lymphocyte = c("CD107a",'CD57','FOXP3','CD21','CD38','GZMB','TOX','TCF-1','CD79a','CD39'),
    immune_checkpoint = c('IFNG','IDO1','PD-1','ICOS','PD-L1','VISTA','LAG3'),
    myeloid = c('iNOS','CD66','MPO','CD163','CD11b','CD206','CD209'),
    cc = c('HistoneH3-pSer28','PCNA','Ki67'),
    collagen_actin =c("CD31",'E-cadherin','SMA','Vimentin','CD34','Beta-actin','Podoplanin','CollagenIV','Caveolin','b-Catenin1'), #CTNB1
    tumor = c('Keratin14','SOX2','Pan-Cytokeratin','Keratin8-18',"ER",'Bcl-2','EpCAM','GP100','TP63','Keratin5'),
    other = c("Osteopontin",'PLVAP-PV-1','COL-1','Periostin','ISG15','OX40')
    
  )
  # setdiff(rownames(sce),unlist(IC66))
  p <- DotPlot(sce, assay = "Akoya", features =IO66, 
               cols = c("#bababa", "brown"),
               group.by = group.by#,
               # scale = F,
               # col.min = 0
  )+
    #scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(position = "bottom") +
    scale_fill_manual( breaks = rev(names(colors)),values = colors) +
    theme(plot.margin=margin(t = 0, b = 0, unit="pt"),
          plot.subtitle = element_text(family = "serif",  colour = "gray0"),
          #plot.background = element_rect(fill = "aliceblue"),
          plot.title = element_text(face = "bold", family = "serif", hjust = 0.5, vjust = 0, colour = "black"),
          #  背景板
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          # 坐标轴标题
          axis.title.x = element_text(family = "serif",face = "bold",colour = "chocolate4", hjust = 0.5, vjust = 0),
          axis.title.y = element_text(family = "serif",face = "bold",colour = "chocolate4", hjust = 0.5, vjust = 0.5),
          # 坐标轴刻度
          axis.text.x = element_text(family = "serif", vjust = 0.5, angle = 90),
          axis.text.y = element_text(family = "serif" ),
          axis.ticks = element_line(colour = "gray0", size = 0.9, linetype = "blank"),
          # 图例
          legend.title = element_text(hjust = 0.55,face = "bold", colour = "chocolate4",family = "serif"),
          legend.direction = "vertical",
          legend.text = element_text(family = "serif"),
          legend.key = element_rect(fill = "aliceblue"),
          #legend.background = element_rect(fill = "aliceblue")
    )
  
  print(p)
  
}

check_IO66(seurat = sce,group.by = "Akoya_snn_res.0.1")
VlnPlot(sce,features = IO66[[1]],pt.size = 0,sort=T)

VlnPlot(sce,features = c("CD45","CD45RO",'CD4',"CD8"),
        group.by = "Akoya_snn_res.0.2",
        pt.size = 0,sort=T,ncol = 2)
FeaturePlot(sce,features = c("CD45","CD45RO"),order = T)

####CD4 19
####CD8 10
head(sce@meta.data)
VlnPlot(sce,features = c("nCount_Akoya" ,"nFeature_Akoya"),pt.size = 0,group.by = "Sample")
qc.Count1 <- sce@meta.data$nCount_Akoya < 1000
table(qc.Count1)
# drop low-quality1, produce clean seurat
cells.to.keep <- rownames(sce@meta.data)[qc.Count1]
sub <- subset(sce, cells=cells.to.keep)

DimPlot(sce,
        group.by = "Akoya_snn_res.0.1",label = T,repel = T)
FeaturePlot(sub,features = c("CD3e","CD4"),order = T)

VlnPlot(sub,features = c("CD45","CD3e",'CD4',"CD8","FOXP3","PD-1"),
        group.by = "Akoya_snn_res.0.2",
        pt.size = 0,sort=T,ncol = 2)

sub
sub = AddModuleScore(sub,features = IO66,ctrl = 2)
head(sub@meta.data)
colnames(sub@meta.data)
# [1] "orig.ident"                "nCount_Akoya"              "nFeature_Akoya"            "cell_id"                  
# [5] "Studylevel1"               "Studylevel2"               "Name"                      "Image"                    
# [9] "LayerData"                 "x"                         "y"                         "Type"                     
# [13] "Objectinfo.LabelType"      "Objectinfo.LabelName"      "Objectinfo.ROIType"        "Objectinfo.ROIName"       
# [17] "Objectinfo.Envelopeleft"   "Objectinfo.Envelopetop"    "Objectinfo.Enveloperight"  "Objectinfo.Envelopebottom"
# [21] "Akoya_snn_res.2"           "seurat_clusters"           "Akoya_snn_res.0.5"         "Akoya_snn_res.0.1"        
# [25] "Akoya_snn_res.0.2"         "Akoya_snn_res.0.3"         "Akoya_snn_res.0.4"         "Akoya_snn_res.0.8"        
# [29] "Sample"                    "SampleID"                  "Tissue"                    "Pathological"             
# [33] "Patient"                   "Cluster1"                  "Cluster2"                  "Cluster3"                 
# [37] "Cluster4"                  "Cluster5"                  "Cluster6"                  "Cluster7"                 
# [41] "Cluster8"    
colnames(sub@meta.data)[34:41] = names(IO66)
VlnPlot(sub,features = names(IO66),group.by = "Akoya_snn_res.0.2",pt.size = 0,sort = T,ncol = 2)
ggsave(filename = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/Figure/IO66_122_Akoya_snn_res.0.2.IO66_VlnPlot.pdf",
       width = 12,height=20 )
deg  = FindAllMarkers(sub,group.by = "Akoya_snn_res.0.2",only.pos = T)
deg
deg %>% group_by(cluster) %>%
  top_n(5,avg_log2FC) -> top5

deg %>%
  filter(cluster == 4 )%>%
  arrange(desc(avg_log2FC))

write.xlsx(deg,file = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/codex_Akoya_snn_res.0.2_FindAllMarkers_deg.xlsx")

library(COSG,lib.loc = "/home/users/zhangzhichao/R/x86_64-pc-linux-gnu-library/4.1")
marker_cosg <- cosg(
  sce,
  groups='all',
  assay='Akoya',
  slot='data',
  mu=1,
  n_genes_user=100)
marker_cosg[["names"]][["1"]]
write.xlsx(marker_cosg,file = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/codex_Akoya_snn_res.0.2_COSG_deg.xlsx")
VlnPlot(sce,features =c("nCount_Akoya","nFeature_Akoya"),pt.size = 0)

tabcor_data<-AggregateExpression(sce,group.by="Akoya_snn_res.0.2",assays="Akoya",return.seurat=TRUE)
cor_data<-cor_data[["Akoya"]]$data
pheatmap::pheatmap(cor(cor_data,method = "spearman"))
table(sce$nCount_Akoya<500)


{
  celltype=data.frame(ClusterID=0:32,
                      celltype=0:32) 
  celltype[celltype$ClusterID %in% c(0,4,5,7,11,12,13,15,16,17,18,20,21,22,24,25,26,28,29,30),2]='Epithelium'
  celltype[celltype$ClusterID %in% c(4,15,27),2]='Epithelium'
  celltype[celltype$ClusterID %in% c( 19),2]='Tcells'
  celltype[celltype$ClusterID %in% c( 3),2]='Tcells'
  celltype[celltype$ClusterID %in% c( 23,32),2]='Tcells'
  celltype[celltype$ClusterID %in% c(  10),2]='Bcells'  
  celltype[celltype$ClusterID %in% c( 1,8,9),2]='Fibroblast' 
  celltype[celltype$ClusterID %in% c( 2,31),2]='Myeloids'
  celltype[celltype$ClusterID %in% c( 6,14),2]='Endothelium'
  
  head(celltype)
  celltype
  table(celltype$celltype)
  
  sce@meta.data$Celltype = "NA"
  for(i in 1:nrow(celltype)){
    sce@meta.data[which(sce@meta.data$Akoya_snn_res.0.2 == celltype$ClusterID[i]),'Celltype'] <- celltype$celltype[i]}
  table(sce@meta.data$Celltype)
  DimPlot(sce,group.by = "Akoya_snn_res.0.2",split.by = "Celltype",ncol = 4,label = T,repel = T)
 
  check_IO66(seurat = sce,group.by = "Celltype")
 VlnPlot(sce,features = c("E-cadherin","Vimentin","CD31","CD34"),
         ncol = 2,
         sort = T,pt.size = 0,group.by = "Akoya_snn_res.0.2")
}

sce$Sample = gsub("A","B",
                  gsub("B","C",
                       gsub("C","D",
                            gsub("D","E",
                                 gsub("E","F",
                                      gsub("F","G",
                                           gsub("G","H",
                                                gsub("H","I",
                                                     gsub("I","J",sce$Studylevel2)))))))))

sce$SampleID = NULL
sce$Tissue = NULL

meta_info = read.csv("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/meta_info/codex_meta_info.csv")
sce@meta.data %>% 
  mutate(cellID = rownames(sce@meta.data)) %>%
  left_join(meta_info,by = "Sample") %>%
  column_to_rownames("cellID") -> sce@meta.data
table(sce$Pathological.y)
del = c("B04","B05","H02",'I02',"J03","C04","D04","E04","F04","G04","H05","J05",
        "H06","I06","F07","G07","D08","H08","F09","G09",'G10',"G11","H11","I11",
        "D12","E12","F12","G12","E13",'J13',"B14",'B16','C16','D16','E16','G16',
        'H16','E17',"B18","D18","E18",'F18',"G18","H18"
        )
sub = subset(sce,Sample %in% del,invert =T)
table(sub$Pathological.y)
head(sub$Tissue) 
for (ct in unique(sce$Celltype)) {
  sub = subset(sce,Celltype==ct)
  save(sub,file = paste0("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/codex/sub_celltype/",ct,'/',ct,"_codex.Rdata"))
  
}

# Tcells -----
load("~/ESCC_codex/sub_celltype/Tcells/Tcells_codex.Rdata")
sub <- NormalizeData(object = sub, normalization.method = "CLR", margin = 2)
sub <- ScaleData(sub, verbose = FALSE)
sub <- RunPCA(sub, npcs = 20, verbose = FALSE)
library(harmony)
sub <- RunHarmony(sub, "Studylevel2",dims.use = 1:20, verbose = TRUE)
sub <- FindNeighbors(sub, reduction = "harmony",
                     dims = 1:20)
sub <- FindClusters(sub, resolution = c(0.1,0.2,0.3))
sub <- FindClusters(sub, resolution = c(0.5))
sub <- RunUMAP(sub,  dims = 1:20, 
               reduction = "harmony")
FeaturePlot(sub,features = "CD3e")
IO66 = list(
  immune = c("CD68",'CD4','CD44','CD45RO','CD45','CD11c','CD8','HLA-A','HLA-DR','CD14','CD20','CD3e','CD56'),
  lymphocyte = c("CD107a",'CD57','FOXP3','CD21','CD38','GZMB','TOX','TCF-1','CD79a','CD39'),
  immune_checkpoint = c('IFNG','IDO1','PD-1','ICOS','PD-L1','VISTA','LAG3'),
  myeloid = c('iNOS','CD66','MPO','CD163','CD11b','CD206','CD209'),
  cc = c('HistoneH3-pSer28','PCNA','Ki67'),
  collagen_actin =c("CD31",'E-cadherin','SMA','Vimentin','CD34','Beta-actin','Podoplanin','CollagenIV','Caveolin','b-Catenin1'), #CTNB1
  tumor = c('Keratin14','SOX2','Pan-Cytokeratin','Keratin8-18',"ER",'Bcl-2','EpCAM','GP100','TP63','Keratin5'),
  other = c("Osteopontin",'PLVAP-PV-1','COL-1','Periostin','ISG15','OX40')
  
)
all = unlist(IO66)
Idents(sub) = sub$Akoya_snn_res.0.5
avg = AverageExpression(sub,features = all,slot = "data" )
avg = avg$Akoya

immune_cells = list(
  immune = c('CD4','CD45RO','CD45','CD8','HLA-A','HLA-DR','CD3e','CD56'),
  lymphocyte = c('CD57','FOXP3','CD38','GZMB','TOX','TCF-1','CD39'),
  immune_checkpoint = c('IFNG','IDO1','PD-1','ICOS','VISTA','LAG3',"OX40","ISG15")
)
immune_cells = unlist(immune_cells)
avg = avg[rownames(avg) %in% immune_cells,]
pheatmap::pheatmap(t(scale(t(avg))),scale = 'none')

pdf('~/ESCC_codex/sub_celltype/Tcells/Tcells_heatmap_Akoya_snn_res.0.5.pdf')
pheatmap::pheatmap(t(scale(t(avg))),scale = 'none',cellwidth = 15,cellheight  = 15,)
dev.off()


VlnPlot(sub,features = c("CD3e","CD4",'FOXP3',"OX40",'CD39',"IFNG","CD44"),pt.size = 0,group.by = "Akoya_snn_res.0.5",sort = T)

VlnPlot(sub,features = c("CD3e","CD8","CD45RO","GZMB","TOX","LAG3",'PD-1',"TCF-1"),pt.size = 0,group.by = "Akoya_snn_res.0.5",sort = T)

VlnPlot(sub,features = c("CD3e","CD8","CD56"),pt.size = 0,group.by = "Akoya_snn_res.0.5",sort = T)

sub = AddModuleScore(sub,features = list(c("CD3e",'CD8',"PD-1","LAG3","TOX")),ctrl = 2)
VlnPlot(sub,features = c("Cluster1"),pt.size = 0,group.by = "Akoya_snn_res.0.5",sort = T)

DotPlot(sub,features = unique(top5$gene),group.by = "Akoya_snn_res.0.3")
FeaturePlot(sub,features = c("PD-1","TOX"),label = T,order = T)
table(sub$Akoya_snn_res.0.3,sub$Mimer_region)
FeatureScatter(sub,feature1 = "PD-1",feature2 = "TOX")

celltype=data.frame(ClusterID=0:21,
                    celltype=0:21) 
celltype[celltype$ClusterID %in% c(8,9,15),2]='CD4T-Treg' 
celltype[celltype$ClusterID %in% c(0),2]='CD4T-aTreg' 
celltype[celltype$ClusterID %in% c(6,13,20,4,14),2]='CD4Th'
celltype[celltype$ClusterID %in% c(13),2]='ISG15_CD4T'

celltype[celltype$ClusterID %in% c(11),2]='CD4TRTE'
celltype[celltype$ClusterID %in% c(17),2]='CD8T'
celltype[celltype$ClusterID %in% c(10,19,21),2]='LAG3+Tex'
celltype[celltype$ClusterID %in% c(18),2]='PDCD1+Tex'
celltype[celltype$ClusterID %in% c(1),2]='PDCD1+Tex2'

celltype[celltype$ClusterID %in% c(2,9,7,12,17),2]='CD8Tm'
celltype[celltype$ClusterID %in% c(9),2]='Tcells.cc'
celltype[celltype$ClusterID %in% c( 3,5,16),2]='NKT'

sub@meta.data$subCelltype = "NA"
for(i in 1:nrow(celltype)){
  sub@meta.data[which(sub@meta.data$Akoya_snn_res.0.5 == celltype$ClusterID[i]),'subCelltype'] <- celltype$celltype[i]}
table(sub@meta.data$subCelltype)

DimPlot(sub,group.by = "Akoya_snn_res.0.5",label = T,repel = T,raster = T)|DimPlot(sub,group.by = "subCelltype",label = T,repel = T,raster = T)+scale_color_igv()
Tcells = sub@meta.data
avg = AverageExpression(sub,features = all,slot = "data",group.by = "subCelltype")
avg = avg$Akoya
avg = avg[rownames(avg) %in% c(immune_cells,"Ki67",'PCNA',"CD31"),]
pdf('~/ESCC_codex/sub_celltype/Tcells/Tcells_heatmap.pdf')
pheatmap::pheatmap(t(scale(t(avg))),scale = 'none',cellwidth = 15,cellheight  = 15,)
dev.off()

#Bcells -----
celltype=data.frame(ClusterID=0:9,
                    celltype='Bcells') 
celltype[celltype$ClusterID %in% c(4),2]='Ki-67+Bcells' #0 3 8 9 15 
celltype[celltype$ClusterID %in% c(5),2]='ISG15+Bcells'
celltype[celltype$ClusterID %in% c(1,3,7),2]='ImmatureBcells' # CD20 hi CD38 hi
celltype[celltype$ClusterID %in% c(2),2]='matureBcells'
celltype[celltype$ClusterID %in% c(0,6),2]='activateBcells'

sub@meta.data$subCelltype = "NA"
for(i in 1:nrow(celltype)){
  sub@meta.data[which(sub@meta.data$Akoya_snn_res.0.2 == celltype$ClusterID[i]),'subCelltype'] <- celltype$celltype[i]}
# Mye -----
celltype=data.frame(ClusterID=0:13,
                    celltype=0:13) 
celltype[celltype$ClusterID %in% c(1),2]='MDSC' 
celltype[celltype$ClusterID %in% c(5),2]='ISG15+Mac' 
celltype[celltype$ClusterID %in% c(7),2]='M1like-Mac'
celltype[celltype$ClusterID %in% c(0,3,10,11),2]='M2like-Mac'
celltype[celltype$ClusterID %in% c(0),2]='SPP1+Mac'
celltype[celltype$ClusterID %in% c(2,12,13),2]='Neutrophil'
celltype[celltype$ClusterID %in% c(9),2]='Monocyte'
celltype[celltype$ClusterID %in% c(8),2]='DC'
celltype[celltype$ClusterID %in% c( 6),2]='Ki-67+Mac'

sub@meta.data$subCelltype = "NA"
for(i in 1:nrow(celltype)){
  sub@meta.data[which(sub@meta.data$Akoya_snn_res.0.3 == celltype$ClusterID[i]),'subCelltype'] <- celltype$celltype[i]}
#fibro 
celltype=data.frame(ClusterID=0:7,
                    celltype=0:7) 
celltype[celltype$ClusterID %in% c(1,2),2]='MyoFibroblast'
celltype[celltype$ClusterID %in% c(6,5,7,3,4),2]='COL1A1+CAF'
celltype[celltype$ClusterID %in% c(7),2]='ISG15+CAF'
celltype[celltype$ClusterID %in% c(0),2]='POSTN+CAF'
sub@meta.data$subCelltype = "NA"
for(i in 1:nrow(celltype)){
  sub@meta.data[which(sub@meta.data$Akoya_snn_res.0.1 == celltype$ClusterID[i]),'subCelltype'] <- celltype$celltype[i]}
#endo
celltype=data.frame(ClusterID=0:9,
                    celltype=0:9) 

celltype[celltype$ClusterID %in% c(0,2,7,6,8),2]='PLVAP+Endo' 
celltype[celltype$ClusterID %in% c(1,9),2]='Vimentin+Endo'
celltype[celltype$ClusterID %in% c(5),2]='PDPN+Endo'
celltype[celltype$ClusterID %in% c(3),2]='CollagenIV+Endo'
celltype[celltype$ClusterID %in% c(4),2]='Ki-67+Endo'
sub@meta.data$subCelltype = "NA"
for(i in 1:nrow(celltype)){
  sub@meta.data[which(sub@meta.data$Akoya_snn_res.0.2 == celltype$ClusterID[i]),'subCelltype'] <- celltype$celltype[i]}
table(sub@meta.data$subCelltype)
#epi
celltype=data.frame(ClusterID=0:21,
                    celltype= 0:21) 
celltype[celltype$ClusterID %in% c(0),2]='non-canonical_epi1' 
celltype[celltype$ClusterID %in% c(12),2]='non-canonical_epi2'
celltype[celltype$ClusterID %in% c(21),2]='non-canonical_epi3'
celltype[celltype$ClusterID %in% c(20),2]='GP100+tumor'
celltype[celltype$ClusterID %in% c(13,19,7,10,8),2]='Ki-67+cancer'
celltype[celltype$ClusterID %in% c(2,4,6,15,9,18,17,11,3,5),2]='differentiating_cancer'
celltype[celltype$ClusterID %in% c(16),2]='basal-like_epi'
celltype[celltype$ClusterID %in% c(1,14),2]='differentiating_epi'
sub@meta.data$subCelltype = "NA"
for(i in 1:nrow(celltype)){
  sub@meta.data[which(sub@meta.data$Akoya_snn_res.0.1 == celltype$ClusterID[i]),'subCelltype'] <- celltype$celltype[i]}

#Extended Data Fig 3c -----
Idents(codex) = 'subCelltype'
codex = subset(codex,subCelltype=='Epithelium',invert=T)
codex <- NormalizeData(object = codex, normalization.method = "CLR", margin = 2)
codex  = ScaleData(codex)
avg = AverageExpression(codex,layer='data')
avg = avg$Akoya
rank = rank %>% filter(subCelltype!='Epithelium')
avg = avg[,match(rank$subCelltype,colnames(avg))]

palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-4, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05,4, length.out=floor(palette_length/2)))
anno_cols = rank
names(major) = unique(anno_cols$Celltype.y)
annotation_colors = list(
  subCelltype = annotation_colors$subCelltype[names(annotation_colors$subCelltype) !="Epithelium"],
  Celltype.y = annotation_colors$Celltype.y
  
)
rownames(anno_cols)=anno_cols$subCelltype
anno_cols = anno_cols[,c(2,1)]
save(cell_color,major,file = '~/ESCC_codex/codex_color.Rdata')
pdf("~/ESCC_codex/CODEX_Average_Scaled_Celltype_Marker_Expression.pdf",width = 10,height = 15)
pheatmap::pheatmap(

  t(scale(t(avg))),
  scale = 'none',
  cluster_cols = F,
  color = my_color,
  annotation_col = anno_cols,
  annotation_colors = annotation_colors,
  # display_numbers = T,
  breaks = my_breaks,
  cellwidth = 8,
  cellheight = 8,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "CODEX Average Scaled Celltype Marker Expression",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean")
dev.off()
#Extended Data Fig 3d -----
load("/BGFS1/home/zhangzc/ESCC_codex/cor/scRNA_avg.Rdata")
p2r = openxlsx::read.xlsx('/BGFS1/home/zhangzc/ESCC_codex/cor/rna_protein_df_edit.xlsx')

head(p2r)
table(p2r$Protein %in% rownames(codex_avg))
p2r$Protein[p2r$Protein %in% rownames(codex_avg)]
setdiff(p2r$Protein,rownames(codex_avg))
setdiff(rownames(codex_avg),p2r$Protein)
p2r$Protein = gsub("Pan.Cytokeratin","Pan-Cytokeratin",
                   gsub("HistoneH3.pSer28","HistoneH3-pSer28",
                        gsub("Keratin8.18","Keratin8-18",
                             gsub("b-Catenin","b-Catenin1",
                                  gsub('COL.1',"COL-1",
                                       gsub("PLVAP.PV.1","PLVAP-PV-1",
                                            gsub("Keratin8.18","Keratin8-18",
                                                 gsub("E.cadherin","E-cadherin",
                                                      gsub("Bcl.2","Bcl-2",p2r$Protein)))))))))
scRNA_avg
scRNA_avg %>%
  mutate(RNA = rownames(.)) %>%
  left_join(p2r,by = "RNA") %>%
  filter(RNA !='KRT18') %>%
  column_to_rownames("Protein") %>%
  dplyr::select(-RNA) -> scRNA_avg2
all(rownames(scRNA_avg2)==rownames(codex_avg))
setdiff(rownames(scRNA_avg2),rownames(codex_avg))
setdiff(rownames(codex_avg),rownames(scRNA_avg2))
scRNA_avg2 = scRNA_avg2[match(rownames(codex_avg),rownames(scRNA_avg2)),]
scRNA_avg3 = t(scale(t(scRNA_avg2)))
codex_avg2 = t(scale(t(codex_avg)))
num_cols_scRNA <- ncol(scRNA_avg3)
num_cols_codex <- ncol(codex_avg2)
cross_correlation_matrix <- matrix(NA, nrow = num_cols_scRNA, ncol = num_cols_codex)
rownames(cross_correlation_matrix) <- colnames(scRNA_avg3)
colnames(cross_correlation_matrix) <- colnames(codex_avg2)

for (i in 1:num_cols_scRNA) {
  for (j in 1:num_cols_codex) {
    # 计算 matrix_A 的第 i 列和 matrix_B 的第 j 列之间的相关性
    cross_correlation_matrix[i, j] <- cor(scRNA_avg3[, i], codex_avg2[, j],method = 'pearson')
  }
}

scRNA_cols = rep('grey',num_cols_scRNA)
names(scRNA_cols) =  colnames(scRNA_avg3)
scRNA_cols[names(scRNA_cols) %in% c("CD4_C7_OX40","CD8_C6_CD39",'FB_C3_COL1A1','Mac_C2_SPP1','Endo_C3_RGCC')] = 'brown'
codex_cols = rep('grey',num_cols_codex)
Mimer = c("SPP1+Mac",'ISG15+Mac','Ki-67+Mac',
          'CD4T-aTreg','Tcells.cc',
          'PDCD1+Tex','PDCD1+Tex2','LAG3+Tex',
          'COL1A1+CAF','POSTN+CAF',
          'PLVAP+Endo','Ki-67+Endo')
names(codex_cols) =  colnames(codex_avg2)
codex_cols[names(codex_cols) %in% Mimer] = 'brown'
anno_cor_colors = list(
  codex_mimer = codex_cols,
  scRNA_mimer = scRNA_cols
)
anno_cols_codex = data.frame(
  row.names = names(codex_cols) ,
  codex_mimer = names(codex_cols)
)
anno_cols_scRNA = data.frame(
  row.names = names(scRNA_cols) ,
  scRNA_mimer = names(scRNA_cols)
)

palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-1, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05,1, length.out=floor(palette_length/2)))

pdf('scRNA_codex_cor_pheatmap.pdf',width = 8,height = 10 )
pheatmap::pheatmap(cross_correlation_matrix,scale = 'none',
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   color = my_color,
                   annotation_col = anno_cols_codex,
                   annotation_row = anno_cols_scRNA,
                   annotation_colors = anno_cor_colors,
                   breaks = my_breaks,
                   fontsize_row = 6,fontsize_col = 6,
                   cellwidth  = 6,cellheight  = 6)
dev.off()
#Extended Data Fig 3e -----

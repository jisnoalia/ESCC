#Figure 2a ------


#Figure 2c -----
DefaultAssay(sp) = "rctd_full"
celltypes = c("Epithelium","Cancer","FB.C1.CFD","FB.C2.IGF1" ,"FB.C3.COL1A1" ,"FB.C4.APOE","FB.C5.PDGFRB","FB.C6.ACTA2",
              "Endo.C1.ACKR1", "Endo.C2.FBLN5", "Endo.C3.RGCC","Endo.C4.CCL21",
              "Mac.C1.NLRP3","Mac.C2.SPP1" ,"Mac.C3.C1QC" ,"Mac.C4.LYVE1","Mac.C5.CXCL10","Mac.C6.MKI67",
              "Neutrophil",'CD8.C6.CD39','Platelet','T.C2.MKI67','DC.C3.CD1A',"CD4.C7.OX40"
)
sp_cc_sel = subset(sp,cc_15 %in% c("CC11","CC12","CC9","CC5","CC3","CC10"))

VlnPlot(sp_cc_sel,features = celltypes ,group.by = "cc_15",pt.size = 0,stack = T,flip = T)
# DotPlot(sp_cc_sel,features = celltypes ,group.by = "cc_15")
# Create heatmap
palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-1, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 1, length.out=floor(palette_length/2)))

anno_cols = data.frame(ct = rownames(t(top_acts_mat)),
                       Group = 'grey')
anno_cols$Group = ifelse(anno_cols$ct %in% c("CD4.C7.OX40","CD8.C6.CD39","Endo.C3.RGCC","FB.C3.COL1A1","Mac.C2.SPP1"),"MIMER","other")
rownames(anno_cols) = anno_cols$ct
anno_cols$ct = NULL
annotation_colors = list(
  Group = c('MIMER' = "brown","other" = 'grey')
)
pdf("Spatial_mapping_of_the_cellular_compostion_cc_15_tmp.pdf", width = 7, height = 15)
pheatmap(
  t(top_acts_mat),
  border_color = "white",
  color = my_color,
  breaks = my_breaks,
  cellwidth = 10,
  cellheight = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  annotation_row = anno_cols,
  annotation_colors = annotation_colors,
  main = "Spatial mapping of \nthe cellular compostion",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_colnames = TRUE,
  show_rownames = TRUE,
  cluster_cols = T,
  angle_col = 90
)
dev.off()

#Figure 2e -----

ligand_receptor_df = 
  data.frame(ligand  = c("ICAM1","THY1","VCAM1","LGALS9","BTN3A3","CD24","SELP","KIT","EGFR","MTDH",'ST6GALNAC1','HLA-A','CD274',"CD47","EPHA4","VEGFA"),
             receptor = c("ITGAM","ITGAM","ITGA4","CLEC7A","CLEC4G","SIGLEC10","SELPLG","KITLG","EGF","CEACAM1","SIGLEC15",'LILRB1','PDCD1',"SIRPA","EFNA1","KDR")
             )


colnames(sp@meta.data)

DefaultAssay(sp) ='spatial'
VlnPlot(sp,features = rownames(sp)[grep("^EFNA",rownames(sp))],group.by = "cc_15",pt.size = 0,ncol = 1)
table(ligand_receptor_df$receptor %in% rownames(sp))
patient = unique(sp$file)

library(tidyverse)
library(ggpubr)
for (p in patient) {
  
  s1 = subset(sp,development =='unknown',invert=T)
  # slide = df[df$file == p,"image"]
  # # sample = df[df$file == p,"image"]
  # s1@images = s1@images[slide]
  
  ST.exp = s1@assays$spatial@data %>% as.data.frame()
  spot_abundance_category = s1@meta.data %>% as.matrix()
  
  ligand_receptor_df1 =  ligand_receptor_df
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
    # i =1 
    ct = ligand_receptor_df1[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  
  ligand_receptor_df1$lrpair <- paste0(ligand_receptor_df1$ligand, "_", ligand_receptor_df1$receptor)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df1, seq(nrow(ligand_receptor_df1)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) # step 2
  rownames(LigandR_mean.m) -> turn1
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("cc_15")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "cc_15"
  data_long1 <- reshape2::melt(data.plot, id.vars = c("cellID", "cc_15"), measure.vars = ligand_receptor_df1$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long1[, "cc_15"] <- factor(data_long1[, "cc_15"])
  
  
  data_long = data_long1
  turn = gsub("\\_$","",
              paste0(ligand_receptor_df$ligand,"_",ligand_receptor_df$receptor))
  data_long$LigandReceptor = factor(data_long$LigandReceptor,
                                    levels = turn)
  # data_long$cc_15 = ifelse(!data_long$cc_15 %in% c("CC3","CC5","CC9","CC10"),'other',data_long$cc_15)
  
  data_long_fil = data_long %>%  filter(cc_15 %in% c("CC3","CC5","CC9","CC10"))
  data_long_fil$cc_15 = factor(data_long_fil$cc_15,
                               levels = c("CC3","CC5","CC9","CC10",'other'))
  all_ccs <- unique(sp$cc_15)
  color_mapping15 <- setNames(dittoSeq::dittoColors()[1:15][1:length(all_ccs)], sort(all_ccs))
  
  my_comparsion = list(
    c("CC3","CC5"),
    c("CC5","CC9"),
    c("CC9","CC10")
  )
  g <- ggplot(data_long_fil,aes(x=LigandReceptor,y= LRmean,fill=cc_15))+ 
    # geom_point()+
    geom_boxplot(width = 0.4, outliers  = F, position = position_dodge(width = 0.7)) +
    # facet_grid(~LigandReceptor) + 
    scale_fill_manual(values = color_mapping15)+
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    stat_compare_means(mapping = aes(group = cc_15),label.y = 2,show.legend = F,label = "p.signif",size=2,step.increase = 0.05)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 8,color = "black"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  g
  
  # figdir = "./Figure/cci"
  ggsave(filename = "~/workspace/esca/picture/Development_CC3_5_9_10_co_expression.pdf", width = 8, height = 3)
  head(data_long_fil)
  data_long_fil %>%
    group_by(cc_15,LigandReceptor) %>%
    summarise(mean = mean(LRmean)) %>%
    pivot_wider(names_from =  LigandReceptor ,values_from = mean) %>%
    column_to_rownames("cc_15") -> mean_data_long_fill
  pdf(file =  "~/workspace/esca/picture/Development_CC3_5_9_10_co_expression_mean_pheatmap.pdf")
  pheatmap::pheatmap(t(scale(mean_data_long_fill)),
                     scale = "none",
                     cellwidth = 15,
                     cellheight = 15)
  dev.off()  

  }

{
  s1 = subset(sp,ME =='unknown',invert=T)
  
  ST.exp = s1@assays$spatial@data %>% as.data.frame()
  spot_abundance_category = s1@meta.data %>% as.matrix()
  
  ligand_receptor_df1 =  ligand_receptor_df
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
    # i =1 
    ct = ligand_receptor_df1[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  
  ligand_receptor_df1$lrpair <- paste0(ligand_receptor_df1$ligand, "_", ligand_receptor_df1$receptor)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df1, seq(nrow(ligand_receptor_df1)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) # step 2
  rownames(LigandR_mean.m) -> turn1
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("cc_15")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "cc_15"
  data_long1 <- reshape2::melt(data.plot, id.vars = c("cellID", "cc_15"), measure.vars = ligand_receptor_df1$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long1[, "cc_15"] <- factor(data_long1[, "cc_15"])
  
  
  data_long = data_long1
  turn = gsub("\\_$","",
              paste0(ligand_receptor_df$ligand,"_",ligand_receptor_df$receptor))
  data_long$LigandReceptor = factor(data_long$LigandReceptor,
                                    levels = turn)
  
  data_long_fil = data_long %>%  filter(cc_15 %in% c("CC3","CC5","CC9","CC10"))
  data_long_fil$cc_15 = factor(data_long_fil$cc_15,
                               levels = c("CC3","CC5","CC9","CC10",'other'))
  all_ccs <- unique(sp$cc_15)
  color_mapping15 <- setNames(dittoSeq::dittoColors()[1:15][1:length(all_ccs)], sort(all_ccs))
  
  my_comparsion = list(
    c("CC3","CC5"),
    c("CC5","CC9"),
    c("CC9","CC10")
  )
  g <- ggplot(data_long_fil,aes(x=LigandReceptor,y= LRmean,fill=cc_15))+ 
    # geom_point()+
    geom_boxplot(width = 0.4, outliers  = F, position = position_dodge(width = 0.7)) +
    # facet_grid(~LigandReceptor) + 
    scale_fill_manual(values = color_mapping15)+
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    stat_compare_means(mapping = aes(group = cc_15),label.y = 2,show.legend = F,
                       # label = "p.signif",
                       size=2,step.increase = 0.05
                       ,method = 'wilcox.test'
                       )+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 8,color = "black"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  g
  
  ggsave(filename = "~/workspace/esca/picture/ME_CC3_5_9_10_co_expression.pdf", width = 8, height = 3)
  head(data_long_fil)
  data_long_fil %>%
    group_by(cc_15,LigandReceptor) %>%
    summarise(mean = mean(LRmean)) %>%
    pivot_wider(names_from =  LigandReceptor ,values_from = mean) %>%
    column_to_rownames("cc_15") -> mean_data_long_fill
  pdf(file =  "~/workspace/esca/picture/ME_CC3_5_9_10_co_expression_mean_pheatmap.pdf")
  pheatmap::pheatmap(t(scale(mean_data_long_fill)),
                     scale = "none",
                     cellwidth = 15,
                     cellheight = 15)
  dev.off()  
}

#Figure 2f -----
 DefaultAssay(sp) = 'spatial'
  df[df$file=='QN','image']
  location = GetTissueCoordinates(sp@images$TF01ST021104)
  # meta1 = t(meta1)
  count = GetAssayData(sp,slot = "data")
  count1 = count[,colnames(count) %in% rownames(location)]
  g1 = plot.LR(count1, location = location, LRpair = c('HLA-A','LILRB1'), pt.size = 0.8)
  g2 = plot.LR(count1, location = location, LRpair = c('CD47','SIRPA'), pt.size = 0.8)
  print(g1)
  
  print(g2)
  
  print(g3)

  common_xlim <- c(50, 550)
  common_xbreaks <- seq(50, 550, by = 100)
  all_ccs <- unique(sp$cc_15)
  color_mapping15 <- setNames(dittoSeq::dittoColors()[1:15][1:length(all_ccs)], sort(all_ccs))
  sp$new_mac_cc = ifelse(sp$cc_15 %in% c("CC3","CC5","CC9","CC10"),sp$cc_15,'other')

  ST.exp = sp@assays$spatial@data %>% as.data.frame()

  ligand_receptor_df1 =  ligand_receptor_df
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
    # i =1 
    ct = ligand_receptor_df1[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  ligand_receptor_df1$lrpair <- paste0(ligand_receptor_df1$ligand, "_", ligand_receptor_df1$receptor)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  
  ligandReceptor.ls <- split(ligand_receptor_df1, seq(nrow(ligand_receptor_df1)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) 
  LigandR_mean.m[1:4,1:4]
  
  LR_sel = LigandR_mean.m[rownames(LigandR_mean.m) %in% c("HLA-A_LILRB1",'CD47_SIRPA'),] %>% t() %>% as.data.frame() %>%
    mutate(cellID = rownames(.))
  all(sp$cellID==LR_sel$cellID)
  sp$HLA_A_LILRB1 = LR_sel$`HLA-A_LILRB1`
  sp$CD47_SIRPA = LR_sel$CD47_SIRPA
  sp$HLA_A_LILRB1 = ifelse(sp$new_mac_cc == 'other',NA,sp$HLA_A_LILRB1)
  sp$CD47_SIRPA = ifelse(sp$new_mac_cc == 'other',NA,sp$CD47_SIRPA)
  
  files = c('D_JT2','IT1','IT2','IT3','KT2','PT1','ST1','TT2',"QN","QP","QT1",'MLN_P','PT3LN_P_2','ZLN_P','YLN_3LN_P','QLn_PLn_N',
            'YN','YP','YT1',"YT2")
  LR_list = list()
  for (file  in files ) {
    p0 <- SpatialDimPlot(sp,group.by = 'new_mac_cc',images = df[df$file==file,'image'],stroke = NA,pt.size.factor = 1.4,image.alpha = 0)+
      scale_fill_manual(values = c(color_mapping15[ c("CC3","CC5","CC9","CC10")],'other'='lightgrey'),name='')+
      theme_minimal() +
      coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
      scale_x_continuous(breaks = common_xbreaks) +
      scale_y_continuous(breaks = common_xbreaks) +
      theme(
        aspect.ratio = 1,legend.position = 'top',legend.direction = 'horizontal')
    
    p1 <- SpatialFeaturePlot(sp, features =  "HLA_A_LILRB1", images = df[df$file==file,"image"],stroke = NA,
                             # alpha = c(1,0.6),
                             image.alpha = 0,pt.size.factor = 1.4) +
      theme_minimal() +
      scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),na.value='lightgrey',name = 'HLA-A_LILRB1'
      )+
      coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
      scale_x_continuous(breaks = common_xbreaks) +
      scale_y_continuous(breaks = common_xbreaks) +
      theme(
        aspect.ratio = 1,
        legend.title.position = "top" ,
        legend.justification = "center",
        legend.position = "top")
    
    p2 <- SpatialFeaturePlot(sp, features =  "CD47_SIRPA", images = df[df$file==file,"image"],stroke = NA,
                             # alpha = c(1,0.6),
                             image.alpha = 0,pt.size.factor = 1.4) +
      theme_minimal() +
      scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),na.value = "lightgrey"
      )+
      coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
      scale_x_continuous(breaks = common_xbreaks) +
      scale_y_continuous(breaks = common_xbreaks) +
      theme(
        aspect.ratio = 1,
        legend.title.position = "top" ,
        legend.justification = "center",
        legend.position = "top")
   p0 = p0 +labs(title = file)
   p =  p0|p1|p2 
   LR_list[[file]] = p 
   }
  
  for (file  in files ) {
  p_all = LR_list[[file]]
  ggsave(plot = p_all,filename = paste0(file,"_spatial_LR.pdf"),width = 12,height = 4)
  
  }
  
#Figure 2h -----
load('./spatial_pathology/sp_correct.RData')
df<-read.table("./RCTD/sp_table.txt", sep='\t',header = T)
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(Seurat)
human_H<-msigdbr(species = "human", category = "H")
hall<-as.data.frame(table(human_H$gs_name))
for (i in hall$Var1){
  gene<-human_H[human_H$gs_name==i,]$gene_symbol
  sp<-AddModuleScore(sp,features = list(gene),name = i)
} 
colnames(sp@meta.data)[84:133]<-as.character(hall$Var1)
sel_pathway = colnames(sp@meta.data)[c(95,109,110,94,
                         120,105,84,118,
                         96,101,87,97)]

sp@meta.data %>%
  dplyr::select(c(sel_pathway,file)) %>%
  filter(file %in% c("QP","QT1")) -> hall_score
apply(hall_score[1:12],2,fivenum) -> hall_range
rownames(hall_range) = c("min","Q1","mean","Q3","max")
plot_list = list()
for(image_name in c("QP","QT1")){
  
  for (pathway in sel_pathway) {
    p0 <- SpatialFeaturePlot(sp, features =  pathway, images = df[df$file==image_name,"image"],stroke = NA,
                       # alpha = c(1,0.6),
                       image.alpha = 0,pt.size.factor = 1.8) +
      theme_minimal() +
      scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu"))
                           ,limits = c(hall_range[1,pathway],
                                       hall_range[5,pathway])
      )+
      theme(aspect.ratio = 1,
            legend.title.position = "top" ,
            legend.justification = "center",
            legend.position = "top")
    p0[[1]][['labels']][['fill']] = gsub("HALLMARK_","",
                                         p0[[1]][['labels']][['fill']])
    plot_list[[image_name]][[pathway]] = p0
  }
  
}


wrap_plots(plot_list[["QP"]], ncol = 4) + plot_annotation(
  title = 'QP'
)
ggsave(filename = './hallmark/QP_hallmark.pdf',width = 15,height =15)
wrap_plots(plot_list[["QT1"]], ncol = 4) + plot_annotation(
  title = 'QT1'
)
ggsave(filename = './hallmark/QT1_hallmark.pdf',width = 15,height =15)
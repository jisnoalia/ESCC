#Extended  Data Fig8 -----
#Extended  Data Fig8a -----
library(scRNAtoolVis)
load("CD8_C6_CD39.Rdata");deg$cluster<-'CD8_C6_PDCD1';deg$gene<-rownames(deg);com<-deg
load("CD4_C7_OX40.Rdata");deg$cluster<-'CD4_C7_OX40';deg$gene<-rownames(deg);com<-rbind(com,deg)
load("FB_C3_COL1A1.Rdata");deg$cluster<-'FB_C3_COL1A1';deg$gene<-rownames(deg);com<-rbind(com,deg)
load("Mac_C2_SPP1.Rdata");deg$cluster<-'Mac_C2_SPP1';deg$gene<-rownames(deg);com<-rbind(com,deg)
load("Endo_C3_RGCC.Rdata");deg$cluster<-'Endo_C3_RGCC';deg$gene<-rownames(deg);com<-rbind(com,deg)


mygene <- c('TOX','PDCD1','LAYN','GEM','ENTPD1','DUSP4','CXCL13',
            'CXCL8','TGFB1','TGFB3','WNT2','WNT5A','FAP','MMP14','COL12A1','COL5A1','COL11A1',
            'COL3A1','COL1A1','MMP1','MMP11','POSTN','ANTXR1','CXCL12','TGFBR3','IL6',
            'SPP1','IL1RN','OLR1','MARCO','EREG','CCL20','MMP7','VEGFA','MMP12','MMP9','CXCL8','CCL2','AREG',
            'C1QB','C1QA','C1QC','CCL18','IL18',
            'FOXP3','TNFRSF18','TNFRSF4','TNFRSF9','TNFRSF1B','IL32','TIGIT','BATF','DUSP4','LAIR2','CTLA4','IL2RA',
            "PLVAP", "COL4A1", "COL4A2", "HSPG2", "VWF", "IGFBP7", "PECAM1", "SERPINE1", "SPARC", "INSR"
            )
mygene<-mygene[!duplicated(mygene)]
com = com[com$avg_log2FC >0,]
com$cluster <- factor(com$cluster,
                      levels = c("CD8_C6_PDCD1","CD4_C7_OX40","FB_C3_COL1A1","Endo_C3_RGCC","Mac_C2_SPP1"))
plot<-jjVolcano(diffData = com, myMarkers = mygene,size  = 4,
                fontface = 'italic',log2FC.cutoff = 0.25,
                
                tile.col = c("#5aa1a3","#20b2aa","#ba6c35","#7d5599","#c9c193"),
                aesCol = c('blue','red'),pSize = 1)
plot
#Extended  Data Fig8b -----
rm(list=ls());gc()
library(GseaVis)
library(clusterProfiler)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
library(fgsea)
library(msigdbr)
msigdbr_collections()
human_H<-msigdbr(species = "human", category = "H")
human_C2<-msigdbr(species = "human", category = "C2")
human_C5<-msigdbr(species = "human", category = "C5",subcategory = "BP")

human_C2<-human_C2[human_C2$gs_subcat%in%c('CP:BIOCARTA','CP:KEGG','CP:REACTOME'),]
human<-rbind(human_C2,human_H,human_C5)
gmt<-human %>% dplyr::select(gs_name,gene_symbol)
deg<- deg_fb
rank<-deg[,c('Gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('KEGG_ECM_RECEPTOR_INTERACTION','KEGG_FOCAL_ADHESION','REACTOME_SIGNALING_BY_MET','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX')
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_fb
rank<-deg[,c('Gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY','GOBP_LEUKOCYTE_DIFFERENTIATION','GOBP_LEUKOCYTE_MEDIATED_IMMUNITY')
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_cd8
rank<-deg[,c('gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY','GOBP_REGULATION_OF_T_CELL_DIFFERENTIATION')
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_cd4
deg<-deg_mac
rank<-deg[,c('gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_ANGIOGENESIS')
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_mac
deg<-deg_endo
rank<-deg[,c('gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c("REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS","REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS","REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES")
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_endo
#####################
gs_endo$class<-"Endo_C3_RGCC"
gs_cd4$class<-"CD4_C7_OX40"
gs_cd8$class<-"CD8_C6_CDPD1"
gs_mac$class<-"Mac_C2_SPP1"
gs_fb$class<-"FB_C3_COL1A1"

dat<-rbind(gs_cd8,gs_cd4,gs_endo,gs_fb,gs_mac)
dat<-dat[,c('ID','NES','pvalue','class')]

library(forcats)
library(ggplot2)
dat$ID <- as.factor(dat$ID)
dat$ID <- fct_inorder(dat$ID)
dat$class <- as.factor(dat$class)
dat$class <- fct_inorder(dat$class)
dat$log10p <- -log(dat$pvalue+0.001,10)

ggplot(dat, aes(class, ID)) +
  geom_point(aes(fill=NES, size=-log10(pvalue)),shape=21,color='black')+theme_bw()+
  scale_fill_gradientn(colours = rev(getPalette(10)),limits=c(1,2.7),)+

  theme(axis.text.x = element_text(size=14,angle = 90,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_text(size=10),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
ggsave(file="/realspace/project/proj_ESCC_STW_ZWM_2022_01/liuliqiu/figs1/function.pdf",height = 6,width = 9)

#Extended  Data Fig8c -----
sub = subset(sp,development == 'unknown',invert=T)
table(sub$new_mimer,sub$file) %>% as.data.frame.array() %>% as.data.frame() -> plot_df


melted_data <- plot_df %>% 
  rownames_to_column("rownames") %>%
  melt()
melted_data %>% 
  group_by(variable) %>%
  mutate(sum = sum(value),
         ratio = value/sum) -> spot_type_ratio 
spot_type_ratio %>% filter(rownames=="mimer") %>% arrange(desc(ratio)) -> rank_mimer
spot_type_ratio$variable = factor(spot_type_ratio$variable,levels = rank_mimer$variable)
spot_type_ratio

spot_type_ratio$rownames = factor(spot_type_ratio$rownames,levels = c("mimer","CD8.C6.PDCD1",'CD4.C7.OX40','FB.C3.COL1A1','Endo.C3.RGCC',"Mac.C2.SPP1","other"))
cols_mimer = c("mimer" = 'brown','other'='grey','CD4.C7.OX40' = "#1fb1aa",
               "Mac.C2.SPP1" = "#c2bc8d","CD8.C6.PDCD1" = "#569c9d",
               "Endo.C3.RGCC" = "#775095","FB.C3.COL1A1" = "#b96c35")

theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   # 图例键为空
    legend.text = element_text(color="black",size=10), # 定义图例文本
    legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
    legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
    legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
    legend.background=element_blank()) 
}
spot_type_ratio %>% 
  # filter(rownames!="other") %>%
  ggplot(.,aes(x=variable,y = ratio))+
  geom_col(aes(fill = rownames),size=2,width = 0.7#,theta = "x",start=0
  ) +
  scale_fill_manual(values = cols_mimer)+
  # coord_polar()+
  # theme_minimal()+
  theme(
    axis.title = element_blank(),
    # axis.ticks = element_blank(),
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.direction = 'vertical',
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background = element_blank()
  )+labs(title ="Epithelium/cancer dominted area",fill = 'spot type') +
  theme_niwot()-> bar
bar

sub@meta.data %>% 
  mutate(new_dev = paste0(file,'_',development)) %>%
  distinct(new_dev,.keep_all = TRUE) %>% 
  dplyr::select(file,development,new_dev) -> m1 
m1$file = factor(m1$file,
            levels = rank_mimer$variable)  
m1$development = factor(m1$development,
                        levels = c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA'))
m1 %>%
  mutate(group = "development") %>%
  ggplot(aes(file, group, fill = development)) +
  geom_tile(aes(width = 0.7, height = 0.2), position = position_dodge(width = 0.5)) +
  scale_y_discrete(expand = c(0, 0), position = "right") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = domin_cols)+
  # theme_void() +
  theme(axis.title = element_blank(),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.direction = 'vertical',
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank()) +
  theme_niwot() -> group3

library(aplot)
bar %>% insert_bottom(group3,height = 0.02) 
ggsave("spot_ratio_in_slide_dominted.pdf",width = 9,height = 4)

{
  table(sub$new_mimer,sub$development) %>% as.data.frame.array() %>% as.data.frame() -> plot_df

  melted_data <-plot_df[,1:7] %>% 
    rownames_to_column("rownames") %>%
    melt()
  melted_data %>% 
    # filter(rownames != 'other') %>%
    group_by(variable) %>%
    mutate(sum = sum(value),
           ratio = value/sum) -> spot_type_ratio 
  spot_type_ratio$rownames = factor(spot_type_ratio$rownames,levels = c("mimer","CD8.C6.PDCD1",'CD4.C7.OX40','FB.C3.COL1A1','Endo.C3.RGCC',"Mac.C2.SPP1","other"))

  spot_type_ratio %>% 
    # filter(rownames!="other") %>%
    ggplot(.,aes(x=variable,y = ratio))+
    geom_hline(
      aes(yintercept = y), 
      data.frame(y = c(0:4) * 0.25),
      color = "lightgrey"
    ) +
    geom_col(aes(fill = rownames),size=2,width = 0.7#,theta = "x",start=0
    ) +
    geom_segment(
      aes(
        x = reorder(variable, ratio),
        y = 0,
        xend = reorder(variable, ratio),
        yend = 1
      ),
      linetype = "dashed",
      color = "gray12"
    )+
    scale_fill_manual(values =cols_mimer)+
    coord_polar(start = 0)+
    annotate(
      "text",
      x = rep(0.5,5),
      y = seq(0,1,0.25),
      label = seq(0,1,0.25))+
    # coord_radial(start =0,end=1.8*pi ,r_axis_inside=TRUE)+
    theme_minimal()+
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = "gray12", size = 12),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )+labs(title = "Epithelium/cancer\ndominted area",fill = 'spot type')
  ggsave("Dominted_spot_ratio_bar_with_other.pdf",width = 5,height = 5)
#Extended  Data Fig8d -----
sub = subset(sp,ME == 'unknown',invert=T)

for (i in 1:nrow(df)) {
  SpatialDimPlot(sub,group.by = "new_mimer",
                 images = df$image[i],crop=F,image.alpha = 0.7,
                 # pt.size.factor = 1.5,
                 stroke = NA) + labs(title = df$file[i]) +
    scale_fill_manual(values = cols_mimer ) &
    guides(fill = guide_legend(override.aes = list(size=5)))
  ggsave(filename = paste0("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/Figure/Mimer/mimer_spot_vs_only/ME/",df$file[i],"_spot_single.pdf"),width = 4,height = 5)
}

sub = subset(sp,development == 'unknown',invert=T)
table(sub$new_mimer,sub$development) %>% as.data.frame.array() %>% as.data.frame() -> plot_df


melted_data <-plot_df[,1:7] %>% 
  rownames_to_column("rownames") %>%
  melt()
melted_data %>% 
  group_by(variable) %>%
  mutate(sum = sum(value),
         ratio = value/sum) -> spot_type_ratio 
spot_type_ratio$rownames = factor(spot_type_ratio$rownames,levels = c("mimer","PDCD1_only",'OX40_only','COL1A1_only','RGCC_only',"SPP1_only","other"))
spot_type_ratio %>% 
  filter(rownames!="other") %>%
  ggplot(.,aes(x=variable,y = ratio))+
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:3) * 0.1),
    color = "lightgrey"
  ) +
  geom_col(aes(fill = rownames),size=2,width = 0.7#,theta = "x",start=0
  ) +
  geom_segment(
    aes(
      x = reorder(variable, ratio),
      y = 0,
      xend = reorder(variable, ratio),
      yend = 0.4
    ),
    linetype = "dashed",
    color = "gray12"
  )+
  coord_polar(start = 0)+
  annotate(
    "text",
    x = rep(0.5,5),
    y = seq(0,0.4,0.1),
    label = seq(0,0.4,0.1))+
  # coord_radial(start =0,end=1.8*pi ,r_axis_inside=TRUE)+
  theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "gray12", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )+labs(title = "Epithelium/cancer\ndominted area",fill = 'spot type')
ggsave("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/Figure/Mimer/mimer_spot_vs_only/dominated_spot_ratio_bar.pdf",width = 5,height = 5)
spot_type_ratio %>% 
  # filter(rownames!="other") %>%
  ggplot(.,aes(x=variable,y = ratio))+
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:4) * 0.25),
    color = "lightgrey"
  ) +
  geom_col(aes(fill = rownames),size=2,width = 0.7#,theta = "x",start=0
  ) +
  geom_segment(
    aes(
      x = reorder(variable, ratio),
      y = 0,
      xend = reorder(variable, ratio),
      yend = 1
    ),
    linetype = "dashed",
    color = "gray12"
  )+
  scale_fill_manual(values =cols_mimer)+
  coord_polar(start = 0)+
  annotate(
    "text",
    x = rep(0.5,5),
    y = seq(0,1,0.25),
    label = seq(0,1,0.25))+
  # coord_radial(start =0,end=1.8*pi ,r_axis_inside=TRUE)+
  theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "gray12", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )+labs(title = "Epithelium/cancer\ndominted area",fill = 'spot type')
ggsave("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/Figure/Mimer/mimer_spot_vs_only/dominated_spot_ratio_bar_with_other.pdf",width = 5,height = 5)
#Extended  Data Fig8e -----
 p = subset(sp,development=='unknown',invert=T)
  p$new_mimer_region = ifelse(p$new_mimer %in% 'mimer','mimer','other')
  table(p$new_mimer_region,p$file) %>% as.data.frame.array() ->plot_df
  plot_df[plot_df!=0] = 1
  plot_df
  meta_plot = p@meta.data %>% dplyr::distinct(file,.keep_all = T) %>% 
    select(file,patient) %>% arrange(patient)
  rownames(meta_plot) = meta_plot$file
  plot_df = plot_df[,match(rownames(meta_plot),colnames(plot_df))]
  meta_plot$file = NULL
  anno_cols = c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", 
                "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", 
                "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", 
                "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", 
                "#D7B5A6")[1:15]
  names(anno_cols)=unique(meta_plot$patient)
  pdf("~/workspace/esca/picture/pheatmap_number_of_MIMER_spot_in_slide_development.pdf",
      width = 20 ,height = 10)
  pheatmap::pheatmap(plot_df,cluster_rows = F,cluster_cols = F,
                     cellwidth = 20,cellheight = 20,color = c("white","brown"),
                     annotation_col = meta_plot,
                     gaps_col = c(4,8,12,16,20,24,29,34,37,41,45,49,53,58),
                     annotation_colors = list(patient = anno_cols ),
                     breaks = c(0, 0.5, 1))
  dev.off()
  
  
  p = subset(sp,ME=='unknown',invert=T)
  p$new_mimer_region = ifelse(p$new_mimer %in% 'mimer','mimer','other')
  table(p$new_mimer_region,p$file) %>% as.data.frame.array() ->plot_df
  plot_df[plot_df!=0] = 1
  plot_df
  meta_plot = p@meta.data %>% dplyr::distinct(file,.keep_all = T) %>% 
    select(file,patient) %>% arrange(patient)
  rownames(meta_plot) = meta_plot$file
  plot_df = plot_df[,match(rownames(meta_plot),colnames(plot_df))]
  meta_plot$file = NULL
  anno_cols = c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", 
                "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", 
                "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", 
                "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", 
                "#D7B5A6")[1:14]
  names(anno_cols)=unique(meta_plot$patient)
  pdf("~/workspace/esca/picture/pheatmap_number_of_MIMER_spot_in_slide_ME.pdf",
      width = 20 ,height = 10)
  pheatmap::pheatmap(plot_df,cluster_rows = F,cluster_cols = F,
                     cellwidth = 20,cellheight = 20,color = c("white","brown"),
                     annotation_col = meta_plot,
                     gaps_col = c(4,8,12,16,20,24,28,33,36,40,44,48,52),
                     annotation_colors = list(patient = anno_cols ),
                     breaks = c(0, 0.5, 1))
  dev.off()
#Extended  Data Fig8f -----
 p0 <- SpatialDimPlot(sp,group.by = 'new_cc_15',images = df[df$file=='KT1_2','image'],stroke = NA,pt.size.factor = 1.4,image.alpha = 0)+
    scale_fill_manual(values = c(color_mapping15,'other'='lightgrey'),name='')+
    theme_minimal() +
    coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
    scale_x_continuous(breaks = common_xbreaks) +
    scale_y_continuous(breaks = common_xbreaks) +
    theme(
      aspect.ratio = 1,legend.position = 'top',legend.direction = 'horizontal')
  
  
  cols_mimer
  colnames(sp@meta.data)
  p1 <- SpatialDimPlot(sp,group.by = 'new_mimer',images = df[df$file=='KT1_2','image'],stroke = NA,pt.size.factor = 1.4,image.alpha = 0)+
    scale_fill_manual(values = cols_mimer,name='')+
    theme_minimal() +
    coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
    scale_x_continuous(breaks = common_xbreaks) +
    scale_y_continuous(breaks = common_xbreaks) +
    theme(
      aspect.ratio = 1,legend.position = 'top',legend.direction = 'horizontal')
  p1
#Extended  Data Fig8g -----
load("KL.RData")
load("class.RData")
brain_cell_class_new
brain_sgraph_KL_mst_cons
# 计算均值
require(purrr)
brain_sgraph_KL_mst_cons = reduce(brain_sgraph_KL_mst_cons, `+`) / length(brain_sgraph_KL_mst_cons)

library(tidyverse)
brain_cell_class_new<-brain_cell_class_new %>% reduce(inner_join, by = "id") 
brain_cell_class_new$freq.Freq<-rowMeans(brain_cell_class_new[,2:46])

long <- reshape2::melt(brain_sgraph_KL_mst_cons,id.vars= rownames(brain_sgraph_KL_mst_cons))
long
mst_cons_am <- brain_sgraph_KL_mst_cons
mst_cons_node <- data.frame(id=rownames(mst_cons_am), label=rownames(mst_cons_am))
directed = FALSE
if (!directed) mst_cons_am[upper.tri(mst_cons_am, diag = T)] <- NA

mst_cons_am <- data.frame(id=rownames(mst_cons_am), mst_cons_am, check.names=F)
mst_cons_edge <- reshape2::melt(mst_cons_am) %>% na.omit() %>% magrittr::set_colnames(c('from', 'to', 'value'))
mst_cons_edge
#节点数据
nodes <- data.frame(name = unique(union(mst_cons_edge$from, mst_cons_edge$to)))

nodes$number = brain_cell_class_new[brain_cell_class_new$id %in% nodes$name,"freq.Freq"]
nodes
class(nodes)
#边数据
rownames(mst_cons_edge) <- 1:nrow(mst_cons_edge)
edges <- mst_cons_edge
colnames(edges) <- c("from","to","weighted")
class(edges)
edges
edges = edges[!edges$weighted == 0.00,]
library(tidygraph)
igraph::graph_from_data_frame(edges, vertices = nodes) %>% as_tbl_graph() -> g
gr1_layout2 <- create_layout(g, layout = "kk")
gr1_layout2
gr1_layout2[1,1:2] <- c(0.9,0.1)
gr1_layout2[2,1:2] <- c(-0.8,0.3)
gr1_layout2[3,1:2] <- c(0.15,0)
gr1_layout2[4,1:2] <- c(0.2,1)
gr1_layout2[5,1:2] <- c(-0.13,-0.7)

use_color=cell_color[cell_color$ct %in% gr1_layout2$name,"color"]

p<-ggraph(gr1_layout2) +
  scale_color_manual(values = use_color) +
  geom_edge_hive(aes(width = weighted),colour ="#6ea6cd")+
  geom_node_point(aes(size = number,colour = name))+ 
  scale_size(range = c(5,20),breaks=seq(50,450,100),limits = c(50,450))+
  geom_node_text(aes(label=name),size=3) +
  scale_edge_colour_gradientn(
    colours =  c("#4575b4","#f1f9d8","#f0663f")
  )+   guides(size= guide_legend(ncol = 2))+
  theme_graph() + expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
#Extended  Data Fig8h -----
#leftz
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% dplyr::inner_join(lr_network %>% dplyr::distinct(from,to), by = c("from","to"))
# filter out ligands not in ligand_target_matrix
ligands = lr_network$from %>% unique()
ligands = intersect(ligands, colnames(ligand_target_matrix))
receptors = lr_network$to %>% unique()
lr_network <- lr_network %>% filter(from %in% ligands & to %in% receptors) 
## receiver
xx<-subset(sc,tissue=='Tu');Idents(xx)<-xx$L4_C
receiver = "FB_C3_COL1A1"
expressed_genes_receiver = get_expressed_genes(receiver, xx, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Mac_C2_SPP1")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, xx, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
# DEG
my<-subset(xx,L3_C=='Fibroblast')
DE_table_receiver = FindMarkers(object = my, ident.1 = "FB_C3_COL1A1", min.pct = 0.10,only.pos = TRUE) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
# define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
# ligand activity
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
# activate target gene
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
save(vis_ligand_target,file = "Mac_CAF_NICHENET.RData")
new<-vis_ligand_target[c('SPP1','MMP9','IL1RN','TGFB1'),1:33]
new %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text = element_text(size=14,colour = "black")) + scale_fill_gradientn(colors =rev(getPalette(10)))

#right
## receiver
xx<-subset(sc,tissue=='Tu');Idents(xx)<-xx$L4_C
receiver = "Mac_C2_SPP1"
expressed_genes_receiver = get_expressed_genes(receiver, xx, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("FB_C3_COL1A1")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, xx, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
# DEG
my<-subset(xx,L3_C=='Macrophage')
DE_table_receiver = FindMarkers(object = my, ident.1 = "Mac_C2_SPP1", min.pct = 0.10,only.pos = TRUE) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
# define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
# ligand activity
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
# activate target gene
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

new<-vis_ligand_target[c('TGFB3','COL5A3','COL1A1','CCL11','IGF2'),4:36]
new %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text = element_text(size=14,colour = "black")) + scale_fill_gradientn(colors =rev(getPalette(10)),breaks=c(0,0.002,0.004))
#Extended  Data Fig8i -----
VlnPlot(sp,features = "MIMER_signature_score",group.by = "cc_15",sort = T,pt.size = 0) +
  scale_color_manual(values = color_mapping15) +
  scale_fill_manual(values = color_mapping15)
#Extended  Data Fig8j -----
sp@meta.data$cc_15_combine = sp@meta.data$cc_15
sp@meta.data$cc_15_combine = ifelse(sp$cc_15_combine %in% c("CC4","CC7","CC1","CC6","CC5","CC14"),
                                    "Mimer_resident_niches",sp$cc_15_combine)

DimPlot(sp,group.by='cc_15_combine',raster=T,reduction = "umap.ischia15",label = T,repel = T)
development = subset(sp,development =='unknown',invert=T)
SpatialDimPlot(development,group.by  = c("development"),images = df[df$file=='QP','image'],stroke = NA,crop = F,pt.size.factor = 1.3)+
   scale_fill_manual(values = c("SD&CA"="#66A61E"))
ggsave(filename = "development_QP_SpatialDimPlot.pdf",width = 6,height = 6)


development@meta.data %>%
  group_by(file,development,cc_15_combine) %>% 
  summarise(number = n()) %>%
  ungroup() %>% 
  group_by(file,development) %>%
  mutate(sum = sum(number),
         ratio = number/sum) -> dev_df

head(dev_df)
library(ggrastr)
my_levels<-c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')
# Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor','Hyp'),c('Hyp','MiD'),c('MiD','MoD'),c('MoD','SD&CA'),c('SD&CA','ICA'),c('ICA','MCA'))
dev_df$development = factor(dev_df$development,
                            levels = my_levels)

dev_df %>% 
  filter(file != 'N3NT') %>% 
  filter(cc_15_combine =="Mimer_resident_niches") %>%
  ggplot(.,aes(x=development,y=ratio,color=development))+
  geom_boxplot()+
  geom_jitter(width=0.2)+ 
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  stat_compare_means(comparisons =my_comparison )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0,size = 14), 
        legend.position= "none",
        strip.background = element_rect(color="black",fill =  "#76B7B2", linetype = "solid",size=1),
        strip.text = element_text( color = "black",hjust = 0.5, size = 18,),
        plot.title=element_text(size="20", color="brown",hjust = 0.5),
        axis.text.y = element_text(size = 14) ,
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", size=1, linetype="solid")
  )+labs(x='',y="The proportion of Mimer resident niches in slide")

ggsave(filename = 'Mimer_resident_niches_ratio_in_slide.pdf',width = 6,height = 5)


ME = subset(sp,ME=='unknown',invert=T)
ME@meta.data %>%
  group_by(file,ME,cc_15_combine) %>% 
  summarise(number = n()) %>%
  ungroup() %>% 
  group_by(file,ME) %>%
  mutate(sum = sum(number),
         ratio = number/sum) -> dev_df

head(dev_df)
library(ggrastr)
my_levels<-c('Nor-ME','Hyp-ME','MiD-ME','MoD-ME','SD&CA-ME','ICA-ME','MCA-ME')
# Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor-ME','Hyp-ME'),c('Hyp-ME','MiD-ME'),c('MiD-ME','MoD-ME'),c('MoD-ME','SD&CA-ME'),c('SD&CA-ME','ICA-ME'),c('ICA-ME','MCA-ME'))
dev_df$ME = factor(dev_df$ME,
                            levels = my_levels)
dev_df %>%
  filter(file != 'N3NT') %>% 
  filter(cc_15_combine =="Mimer_resident_niches") %>%
  ggplot(.,aes(x=ME,y=ratio,color=ME))+
  geom_boxplot(aes(fill=ME))+geom_point()+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  stat_compare_means(comparisons =my_comparison )

dev_df %>% 
  filter(file != 'N3NT') %>% 
  filter(cc_15_combine =="Mimer_resident_niches") %>%
  ggplot(.,aes(x=ME,y=ratio,color=ME))+
  geom_boxplot()+
  geom_jitter(width=0.2)+ 
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  stat_compare_means(comparisons =my_comparison )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,size = 14,vjust = 0.5,hjust = 1), 
        legend.position= "none",
        strip.background = element_rect(color="black",fill =  "#76B7B2", linetype = "solid",size=1),
        strip.text = element_text( color = "black",hjust = 0.5, size = 18,),
        plot.title=element_text(size="20", color="brown",hjust = 0.5),
        axis.text.y = element_text(size = 14) ,
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", size=1, linetype="solid")
  )+labs(x='',y="The proportion of Mimer resident niches in slide")


ggsave(filename = 'Mimer_resident_niches_ratio_in_slide_ME.pdf',width = 6,height = 5)
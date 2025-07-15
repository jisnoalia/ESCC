#####Extended Data Fig. 2b & 2c#####
library(aplot)
meta = sp@meta.data
development = subset(meta,development!='unknown')
library(dittoSeq,lib.loc =  "/home/users/zhangzhichao/R/x86_64-pc-linux-gnu-library/4.1")
all_ccs <- unique(sp$cc_15)
color_mapping15 <- setNames(dittoColors()[1:cut][1:length(all_ccs)], sort(all_ccs))

development$development <- factor(development$development,
                                c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')) 
development %>% 
  mutate(new_group = paste0(file,"_",development)) -> development

file_level = development %>% 
  distinct(new_group,.keep_all = T) %>% 
  dplyr::select(development,file,new_group) %>%
  arrange(development)
head(file_level)
development$new_group = factor(development$new_group,
                          levels = file_level$new_group)
Cellratio <- prop.table(table(development$cc_15,development$new_group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
head(Cellratio)
colnames(Cellratio) <- c("cc_15","file","ratio") #对列名重命名
# Cellratio$development <- factor(Cellratio$development,
#                                 c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')) 
Cellratio = Cellratio[!is.na(Cellratio$file),]
ggplot(Cellratio,aes(x=file,y=ratio,fill=cc_15,stratum=cc_15))+
  geom_col(width = 0.8,color=NA)+
  # geom_flow(width=0.4,alpha=0.2,knot.pos=0)+ # knot.pos参数可以使连线变直
  scale_fill_manual(values=color_mapping15)+
  theme_classic()
Cellratio$cc_15 = factor(Cellratio$cc_15,
                         levels = paste0("CC",1:15))

theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   # 图例键为空
    legend.text = element_text(color="black",size=14), # 定义图例文本
    legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
    legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
    legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
    legend.background=element_blank())
}
# save(Cellratio,file_level,color_mapping15,file = "cc_in_development.Rdata")
ggplot(Cellratio, aes(file,ratio, fill=cc_15)) +geom_col(position = 'stack', width = 0.8)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('development')+ylab('Proportion')+scale_fill_manual(values = color_mapping15) -> bar

group3 <- file_level %>% 
  mutate(group="development") %>%
  as_tibble() %>% 
  mutate(development=as.character(development)) %>% 
  ggplot(aes(new_group,group,fill=development))+
  geom_tile()+
  scale_fill_manual(values = tissue_cols)+
  scale_y_discrete(expand = c(0,0),position="right")+
  scale_x_discrete(expand=c(0,0))+
  theme_void()+
  theme(axis.text.y=element_text(color="black",size=12))+
  theme_niwot()

bar %>% insert_top(group3,height = 0.05) 


#####Extended Data Fig. 2d #####
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
#####Extended Data Fig. 2e #####
sp$Celltype = sp$cc_15
sp$Celltype = ifelse(sp$cc_15 %in% c("CC9","CC11","CC12"),"epithelia-dominant",sp$Celltype)
sp$Celltype = ifelse(sp$cc_15 %in% c("CC3","CC5","CC10"),"tumor-dominant",sp$Celltype)
sp$Celltype = ifelse(sp$cc_15 %in% c("CC15"),"lamina propria-like",sp$Celltype)
sp$Celltype = ifelse(sp$cc_15 %in% c("CC8"),"tertiary lymphoid structure",sp$Celltype)
sp$Celltype = ifelse(sp$cc_15 %in% c("CC1","CC6"),"lymphoid aggregate",sp$Celltype)
sp$Celltype = ifelse(sp$cc_15 %in% c("CC2"),"muscularis-like",sp$Celltype)
sp$Celltype = ifelse(sp$cc_15 %in% c("CC4","CC7","CC14"),"fibroblasts-enriched",sp$Celltype)
sp$Celltype = ifelse(sp$cc_15 %in% c("CC13"),"LYVE1+ macrophage enriched",sp$Celltype)
DimPlot(sp,group.by='Celltype',raster=T,reduction = "umap.ischia15",label = T,repel = T)+
  scale_color_d3()

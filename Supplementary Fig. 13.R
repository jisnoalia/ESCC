#Supplementary Fig. 13 -----
#Supplementary Fig. 13 a -----
co=subset(sc,subset =L4_C=='CD8_C6_CD39')
FeatureScatter(co, pt.size = 2, feature1 = 'GEM', feature2 = 'ENTPD1') +stat_cor(method = "pearson",size=6,vjust=0.5)+
  stat_regline_equation(aes(label=paste(..eq.label..)),size=6,vjust=1.6)+stat_smooth(method='lm', color="black", se=TRUE)+labs(title='')+guides(size=FALSE)+theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  scale_color_manual(values ='#F37F4F')
#Supplementary Fig. 13 b -----
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

#Supplementary Fig. 13 c -----
DefaultAssay(sp)<-'spatial'
p<-subset(sp,CD8A>0&ENTPD1>0&GEM>0&PDCD1>0)
DefaultAssay(p)<-'spatial';p$group<-'CD8+CD39+';Idents(p)<-p$group
FeatureScatter(p, pt.size = 2, feature1 = 'GEM', feature2 = 'PDCD1') +stat_cor(method = "pearson",size=6,vjust=0.5)+
  stat_regline_equation(aes(label=paste(..eq.label..)),size=6,vjust=1.6)+stat_smooth(method='lm', color="black", se=TRUE)+labs(title='')+guides(size=FALSE)+theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  scale_color_manual(values ='#F37F4F')


#Supplementary Fig. 13 f -----
load('/realspace/project/proj_ESCC_STW_ZWM_2022_01/liuliqiu/spatial_pathology/sp_correct_with_rctd_20250523.RData')
cells = subset(sp,GEM >0)
cells$cellID
sp$GEM_spot = ifelse(sp$cellID %in% cells$cellID,"GEM+ Cells",'other Cells')
table(sp$GEM_spot,sp$new_mimer,sp$file)

table(cells$new_mimer,cells$file)
pdcd1_pos = sp@meta.data %>% filter(new_mimer %in% c("mimer","PDCD1_only")) %>% dplyr::select(cellID)
intersect(pdcd1_pos$cellID,cells$cellID)->over

sp$GEM_spot = ifelse(sp$cellID %in% over,"GEM+PDCD1+CD8+ T",sp$GEM_spot)
table(sp$GEM_spot,sp$file)

FetchData(sp,c('GEM','file')) %>% filter(GEM>0) ->gem_df
gem_df %>%
  summarise(mean = mean(GEM)) 
gem_df$GEM_pos = ifelse(gem_df$GEM > 0.89 ,"GEM_over",'GEM+')
gem_df %>% filter(GEM_pos == 'GEM+') -> GEM_over

cols_mimer = c("mimer" = 'brown','other'='grey','OX40_only' = "#1fb1aa",
               "SPP1_only" = "#c2bc8d","PDCD1_only" = "#569c9d",
               "RGCC_only" = "#775095","COL1A1_only" = "#b96c35")

sp$GEM_spot = ifelse(sp$cellID %in% rownames(GEM_over) ,'other Cells',sp$GEM_spot)
p0<-SpatialDimPlot(sp,group.by = "new_mimer",images = df[df$file =='D_JT2','image'],stroke = NA,crop = FALSE,pt.size.factor = 1.3,alpha = 0.7)+
  scale_fill_manual(values = cols_mimer)
p1 <- SpatialDimPlot(sp,group.by = "GEM_spot",images = df[df$file =='D_JT2','image'],stroke = NA,crop = FALSE,pt.size.factor = 1.3,alpha = 0.7)+
  scale_fill_manual(values = c("GEM+PDCD1+CD8+ T" = '#008080',
                               "GEM+ Cells"= '#FFD700',
                               'other Cells' = 'grey'))#+labs(title = 'GEM over (cutoff-1.10: upper-hinge in D_JT2)\nGEM+ (GEM expression >0)')
p0|p1
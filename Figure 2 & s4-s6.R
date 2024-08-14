library(psych)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(Seurat)
library(dplyr)
require("DESeq2")
# Figure 2 a -----
####### calculate cell num ratio and abundance
load('scRNA.RData');
colnames(sc@meta.data)
table(sc$tissue)
primary_all = subset(sc,tissue=='Tu')
# normalize each sample by all cells in that sample, 
# in that dataset and condition
deseq2_norm_single = function(data){
  
  a.mtx = data %>% 
    dplyr::select(sample_id, sub, sub_num) %>%
    pivot_wider(.,names_from = sub, values_from = sub_num) %>%
    replace(., is.na(.),0) %>%
    column_to_rownames(var='sample_id')     %>%
    as.matrix() %>%
    t()
  
  sf <- estimateSizeFactorsForMatrix(a.mtx, type ="poscounts")
  a.norm <- log2(sweep(a.mtx,2,sf,"/")+1)
  
  a.normdf = a.norm %>% 
    t() %>%
    as.data.frame() %>%
    as_tibble(.,rownames='sample_id') %>%
    pivot_longer(!sample_id, names_to = 'sub',values_to = 'abundance')%>%
    mutate(`sf.DESeq` = sf[sample_id])
  return(a.normdf)
} 

load('scRNA_cd4.RData');
my@meta.data %>% select(c(subC)) -> CD4
load('scRNA_cd8.RData');
my@meta.data %>% select(c(subC)) -> CD8
load('scRNA_b.RData');
my@meta.data %>% select(c(subC)) -> Bcells
load('scRNA_endo.RData');
my@meta.data %>% select(c(subC)) -> Endo
load('scRNA_fb.RData');
my@meta.data %>% select(c(subC)) -> Fibro
load('scRNA_myeloid.RData');
my@meta.data %>% select(c(subC)) -> Myeloid
CD4$major = 'CD4'
CD8$major = 'CD8'
Bcells$major = "Bcell"
Endo$major = "Endo"
Fibro$major = "Fibro"
Myeloid$major = "Myeloid"
df = rbind(CD4,CD8,Bcells,Endo,Fibro,Myeloid)
df$cellID = rownames(df)
primary_all$cellID = rownames(primary_all@meta.data)
phe = primary_all@meta.data %>% select(c(cellID,loc,metastasis,tissue,patient))
df2 = df %>% left_join(phe,by = "cellID")
df2 = df2 %>% filter(tissue=="Tu")
rownames(df2) = df2$cellID
my$tissue = stringr::str_split(my$orig.ident,"_",simplify = T)[,2]
sub$cellID = rownames(sub@meta.data)
df2$newC = df2$subC

df2 = df2 %>% filter(newC!="Smooth muscle cell")

dat<-as.data.frame(table(df2$patient,df2$newC));
colnames(dat)<-c('sample_id','sub','sub_num');
dat<-deseq2_norm_single(dat)

# ratio -----
t<-as.data.frame(table(df2$patient,df2$newC));colnames(t)<-c('sample_id','sub','sub_num')
d<-as.data.frame(table(df2$newC,df2$major));colnames(d)<-c('sub','major','freq');
d<-d[d$freq>0,];t<-merge(t,d[,1:2],by='sub')
total<-as.data.frame(table(df2$patient,df2$major));
colnames(total)<-c('sample_id','major','major_num');t<-merge(t,total,by=c('sample_id','major'))
t$ratio<-t$sub_num/t$major_num
dat<-merge(dat,t[,c('sample_id','sub','sub_num','ratio')],by=c('sample_id','sub'))
############### Imm ALL R+NR+U ######################################
samples = dat %>%
    mutate(sub = factor(sub, levels=dat$sub %>% unique()))
samples[is.na(samples$ratio),'ratio']<-0

abun = samples %>% 
    dplyr::select(sub,sample_id,ratio  ) %>%
    pivot_wider(names_from = sub, values_from =ratio ) %>%
    column_to_rownames(var='sample_id')
# cor -----
res = corr.test(abun)
r1 = res$p 
r1[lower.tri(r1,diag = T)] = 0
r1 = r1 + t(r1)
res0 = r1 %>%
    as_tibble() %>%
    mutate_all(function(x)  case_when(
        (x>=0.01 & x<0.05)~'*', 
        (x>=0.001 & x<0.01)~'**',
        x<0.001~'***',
        x>0.05 ~' '))

n_CM = 7
res1 = pheatmap(res$r,clustering_distance_rows = 'correlation',
                clustering_distance_cols = 'correlation',
                clustering_method = 'ward.D2',
                color = colorRampPalette(rev(RColorBrewer::brewer.pal(11,'RdBu')))(100),
                breaks=seq(-1,1,length.out=101),silent = F,border_color = NA,
                # display_numbers =res0,
                # display_numbers = T,
                number_color = 'black',cutree_cols = n_CM,
                cutree_rows =n_CM )
res1

 

modules = data.frame('sub'=names(cutree(res1$tree_row,k=n_CM)),
                     'module'= cutree(res1$tree_row,k=n_CM)) 
modules = modules[res1$tree_row$labels[res1$tree_row$order],]
modules = modules %>%mutate(module=paste0('M',match(modules$module,unique(modules$module))))

rownames(modules)<-NULL;annotation_col = modules %>%
  distinct(sub,module) %>%
  column_to_rownames(var='sub') %>%
  mutate(module = factor(module))
ann_colors = list(
  module = c(M1="#FF6666",M2="#FFB266",M3="#FFCC99",M4="#66B266",M5="#66CCB2",M6="#6699CC",M7="#B266CC"))

pheatmap(res$r,clustering_distance_rows = 'correlation',annotation_col = annotation_col,annotation_colors = ann_colors,
         clustering_distance_cols = 'correlation',clustering_method = 'ward.D2',show_colnames = F,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(11,'RdBu')))(100),
         breaks=seq(-1,1,length.out=101),silent = F,border_color = NA,
         # display_numbers =res0,
         number_color = 'black',cutree_cols = n_CM,

# Figure 2b -----
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

# Figure 2c -----
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/liuliqiu/permutation_All.RData")
all_long <- tidyr::gather(all_sum_mat, key = "sample", value = "closeness", -"cell")


all_long$sample <- lapply(strsplit(all_long$sample, '[.]'),function(x){res<-x[1];return(res);})%>%unlist()
mycell<-"CD4_C7_OX40;CD8_C6_CD39;Mac_C2_SPP1;Endo_C3_RGCC;FB_C3_COL1A1"
all_long$cell_type = all_long$cell
all_long$cell_type[which(all_long$cell!=mycell)]<-"other"
all_long$cell_type[which(all_long$cell==mycell)]<-"mimer"

all = rbind(all_long,all_long_other)
all$cell_type  = factor(all$cell_type,
                             levels = c("mimer",'other','CD4','CD8','Mac','Endo','FB'))
mycolor<-c("#1B9E77","grey","#1e90ff","#20b2aa",'#ba6338',"#ce3d32","#f0e685")
my_comparison<-list(c('mimer','other'))


all_cell_percentage<-as.data.frame(all_cell_percentage)
rownames(all_cell_percentage)<-lapply(strsplit(rownames(all_cell_percentage), '[.]'),function(x){res<-x[1];return(res);})%>%unlist()
ggplot(as.data.frame(all_cell_percentage), aes(x = mimer_percentage))+ geom_density(color = "black", fill = "gray")
drop<-rownames(all_cell_percentage)[all_cell_percentage$mimer_percentage<0.05]
all_sub<-all[all$sample%in%drop==F,]


p<-ggplot(all_sub, aes(x=cell_type, y=closeness, fill = cell_type)) + 
  # geom_violin(color="white")+ 
  # ggrastr::geom_point_rast(colour = "grey", shape=16, alpha = 0.7, stroke = 0, size=3, position = position_jitter(0.12),
  # raster.dpi = getOption("ggrastr.default.dpi", 300),) +
  geom_boxplot(width=0.5, outlier.shape = NA, alpha = 0.8)+ 
  scale_fill_manual(values = mycolor)+ 
  # geom_jitter(shape=16, position = position_jitter(0.1), alpha = 0.6, stroke = 0, size=3, outlier.shape = NA) +
  stat_compare_means(paired = FALSE, comparisons = my_comparison) + 
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(size = 14)) + 
  theme_classic() + 
  NoLegend()

p+coord_fixed(ratio = 0.5)

# Figure 2d -----
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


# Figure 2e -----
load("cpdb_result.Rdata")
ggplot(cpdb,aes(x = LR,y = gene_pair, fill = mean,size = logPval)) +
  geom_point(shape=21,color='black') +
  scale_fill_gradientn(colors = color) +
  scale_y_discrete(limits=rev(levels(plot_df$subC))) +
  scale_size_continuous(range = c(1,8),name = 'Proportion') +
  facet_grid(~type,space = 'free_x',scale = 'free_x',labeller = label_wrap_gen(width=5)) +
  theme_bw() + xlab('') + ylab('') +
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_text(size = 10))

# Figure 2f -----
load('sp_correct.RData')
df<-read.table("sp_table.txt", sep='\t',header = T)
gene<-c("CXCL13","ACP5","LAG3","PHLDA1","HAVCR2","RGS2","FOXP3","TIGIT","BATF","TNFRSF18","TNFRSF4","TNFRSF9","COL1A1","COL3A1","COL1A2","SPARC","FN1","POSTN","PLVAP","COL4A1","COL4A2","HSPG2","VWF","IGFBP7","SPP1","APOC1","MMP12","MMP9","FBP1","APOE")
Idents(sp)<-sp$development;
p<-subset(sp,development=='unknown',invert=T)
Idents(p)<-p$development
p<-AddModuleScore(p,features = list(gene),name = 'MIMER')
length(colnames(p@meta.data))
colnames(p@meta.data)[84]<-'MIMER'
my_levels<-c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')
Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor','Hyp'),c('Hyp','MiD'),c('MiD','MoD'),c('MoD','SD&CA'),c('SD&CA','ICA'),c('ICA','MCA'))

VlnPlot(p,features = c('MIMER'),pt.size = 0,assay = 'spatial',cols = mycolor)+
  geom_signif(comparisons = my_comparison,y_position = 1.4,textsize=4,step_increase=0)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('MIMER expression')

# Figure 2g -----
Idents(sp)<-sp$development;
p<-subset(sp,development=='unknown',invert=T)
Idents(p)<-p$development
p<-AddModuleScore(p,features = list(gene),name = 'MIMER')
length(colnames(p@meta.data))
colnames(p@meta.data)[84]<-'MIMER'
my_levels<-c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')
Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor','Hyp'),c('Hyp','MiD'),c('MiD','MoD'),c('MoD','SD&CA'),c('SD&CA','ICA'),c('ICA','MCA'))

VlnPlot(p,features = c('MIMER'),pt.size = 0,assay = 'spatial',cols = mycolor)+
  geom_signif(comparisons = my_comparison,y_position = 1.4,textsize=4,step_increase=0)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('MIMER expression')

# Figure 2h -----
VlnPlot(p,features = c('MIMER'),pt.size = 0,assay = 'spatial',cols = mycolor)+
  geom_signif(comparisons = my_comparison,y_position = 1.4,textsize=4,step_increase=0)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('EMT score')

# Figure 2j & 2k-----
total_le = rbind(hyp_LE,SDCA_LE)
total_me = rbind(hyp_ME,SDCA_ME)
total_me$distance = gsub("far","Distal",
                         gsub("near","Proximal",
                              total_me$distance))

ggplot(total_me,aes(x=development,y=MIMER,
                          color=distance,
                          fill=distance))+
  geom_violin()+
  stat_compare_means(mapping = aes(condition='distance'),method = "wilcox.test")+NoLegend()+theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=16,color = "black"),
        axis.title.y  = element_text(size=18,color = "black"),
        legend.text  = element_text(size = 12),
        legend.title =  element_text(size = 14))+ggtitle("")+
  ylab('ME\nMIMER expression')

total_major$distance = gsub("far","Distal",
                            gsub("near","Proximal",
                                 total_major$distance))
ggplot(total_major,aes(x=development,y=MIMER,
                             color=distance,
                             fill=distance))+
  geom_violin()+
  stat_compare_means(mapping = aes(condition='distance'),method = "wilcox.test")+NoLegend()+theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=16,color = "black"),
        axis.title.y  = element_text(size=18,color = "black"),
        legend.text  = element_text(size = 12),
        legend.title =  element_text(size = 14))+ggtitle("")

# Figure 2j & 2k-----
load("emt_mimmer_score.Rdata")
ggplot(score, aes(x=EMT, y=MIMER)) + 
  geom_point(aes(color=area),size=4)+ 
  labs(x="OX40+ CD4 in CD4 (%)", y="PDCD1+ CD8 in \n CD8 (%)")+
  stat_cor(method = "pearson",size=6,vjust=0.5,digits = 2)+
  stat_smooth(method='lm', 
              color="black", se=FALSE)+labs(title='')+guides(size=FALSE)+
  theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  theme_classic()+
  theme(axis.text = element_text(size=12,colour = "black"),
        axis.title = element_text(size=15))+ylim(0,50)

# Extended Data Fig. 4:
# Extended Data Fig. 4a  -----

load("scRNA.RData")
#Endothelium
my<-subset(sc,L1_C=='Endothelium');Idents(my)<-my$L4_C
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('Endo_C4_CCL21','Endo_C3_RGCC','Endo_C2_FBLN5','Endo_C1_ACKR1'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="",cellwidth = 60,cellheight = 60)
#Fibroblast
my<-subset(sc,L1_C=='Fibroblast');Idents(my)<-my$L4_C;my<-subset(my,ident='Smooth muscle cell',invert=T)
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('FB_C6_ACTA2','FB_C5_PDGFRB','FB_C4_APOE','FB_C3_COL1A1',
                   'FB_C2_IGF1','FB_C1_CFD'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")
#Bcell
Idents(sc)<-sc$L1_C;my<-subset(sc,ident=c('B','Plasma'));Idents(my)<-my$L4_C
getPalette = colorRampPalette(brewer.pal(6, "YlOrRd"))
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('Plasma','B_C7_MKI67','B_C6_BCL6','B_C5_ISG15','B_C4_TCL1A',
                   'B_C3_GPR183','B_C2_DUSP4','B_C1_CCR7'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")

#Myeloid
my<-subset(sc,ident=c('Myeloid'));Idents(my)<-my$L4_C
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('Neutrophil','Mono_C2_CD16','Mono_C1_CD14','Mac_C6_MKI67','Mac_C5_CXCL10','Mac_C4_LYVE1','Mac_C3_C1QC',
                   'Mac_C2_SPP1','Mac_C1_NLRP3','DC_C5_IL3RA','DC_C4_LAMP3','DC_C3_CD1A','DC_C2_CD1C','DC_C1_CLEC9A'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")
#CD4T
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD4'));Idents(my)<-my$L4_C
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('CD4_C9_IFNG','CD4_C8_ITGB1','CD4_C7_OX40','CD4_C6_CD25','CD4_C5_RTKN2','CD4_C4_GPR183','CD4_C3_ANXA1',
                   'CD4_C2_TCF7','CD4_C13_ISG15','CD4_C12_NKG7','CD4_C11_IL17A','CD4_C10_CXCR5','CD4_C1_CCR7'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")

#CD8T
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD8'));Idents(my)<-my$L4_C
Idents(my)<-my$tissue
mysub <- as.matrix(table(Idents(my),my@meta.data$L4_C))
roe.sub<-t(ROIE(mysub))
roe.sub<-as.data.frame(roe.sub)
roe.sub<-roe.sub[c('CD8_C9_KLRC2','CD8_C8_KIR','CD8_C7_CX3CR1','CD8_C6_CD39','CD8_C5_GZMK','CD8_C4_ANXA1','CD8_C3_TCF7','CD8_C2_IL7R',
                   'CD8_C12_CD161','CD8_C11_ISG15','CD8_C10_CD16','CD8_C1_CCR7'),]
pheatmap(roe.sub, cluster_rows = F, cluster_cols = F,display_numbers = T,
         color = getPalette(10),na_col="white",angle_col = "0",
         fontsize = 15,number_color = "black",main="")



# Extended Data Fig. 4b  -----
library("FactoMineR")
library("factoextra")
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(cluster)
library(ggbreak)
set.seed(123)
produc_pca <- PCA(scale(rall_m), ncp = 3, graph = FALSE)
produc_hcpc <- HCPC(produc_pca, graph = FALSE,method = "ward",min = 7)
fviz_cluster(produc_hcpc,
             repel = TRUE, 
             show.clust.cent = TRUE, labelsize=7,
             palette = "jco",  ggtheme = theme_minimal(),
             main = "")+ 
  theme_classic()+
  scale_x_break(c(4,6.5),
                space = 0.2,
                scales = 0.2)

# Extended Data Fig. 4c  -----
library(ggrepel)
library(ggbreak)
roeMYE$major = "Myeloids"
roeT$major = "CD4T"
roeENDO$major = "Endothelium"
roecd8$major = 'CD8T'
roe_fb$major ="Fibroblast"
roeB$major ="Bcell"

rall <- rbind(roeMYE,roeCD4T,roeENDO,roecd8,roe_fb,roeB)
rall$cluster = rownames(rall)
rall$Nor =rall$Nor+0.001
rall$rank = (rall$Tu*rall$pLN)/(rall$Nor+rall$Adj+rall$nLN)
rall$label = ifelse(rall$cluster %in% c("CD4_C7_OX40","CD8_C6_CD39","Mac_C2_SPP1","FB_C3_COL1A1","Endo_C3_RGCC"),"MIMER","other")
rall %>% 
  ggplot(.,aes(x=rank,y=rank))+
  geom_point()+
  geom_text_repel(aes(color=label),label=rall$cluster)

rall %>% arrange(rank)

plotRSS_oneSet(rall, setName = "rank")

library(ggplot2)
library(ggrepel)

rss = rall
setName = "rank"
rss = rss %>% as.matrix()
rssThisType <- sort(rss[, "rank"], decreasing = TRUE)
thisRss <- data.frame(regulon = names(rssThisType), rank = seq_along(rssThisType), 
                      rss = rssThisType)
thisRss$cluster = rownames(thisRss)
thisRss$regulon = ifelse(thisRss$cluster %in% c("CD4_C7_OX40","CD8_C6_CD39","Mac_C2_SPP1",
                                                "CD4_C13_ISG15","CD8_C11_ISG15",
                                                "FB_C3_COL1A1","Endo_C3_RGCC"),thisRss$cluster,NA)
thisRss$rss = as.numeric(thisRss$rss)
thisRss2 = thisRss %>% left_join(rall,by = "cluster")
ggplot(thisRss2, aes(x = log2(rank.x), y = log2(rss))) +
  geom_point(aes(color=major),shape=4) + 
  ggtitle(setName) +
  geom_text_repel(aes(label = regulon),
                  na.rm = TRUE
  ) + theme_classic()
rownames(thisRss2) <- thisRss2$cluster
thisRss2 %>% dplyr::filter(label %in% "MIMER") ->thisRss3

thisRss2$major <- factor(thisRss2$major,
                         levels = unique(thisRss2$major))

thisRss2$major <-  gsub("CD8T","CD8 T",
                        gsub("CD4T","CD4 T",thisRss2$major ))
ggplot() +
  geom_point(data = thisRss2,
             aes(x = rank.x,
                 y = log2(rss),
                 fill=major,
                 color=major
                 ),
             size=5) +
  geom_text_repel(data=thisRss2,
                    aes(x = rank.x,
                        y = log2(rss),
                        label = regulon),
                  size=4,
                    na.rm = TRUE,hjust=-1,vjust=1,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="grey20",segment.size=0.2,segment.alpha=0.8,nudge_y=1)+
  theme_classic()+
  scale_color_manual(values = c("#5050fe","#1e90ff","#20b2aa","#ce3d32","#f0e685","#ba6338"))+
  scale_fill_manual(values = c("#5050fe","#1e90ff","#20b2aa","#ce3d32","#f0e685","#ba6338"))+
  theme(axis.text = element_text(size = 20,colour = "black"),
        axis.title = element_text(size = 20,colour = "black"),
        legend.position = c(0.3,0.3),
        legend.text = element_text(size=12))+ylab("Tumor enrichment score")+xlab("Rank")+coord_fixed(ratio = 7)

# Extended Data Fig. 4d  -----
library(forcats)
library(ggplot2)
load('gs_5_cluster.Rdata')
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


# Extended Data Fig. 4e  -----
# Immunosuppression
is = c('CD274','PDCD1LG2','ICOSLG','CD276','VTCN1','VSIR','IDO1','TGFB1','IL10','CCL5','CCL17','CCL22','CXCL8','CCL16','CCL18','IL17A','CCL24','IL1R1','IL4','IL13','ARG1','CD40LG','TPH1','XBP1','SOCS1')
my = AddModuleScore(my,features=list(is),name='is')
# Angiogenesis
ag = c('VEGFA','CXCL8','CXCL1','CXCL2','TEK','ANGPT2','CXCL12','FLT1','CXCR4','CTSB','CTSS','SEMA4D','PROK2','CXCR1','S1PR1','PDGFA','FGF2','MMP12','VEGFC','DNTTIP2','PGF','HIF1A','NFE2L2','CCND2','CCNE1','CD44','E2F3','EDN1','EZH2','FGF18','FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2','NOTCH1','PTK2','SPP1','STC1','TNFAIP6','TYMP','VAV2')
my = AddModuleScore(my,features=list(ag),name='ag')
# Invasion and metastasis
im = c('MMP7','MMP2','MMP3','MMP9','CCL2','CCL3','CCL22','TGFB1','EGF','CCR2','FLT1','CTSS','CTSB','CXCR3','WNT7B','ALOX5','CDH1','CCL18','CCL24','S1PR1','CXCL16','MAPK7','CSF1','SPARC','TLR4','VCAM1','CCL20')
my = AddModuleScore(my,features=list(im),name='im')
# Tumorigenesis
tg = c('CXCL8','CXCL12','CSF1R','CCR2','EGF','IL6','IL1B','IL23A','MSR1','CCL5','PTGS2','MFGE8')
my = AddModuleScore(my,features=list(tg),name='tg')
# Tumor suppression
ts = c('IL1B','IL6','IL12A','IL23A','TNF','CXCL9','CXCL10','TLR2','TLR4','CXCL11','IFNG','CD40','FCGR2A','ITGAX','IFNGR1','HLA-DPB1','HLA-DPA1','HLA-DRA','HLA-DRB1','HLA-DQA1','CD74','HLA-DRB5','IRF5')
my = AddModuleScore(my,features=list(ts),name='ts')
# Phagocytosis
pc = c('MRC1','CD163','MERTK','C1QB','FCRLA','CD5L','CD81','GPNMB','CD36')
my = AddModuleScore(my,features=list(pc),name='pc')

my@meta.data = dplyr::rename(my@meta.data,ImmunoSuppression=is1, Angiogenesis=ag1, Metastasis=im1, Tumorigenesis=tg1, TumorSuppression=ts1, Phagocytosis=pc1)

fea <- c('ImmunoSuppression','Angiogenesis','Metastasis','Tumorigenesis','TumorSuppression','Phagocytosis')
myd = as.data.frame(scale(my@meta.data[,fea])); myd = cbind(myd,my@meta.data[,c('tissue','subC')]); myd = myd[,-7]
myd = dplyr::group_by(myd, subC) %>% summarise(across(all_of(fea), mean))
myd = as.data.frame(myd); myd = tibble::column_to_rownames(myd,var='subC'); myd = t(as.matrix(myd))

pheatmap(myd, cellwidth=30, cellheight=30, cluster_cols=TRUE, clustering_method='ward.D', cluster_rows=T, treeheight_row=0,
         kmeans_k=NA, border_color="white", scale="none", drop_levels=T, show_rownames=TRUE, show_colnames=TRUE, name='control',
         color=colorRampPalette(c("navy","white","firebrick3"))(50), fontsize=12, fontsize_row=15, legend=T, fontsize_col=15)

# Extended Data Fig. 4f  -----
# pipeline
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

# Extended Data Fig. 4g  -----
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


# Extended Data Fig. 4h  -----

my<-subset(sc,L3_C=='Endothelium');Idents(my)<-my$L4_C
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
mycolor<-cell_color[cell_color$ct%in%levels(Idents(my)),]$color
VlnPlot(my, feature='ITGA1',pt.size=0, cols = mycolor)+NoLegend()
VlnPlot(my, feature='ITGA2',pt.size=0, cols = mycolor)+NoLegend()
VlnPlot(my, feature='ITGA5',pt.size=0, cols = mycolor)+NoLegend()
VlnPlot(my, feature='KDR',pt.size=0, cols = mycolor)+NoLegend()

# Extended Data Fig. 4i  -----
gene1<-c('CXCL13','ACP5','LAG3','PHLDA1','HAVCR2','RGS2','PLPP1','RHOB','SNX9','CCL5','CD8A','CD3D') # CD39+CD8
gene2<-c('FOXP3','TIGIT','BATF','TNFRSF18','TNFRSF4','TNFRSF9','IL32','CD4','IL10','IL2RA') #Treg
gene3<-c('SPP1','APOC1','MMP12','MMP9','FBP1','APOE','CTSB','CD68','CCL3','TYROBP') # SPP1+MAC
gene4<-c('PLVAP','COL4A1','COL4A2','HSPG2','VWF','IGFBP7','PECAM1','SERPINE1','SPARC','INSR') # Endo
gene5<-c('COL1A1','COL3A1','COL1A2','SPARC','FN1','POSTN','CST1','MMP11','CTHRC1','COL6A3') # CAF
gene6<-c('RNASE1','CCL18','C1QA','C1QB','C1QC','SELENOP','F13A1','PLTP','LGMN','LYVE1','CD68') # LYVE1+Mac
Idents(sp)<-sp$development;p<-subset(sp,ident=c('Nor','SD&CA','ICA'))
my_levels<-c('Nor','SD&CA','ICA')
Idents(p)<-factor(Idents(p),levels = my_levels)
p<-AddModuleScore(p,features = list(gene3),name = 'SPP1.TAM');p<-AddModuleScore(p,features = list(gene1),name = 'CD39.CD8')
p<-AddModuleScore(p,features = list(gene2),name = 'OX40.Treg');p<-AddModuleScore(p,features = list(gene4),name = 'PLVAP.Endo')
p<-AddModuleScore(p,features = list(gene5),name = 'POSTN.CAF');p<-AddModuleScore(p,features = list(gene6),name = 'LYVE1.Mac')
fea<-c('SPP1.TAM1','PDCD1.CD81','OX40.Treg1','PLVAP.Endo1','POSTN.CAF1','LYVE1.Mac1')
myd=as.data.frame(scale(p@meta.data[,fea]))
myd=cbind(myd,p@meta.data[,c('file','development','tissue','metastasis')])
colnames(myd)<-c('SPP1+TAM_score','PDCD1+CD8_score','OX40+Treg_score','PLVAP+Endo_score','POSTN+CAF_score','LYVE1+Mac_score','file','Development','tissue','metastasis')
s1<-myd[myd$Development=='Nor',];s1<-aggregate(s1[,1:6],by=list(Group=s1$file),mean);colnames(s1)[1]<-'file';s1$Development<-'Nor'
s2<-myd[myd$Development=='SD&CA',];s2<-aggregate(s2[,1:6],by=list(Group=s2$file),mean);colnames(s2)[1]<-'file';s2$Development<-'SD&CA'
s3<-myd[myd$Development=='ICA',];s3<-aggregate(s3[,1:6],by=list(Group=s3$file),mean);colnames(s3)[1]<-'file';s3$Development<-'ICA'
new<-rbind(s1,s2);new<-rbind(new,s3)
mycolor<-c("#1B9E77","#66A61E","#E6AB02")
my_comparison<-list(c('Nor','SD&CA'),c('Nor','ICA'));ggboxplot(new,x="Development",y="SPP1+TAM_score",color = "Development",palette = mycolor,add = "jitter",add.params = list(size=0.4))+
  stat_compare_means(comparisons = my_comparison)+NoLegend()

# Extended Data Fig. 4j  -----
p<-subset(sp,file=='IT2');Idents(p)<-p$development;p<-subset(p,ident='unknown',invert=T);p$development<-Idents(p)
p<-AddModuleScore(p,features = list(gene3),name = 'SPP1.TAM');p<-AddModuleScore(p,features = list(gene1),name = 'PDCD1.CD8')
p<-AddModuleScore(p,features = list(gene2),name = 'OX40.Treg');p<-AddModuleScore(p,features = list(gene4),name = 'RGCC.Endo')
p<-AddModuleScore(p,features = list(gene5),name = 'POSTN.CAF');p<-AddModuleScore(p,features = list(gene6),name = 'LYVE1.Mac')
colnames(p@meta.data)[84:89]<-c('SPP1+TAM_score','PDCD1+CD8_score','OX40+Treg_score','RGCC+Endo_score','POSTN+CAF_score','LYVE1+Mac_score')
SpatialFeaturePlot(p,features = c('SPP1+TAM_score'),alpha = c(1,1),ncol=1,images = df[df$file=='IT2',]$image)
SpatialFeaturePlot(p,features = c('CD39+CD8_score'),alpha = c(1,1),ncol=1,images = df[df$file=='IT2',]$image)
SpatialFeaturePlot(p,features = c('OX40+Treg_score'),alpha = c(1,1),ncol=1,images = df[df$file=='IT2',]$image)
SpatialFeaturePlot(p,features = c('PLVAP+Endo_score'),alpha = c(1,1),ncol=1,images = df[df$file=='IT2',]$image)
SpatialFeaturePlot(p,features = c('POSTN+CAF_score'),alpha = c(1,1),ncol=1,images = df[df$file=='IT2',]$image)

# Extended Data Fig. 4K  -----
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
    # low = "#4575b4",
    # high = "#f0663f",
    # mid = "#f1f9d8"
  )+   guides(size= guide_legend(ncol = 2))+
  theme_graph() + expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))

#Extended Data Fig. 5 a b c d -----
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

#Extended Data Fig. 6
library(pheatmap)

total = read.csv("nor.csv",header = T, row.names = 1)

annotation = read.csv("annotation.csv",header = T, row.names = 1)

bk <- c(seq(-3,-0.01,by=0.01),seq(0,3,by=0.01))

pheatmap(total, cluster_cols = FALSE, cluster_rows = FALSE,annotation_row = annotation_row,
color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
colorRampPalette(colors = c("white","red"))(length(bk)/2)),breaks=bk,show_rownames = TRUE,show_colnames = FALSE)

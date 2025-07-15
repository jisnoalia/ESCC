library(psych)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(Seurat)
library(dplyr)
require("DESeq2")
# Extended Data Fig 7 -----
# Extended Data Fig 7a -----
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
         number_color = 'black',cutree_cols = n_CM)

# Extended Data Fig 7b -----
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
# Extended Data Fig 7c -----
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
#Extended Data Fig 9------
#Extended Data Fig 9a------
my = subset(sc, L4_C == 'FB_C3_COL1A1' & tissue == 'Tu');Idents(my)<-my$metastasis
deg<-FindMarkers(my,assay='RNA',logfc.threshold = 0,ident.1 = 'Y',only.pos = F,max.cells.per.ident = 1000)
logFC_cutoff <- 0.3
deg$change = as.factor(
  ifelse(deg$p_val < 0.05 & abs(deg$avg_log2FC) > logFC_cutoff,
         ifelse(deg$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
dat<-deg;dat$gene_name<-rownames(dat)
m2d = function(x){
  mean(abs(x))+2*sd(abs(x))
}  
gene<-c('COL14A1','CST1','CHI3L1','COL11A1','COL1A1','CCL11','WNT5A','MMP11','MMP1','IGF2','LRRC15','MMP3')
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
        legend.text = element_text(colour = "black",size=14),legend.title = element_text(colour = "black",size=14))
#Extended Data Fig 9b------
my = subset(sc, L4_C == 'FB_C3_COL1A1' & tissue == 'Tu');Idents(my)<-my$metastasis
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
VlnPlot(my, feature='COL1A1',pt.size=0, cols = c("#377EB8","#E41A1C"),y.max = 8)+NoLegend()+stat_compare_means(paired = FALSE)

# spatial results
load('sp_correct.RData')
df<-read.table("sp_table.txt", sep='\t',header = T)
Idents(sp)<-sp$development;p<-subset(sp,ident=c('SD&CA'))
Idents(p)<-p$metastasis
mycolor<-c("#377EB8","#E41A1C")
my_comparison<-list(c('Y','N'))
VlnPlot(p,features = c('COL1A1'),pt.size = 0,assay = 'spatial',cols = mycolor,y.max = 4.5)+
  stat_compare_means(comparisons = my_comparison)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.text.x= element_text(color = 'black',angle = 0))
Idents(sp)<-sp$development;p<-subset(sp,ident=c('ICA'))
Idents(p)<-p$metastasis
my_levels<-c('N','Y');Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#377EB8","#E41A1C")
my_comparison<-list(c('Y','N'))
VlnPlot(p,features = c('COL1A1'),pt.size = 0,assay = 'spatial',cols = mycolor,y.max = 6)+
  stat_compare_means(comparisons = my_comparison)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.text.x= element_text(color = 'black',angle = 0))
#Extended Data Fig 9e------
library(ggsignif)
library(ggpubr)
library(ggprism)
mrtm <- data.frame(num = c(27,15,25,9),
                   abund = c("Low","Low","High","High"),
                   metastasis = c("No","Yes","No","Yes"))

MRTM = reshape2::melt(mrtm)

ggplot(data = MRTM,aes(x=abund,y=value,fill=metastasis))+
  geom_bar(stat = "identity",
           position = "stack",
           width =0.4)+
  scale_fill_manual(values = c("#377eb8","brown"))+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  geom_text(aes(label=value),vjust=3,size=6,color="black")+
  annotate(geom = "text",
           cex=6,
           x=1.5, y=70, 
           label=paste0("p = 0.032"),
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=12))+labs(x="Abundance of LYVE1+MRTM",y="Patient numbers",fill = "Lymph node metastasis")
#Extended Data Fig 9f------
my = subset(sc, L3_C == 'Macrophage')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
mycolor<-cell_color[cell_color$ct%in%levels(Idents(my)),]$color
VlnPlot(my, feature='CXCR4',pt.size=0, cols = mycolor)+NoLegend()
VlnPlot(my, feature='IL1B',pt.size=0, cols = mycolor)+NoLegend()
#Extended Data Fig 9g------
my = subset(sc, L4_C == 'Mac_C4_LYVE1' & tissue == 'Tu');Idents(my)<-my$metastasis
deg<-FindMarkers(my,assay='RNA',logfc.threshold = 0,ident.1 = 'Y',only.pos = F)
logFC_cutoff <- 0.3
deg$change = as.factor(
  ifelse(deg$p_val < 0.05 & abs(deg$avg_log2FC) > logFC_cutoff,
         ifelse(deg$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
dat<-deg;dat$gene_name<-rownames(dat)
m2d = function(x){
  mean(abs(x))+2*sd(abs(x))
}  
gene<-c('MERTK','CD163')
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
        legend.text = element_text(colour = "black",size=14),legend.title = element_text(colour = "black",size=14))
#Extended Data Fig 9h------
load('scRNA.RData');
primary_all = subset(sc,tissue=='Adj')
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
df2 = df2 %>% filter(tissue=="Adj")
rownames(df2) = df2$cellID
my$tissue = stringr::str_split(my$orig.ident,"_",simplify = T)[,2]
sub$cellID = rownames(sub@meta.data)

sub$cellID = rownames(sub@meta.data)

all_cell_tumor = df2

df2 = df2 %>% filter(newC!="Smooth muscle cell")

dat<-as.data.frame(table(df2$patient,df2$newC));
colnames(dat)<-c('sample_id','sub','sub_num');
dat<-deseq2_norm_single(dat)
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

n_CM = 4
res1 = pheatmap(res$r,clustering_distance_rows = 'correlation',
                clustering_distance_cols = 'correlation',
                clustering_method = 'ward.D2',
                color = colorRampPalette(rev(RColorBrewer::brewer.pal(11,'RdBu')))(100),
                breaks=seq(-1,1,length.out=101),silent = F,border_color = NA,
                # display_numbers =res0,
                # display_numbers = T,
                number_color = 'black',cutree_cols = n_CM,
                cutree_rows =n_CM )

#Extended Data Fig 9 c&i------
metastasis = development %>%
  distinct(file,.keep_all = T) %>%
  dplyr::select(file,metastasis)

dev_cc_ratio %>% as.data.frame()%>%
  mutate(region_id =  paste0(file ,"_", development)) %>% 
  arrange(development) %>% distinct(region_id,.keep_all = T) %>% 
  dplyr::select(development,region_id,file) -> group_info 

group_info = group_info %>%
  left_join(metastasis,by = "file")
head(group_info)

head(dev_cc_ratio)
ICA_ratio = dev_cc_ratio %>%
  filter(development %in% 'ICA')
head(ICA_ratio)
metastasis = development %>% distinct(file,.keep_all = T) %>% dplyr::select(file,metastasis)
head(metastasis)
df_filtered_iqr_grouped <- ICA_ratio %>% left_join(metastasis,by = "file") %>%
  group_by(metastasis,cc_15) %>% # 按照 'group' 列进行分组
  mutate(
    Q1 = quantile(ratio, 0.25),
    Q3 = quantile(ratio, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(ratio >= lower_bound & ratio <= upper_bound) %>% # 过滤掉超出边界的行
  select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound) # 删除临时计算的列

df_filtered_iqr_grouped$metastasis = factor(df_filtered_iqr_grouped$metastasis,
                                            levels = c("N","Y"))
df_filtered_iqr_grouped %>%
  ggplot(.,aes(color=metastasis,y=ratio,x= cc_15))+
  # geom_jitter_rast()+
  geom_boxplot(width=0.4,
               position = position_dodge(width=0.6))+
  stat_compare_means(
    mapping = aes(group=metastasis),label = "p.format",
                     show.legend = F)+
  scale_color_manual(values =c("Y" = "#e4172c",'N' = '#0c79be') )+
  # scale_fill_manual(values =c("Y" = "#e4172c",'N' = '#0c79be') )+
  theme_bw()+
  theme(axis.text = element_text(size = 12,color = "black"))+
  labs(y = 'The ratio of CC in ICA',x = 'Cell Compostion') -> p 

p
ggsave(filename = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/liuliqiu/RCTD/ISCHIA/ISCHIA/cc_ratio/cc_ratio_in_development_ICA_boxplot.pdf",width = 12,height = 4)
#Extended Data Fig 9j------
my = subset(sc, L3_C == 'Cancer');my$group<-paste(my$metastasis,my$tissue,sep = '_')
Idents(my)<-my$group;my<-subset(my,ident=c('N_Tu','Y_pLN','Y_Tu'))
mhc<-c('HLA-A','HLA-B','HLA-C','TAP1','TAP2','B2M')
my<-AddModuleScore(my,features = list(mhc),name = 'MHCI_score');colnames(my@meta.data)[35]<-'MHCI_score'
my_levels<-c('N_Tu','Y_pLN','Y_Tu');Idents(my)<-factor(Idents(my),levels = my_levels)
my_comparison<-list(c('N_Tu','Y_pLN'),c('Y_pLN','Y_Tu'),c('N_Tu','Y_Tu'))
VlnPlot(my, feature='MHCI_score',pt.size=0, cols = c("#a44a3f","#a44a3f","#a44a3f"),y.max = 6)+NoLegend()+
  stat_compare_means(comparisons = my_comparison)+theme(axis.text.x = element_text(angle = 0,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
                                                        panel.grid.major=element_blank(),panel.grid.minor = element_blank())+ ylab('MHCI Score')+ggtitle('')

#Extended Data Fig 9k------
my = subset(sc, L3_C == 'Cancer' & tissue == 'Tu');Idents(my)<-my$metastasis
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
VlnPlot(my, feature=c('IL1RN','HLA-DRA'),pt.size=0, cols = c("#377EB8","#E41A1C"))+NoLegend()&theme(axis.text.x = element_text(angle = 0,hjust = 1))
#Extended Data Fig 9l-----
my = subset(sc, L3_C == 'Cancer' & tissue == 'Tu');Idents(my)<-my$metastasis
gene<-c('VIM','SPRR1A','SPRR1B','PI3')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 6,scale.min = 0)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())
gene<-c('TUBA1A','TUBB','TUBA1B','TOP2A','STMN1','HIF1A','CAV1','KRT4','FGFBP1','IGFBP3','IGFBP5')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 6,scale.min = 0)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())
#Extended Data Fig 9m------
my = subset(sc, L3_C == 'Cancer' & tissue == 'Tu');Idents(my)<-my$metastasis
gene<-c('IGF2R','SLC1A5')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
p1<-DotPlot(my,features = gene,dot.scale = 6,scale.min = 0)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))+labs(title = "Cancer")

my = subset(sc, L4_C == 'FB_C3_COL1A1' & tissue == 'Tu');Idents(my)<-my$metastasis
gene<-c('IGF2','CCL11','DPP4','LGALS9')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
p2<-DotPlot(my,features = gene,dot.scale = 6,scale.min = 0)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))+labs(title = "POSTN+ CAF")
     

#Figure 4 -----
#Figure 4 a -----
library(fgsea)
library(msigdbr)
load("area_hallmark_score.Rdata")
pheatmap(area_hallmark_score, cellwidth=18, cellheight=15, cluster_rows=T, cluster_cols=T, clustering_method='ward.D2', 
        kmeans_k=NA, border_color="white", scale="row", drop_levels=T, show_rownames=TRUE, 
        show_colnames=TRUE,color=rev(getPalette(10)), fontsize_col=13, fontsize_row=11, legend=T)


#Figure 4 b -----
mfuzz.ggplot2 = function(eset, cl, min.mem = 0, time.points = NULL, time.labels = NULL,
                         highlight_genes=highlight_genes,text.size=0.2,membership.color = rev(getPalette(12)),
                         xlab = "Time", ylab = "Expression changes", centre = FALSE, line.size=0.2,
                         selected.clusters = NULL) {
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  
  if(is.null(selected.clusters)) {
    selected.clusters <- unique(clusterindex) # 若没指定，则包括所有cluster
  }
  
  df <- NULL
  for (j in selected.clusters) {
    tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
    tmpmem <- memship[clusterindex == j, j]
    
    if(length(tmpmem) > 0) {
      for (i in 1:dim(tmp)[1]) {
        ggdata <- data.frame("Gene" = row.names(tmp)[i],
                             "Time" = if(is.null(time.points)) as.numeric(as.character(1:ncol(tmp))) else time.points,
                             "Expression" = tmp[i,],
                             "Membership" = tmpmem[i],
                             "Cluster" = j)
        df <- rbind(df, ggdata)
      }
    }
  }
  highlight_genes <- highlight_genes
  df_highlight <- df %>% 
    filter(Time == max(Time) & Gene %in% highlight_genes)
  df2 = df[order(df$Membership),]
  # df$Cluster = factor(df$Cluster,
  #                     levels = selected.clusters)
  df$Gene = factor(df$Gene,
                   levels = unique(df2$Gene))
  p <-ggplot(df, aes(x = Time, y = Expression, group = Gene, color = Membership)) +
    geom_line(size=line.size) + 
    geom_line(data = df %>% filter(Gene %in% highlight_genes), 
              aes(x = Time, y = Expression, group = Gene), 
              color = "black", 
              size = line.size)+
    
    geom_text_repel(data = df_highlight, 
              aes(label = Gene), color="black",size=text.size,
              point.padding = NA,
              size = 3,
              nudge_x = 100000,
              segment.size = 0.5,
              segment.color = "black",
              direction = "y",
              hjust = 1
              # hjust = -1
              ) +
    facet_wrap(~Cluster,ncol = 1) +
    scale_x_continuous(name = xlab, breaks = if(is.null(time.points)) as.numeric(as.character(1:ncol(eset))) else time.points, labels = if(is.null(time.labels)) colnames(eset) else time.labels) +
    ylab(ylab) +
    scale_color_gradientn(colours = membership.color)+theme_bw()
  return(p)
}

mytheme <- theme(plot.title=element_text(
                                         size="20"),
                 axis.title=element_blank(),
                 axis.text=element_text(size = 20),
                 axis.line = element_line(linetype = 1,colour = "black"),

                 panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                 panel.border = element_blank(),panel.background = element_blank(),
                 legend.position=NULL,
                 strip.text = element_text(colour = "black",
                                           hjust = 0.5, size = 16),
                 strip.background = element_rect(color="black",fill = "white"),
                 plot.margin = margin(1, 2, 1, 2, "cm")
                 )

load("mfuzz_dat.Rdata" )
dat <- Mfuzz::filter.NA(dat)
dat <- Mfuzz::fill.NA(dat, mode = 'mean')
dat <- Mfuzz::filter.std(dat, min.std = 0)
dat <- standardise(dat)
n <- 12
m <- mestimate(dat)
set.seed(123)
cl <- mfuzz(dat, c = n, m = m)
p1 <- mfuzz.ggplot2(dat,cl,time.labels = colnames(cl$cluster),selected.clusters = c(11),
              line.size = 0.5,text.size = 4,
              highlight_genes = c("MMP11","PDPN",'COL5A2',
                                  "EGFR",'MAP3K5','HDAC9',
                                  'VIM','GPC1','COPA',
                                  'GPT2','ALAD','ENDOD1',
                                  "PSMD4","BANF1","BUD31"))+mytheme+NoLegend()
p1+coord_cartesian(xlim = c(0.5, 8), ylim = c(-2, 2.5), expand = c(0, 0))

ggsave(p1,filename = "mfuzz_11.pdf")
p2 <- mfuzz.ggplot2(dat,cl,time.labels = colnames(cl$cluster),selected.clusters = c(3),
                    line.size = 0.5,text.size = 4,
                    highlight_genes = c("MMP11","PDPN",'COL5A2',
                                        "EGFR",'MAP3K5','HDAC9',
                                        'VIM','GPC1','COPA',
                                        'GPT2','ALAD','ENDOD1',
                                        "PSMD4","BANF1","BUD31"))+mytheme+NoLegend()
p2+coord_cartesian(xlim = c(0.5, 8), ylim = c(-2.2, 2.5), expand = c(0, 0))

#Figure 4 c -----
load("ds.RData")
### SCENT 
library(SCENT)
neu<-dsc@assays$spatial@data
gene <- bitr(rownames(neu),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
index<-match(gene[,1],rownames(neu))
rownames(neu)[index]<-gene[,2]
neu<-neu[index,]
ccat.v <- CompCCAT(exp = neu, ppiA = net13Jun12.m);
dsc$development<-factor(dsc$development,levels=unique(dsc$development))
dsc$development<-factor(dsc$development,levels=c("Nor","Hyp","MiD","MoD","SD&CA","ICA","MCA"))
ccat.v_mtx<-data.frame("CCAT"=ccat.v,"development"=dsc$development)
my_comparisons<-list(c("Nor","Hyp"),c("Hyp","MiD"),c("MiD","MoD"),c("MoD","SD&CA"),c("SD&CA","ICA"),c("ICA","MCA"))
ggboxplot(ccat.v_mtx,x="development",y="CCAT",fill = "development",palette = c("#219773","#CF5D15","#6F6BAA","#D62F81","#64A032","#DDA715","#A17424"))+
   theme(
    axis.title = element_text(size=18),
    axis.text.x = element_text(size=18,angle = 0,hjust=0.5,colour = "black"),
    axis.text.y = element_text(size=18,colour = "black"))+
  geom_signif(comparisons = my_comparisons, y_position  = 0.3, textsize=6, step_increase=0)+NoLegend()

#Figure 4 d -----
load('sp_correct.RData')
df<-read.table("sp_table.txt", sep='\t',header = T)
my<-subset(sp,ME=='unknown',invert=T)
library(fgsea)
library(msigdbr)
msigdbr_collections()
human<-msigdbr(species = "human")
a<-human[grep('LEUKOCYTE_DEGRANULATION',human$gs_name),];table(a$gs_name)
gene1<-human[human$gs_name=='REACTOME_SIGNALING_BY_VEGF',]$human_gene_symbol #GOBP_CELLULAR_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES
gene2<-human[human$gs_name=='REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION',]$human_gene_symbol
gene3<-human[human$gs_name=='BIOCARTA_IL6_PATHWAY',]$human_gene_symbol
gene4<-human[human$gs_name=='BIOCARTA_IL10_PATHWAY',]$human_gene_symbol
gene5<-human[human$gs_name=='GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE',]$human_gene_symbol
gene6<-human[human$gs_name=='GOBP_VASCULATURE_DEVELOPMENT',]$human_gene_symbol
gene7<-human[human$gs_name=='WP_TGFBETA_SIGNALING_PATHWAY',]$human_gene_symbol
gene8<-human[human$gs_name=='BIOCARTA_EGF_PATHWAY',]$human_gene_symbol
gene9<-human[human$gs_name=='GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY',]$human_gene_symbol
gene10<-human[human$gs_name=='GOMF_FIBRONECTIN_BINDING',]$human_gene_symbol
gene11<-human[human$gs_name=='REACTOME_TNF_SIGNALING',]$human_gene_symbol
gene12<-human[human$gs_name=='REACTOME_INTERFERON_SIGNALING',]$human_gene_symbol
gene13<-human[human$gs_name=='BIOCARTA_IL2_PATHWAY',]$human_gene_symbol
gene14<-human[human$gs_name=='GOBP_PHAGOCYTOSIS',]$human_gene_symbol
gene15<-human[human$gs_name=='GOBP_T_CELL_MEDIATED_CYTOTOXICITY',]$human_gene_symbol
gene16<-human[human$gs_name=='GOBP_CELL_KILLING',]$human_gene_symbol
df<-data.frame(pathway=c('VEGF_Signaling','ECM_Organization','IL6_Pathway','IL10_Pathway','Inflammatory_Response','Vasculature_Development','TGFB_Pathway','EGF_Pathway','Myeloid_Immunity','Fibronectin_Binding','TNF_Pathway',
                         'Interferon_Signaling','IL2_Pathway','Phagocytosis','T_Cytotoxicity','Cell_Killing'))
my<-AddModuleScore(my,features = list(gene1),name = df[1,'pathway']);my<-AddModuleScore(my,features = list(gene2),name = df[2,'pathway']);my<-AddModuleScore(my,features = list(gene3),name = df[3,'pathway'])
my<-AddModuleScore(my,features = list(gene4),name = df[4,'pathway']);my<-AddModuleScore(my,features = list(gene5),name = df[5,'pathway']);my<-AddModuleScore(my,features = list(gene6),name = df[6,'pathway'])
my<-AddModuleScore(my,features = list(gene7),name = df[7,'pathway']);my<-AddModuleScore(my,features = list(gene8),name = df[8,'pathway']);my<-AddModuleScore(my,features = list(gene9),name = df[9,'pathway'])
my<-AddModuleScore(my,features = list(gene10),name = df[10,'pathway']);my<-AddModuleScore(my,features = list(gene11),name = df[11,'pathway']);my<-AddModuleScore(my,features = list(gene12),name = df[12,'pathway'])
my<-AddModuleScore(my,features = list(gene13),name = df[13,'pathway']);my<-AddModuleScore(my,features = list(gene14),name = df[14,'pathway']);my<-AddModuleScore(my,features = list(gene15),name = df[15,'pathway'])
my<-AddModuleScore(my,features = list(gene16),name = df[16,'pathway'])
colnames(my@meta.data)[84:99]<-df$pathway
myd=as.data.frame(my@meta.data[,df$pathway])
myd=cbind(myd,my@meta.data[,c('ME','file')])
myd<-aggregate(myd[,1:16],by=list(ME=myd$ME),mean);rownames(myd)<-myd$ME;myd<-myd[,-1];myd<-t(myd)
myd<-myd[,c('Nor-ME','Hyp-ME','MiD-ME','MoD-ME','SD&CA-ME','ICA-ME','MCA-ME')]
pheatmap(myd,color=rev(getPalette(10)),border_color = "grey60", fontsize = 12, angle_col = "90",clustering_method = "ward.D2",cluster_cols = F,cluster_rows = T,scale = "row")


#Figure 4 e -----
load('sp_correct.RData')
p<-subset(sp,TLS=='unknown',invert=T);Idents(p)<-p$TLS;DefaultAssay(p)<-'spatial'
gene<-c('CCL5','CCL11','CCL19','CXCL9','CXCL10','CXCL14')
VlnPlot(p, features = gene,stacked=T,pt.size=0,cols = mycolor)+theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5,color = 'black'),axis.text.y = element_text(color = 'black'))+xlab('')


#Figure 4 f ------
load('sp_correct.RData')
p<-subset(sp,TLS=='unknown',invert=T);Idents(p)<-p$TLS;DefaultAssay(p)<-'spatial'
my_levels<-c('Hyp-TLS','MiD-TLS','SD&CA-TLS','ICA-TLS');Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#D95F02","#666666","#66A61E","#E6AB02")
gene<-c('CD19','MS4A1','CD79A')
VlnPlot(p, features = gene,stacked=T,pt.size=0,cols = mycolor)+theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5,color = 'black'),axis.text.y = element_text(color = 'black'))+xlab('')
gene<-c('IGHG2','IGLC2','IGLC7')
VlnPlot(p, features = gene,stacked=T,pt.size=0,cols = mycolor)+theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5,color = 'black'),axis.text.y = element_text(color = 'black'))+xlab('')
gene<-c('CR2','CD22')
VlnPlot(p, features = gene,stacked=T,pt.size=0,cols = mycolor)+theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5,color = 'black'),axis.text.y = element_text(color = 'black'))+xlab('')
gene<-c('CXCL13','CXCR5')
VlnPlot(p, features = gene,stacked=T,pt.size=0,cols = mycolor)+theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5,color = 'black'),axis.text.y = element_text(color = 'black'))+xlab('')
gene<-c('CD3D','CD4','CD8A','PDCD1')
VlnPlot(p, features = gene,stacked=T,pt.size=0,cols = mycolor)+theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5,color = 'black'),axis.text.y = element_text(color = 'black'))+xlab('')
gene<-c('PECAM1','RGCC','ACKR1')
VlnPlot(p, features = gene,stacked=T,pt.size=0,cols = mycolor)+theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5,color = 'black'),axis.text.y = element_text(color = 'black'))+xlab('')
gene<-c('COL1A1','POSTN','CFD')
VlnPlot(p, features = gene,stacked=T,pt.size=0,cols = mycolor)+theme(axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5,color = 'black'),axis.text.y = element_text(color = 'black'))+xlab('')

#Figure 4 g ------
library(reshape)
dec<-read.table("/spatialDWLS/results/TLS_nomean.txt", sep='\t',header = T)
sub<-dec[dec$Cluster%in%c('B','DC','Plasma','T','Endothelium','Fibroblast'),]
mycolor<-c("#D95F02","#666666","#66A61E","#E6AB02")
sub$Area<-factor(sub$Area,levels = c('Hyp-TLS','MiD-TLS','SD&CA-TLS','ICA-TLS'))
ggboxplot(sub,x="Cluster",y="Ratio",color = "Area",palette = mycolor,add = "jitter",add.params = list(size=0.1))
sub<-dec[dec$Cluster%in%c('Cancer','Epithelium','Mast','MoMac','Neutrophil','NK'),] # figs8
mycolor<-c("#D95F02","#666666","#66A61E","#E6AB02")
sub$Area<-factor(sub$Area,levels = c('Hyp-TLS','MiD-TLS','SD&CA-TLS','ICA-TLS'))
ggboxplot(sub,x="Cluster",y="Ratio",color = "Area",palette = mycolor,add = "jitter",add.params = list(size=0.1))

#Figure 4 h ------
## TLS gene scores in each sample
gene1<-c('CD19','MS4A1','CD79A');gene2<-c('IGHG2','IGLC2','IGLC7');gene3<-c('CR2','CD22');gene4<-c('CXCL13','CXCR5')
gene5<-c('CD3D','CD4','CD8A','PDCD1');gene6<-c('PECAM1','RGCC','ACKR1');gene7<-c('COL1A1','POSTN','CFD')
# Hyp-TLS
p<-subset(sp,file=='QN')
p$use <- ifelse(p$TLS %in% "Hyp-TLS",1,0)
col = c("white", "#D95F02")
SpatialFeaturePlot(p,features = c('use'),alpha = c(0,1),ncol=1,images = df[df$file=='QN',]$image )+
  scale_fill_gradientn(colours = col)+NoLegend()
DefaultAssay(p)<-'spatial'
p<-subset(p,TLS=='Hyp-TLS')
p<-AddModuleScore(p,features = list(gene1),name = "B");p<-AddModuleScore(p,features = list(gene2),name = "Plasma");p<-AddModuleScore(p,features = list(gene3),name = "fDC")
p<-AddModuleScore(p,features = list(gene4),name = "Tfh");p<-AddModuleScore(p,features = list(gene5),name = "T")
p<-AddModuleScore(p,features = list(gene6),name = "Endo");p<-AddModuleScore(p,features = list(gene7),name = "CAF")
colnames(p@meta.data)[85:91]<-c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score')
SpatialFeaturePlot(p,features = c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score'),alpha = c(1,0.6),ncol=4,images = df[df$file=='QN',]$image,pt.size.factor = 2.5)&scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(-1, 1.5))
# MiD-TLS
p<-subset(sp,file=='D_JT2')
p$use <- ifelse(p$TLS %in% "MiD-TLS",1,0)
col = c("white", "#666666")
SpatialFeaturePlot(p,features = c('use'),alpha = c(0,1),ncol=1,images = df[df$file=='D_JT2',]$image )+
  scale_fill_gradientn(colours = col)+NoLegend()
DefaultAssay(p)<-'spatial'
p<-subset(p,TLS=='MiD-TLS')
p<-AddModuleScore(p,features = list(gene1),name = "B");p<-AddModuleScore(p,features = list(gene2),name = "Plasma");p<-AddModuleScore(p,features = list(gene3),name = "fDC")
p<-AddModuleScore(p,features = list(gene4),name = "Tfh");p<-AddModuleScore(p,features = list(gene5),name = "T")
p<-AddModuleScore(p,features = list(gene6),name = "Endo");p<-AddModuleScore(p,features = list(gene7),name = "CAF")
colnames(p@meta.data)[85:91]<-c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score')
SpatialFeaturePlot(p,features = c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score'),alpha = c(1,0.6),ncol=4,images = df[df$file=='D_JT2',]$image,pt.size.factor = 5)&scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(-1, 1.5))
# SD&CA-TLS
p<-subset(sp,file=='KT2')
p$use <- ifelse(p$TLS %in% "SD&CA-TLS",1,0)
col = c("white", "#66A61E")
SpatialFeaturePlot(p,features = c('use'),alpha = c(0,1),ncol=1,images = df[df$file=='KT2',]$image )+
  scale_fill_gradientn(colours = col)+NoLegend()
DefaultAssay(p)<-'spatial'
p<-subset(p,TLS=='SD&CA-TLS')
p<-AddModuleScore(p,features = list(gene1),name = "B");p<-AddModuleScore(p,features = list(gene2),name = "Plasma");p<-AddModuleScore(p,features = list(gene3),name = "fDC")
p<-AddModuleScore(p,features = list(gene4),name = "Tfh");p<-AddModuleScore(p,features = list(gene5),name = "T")
p<-AddModuleScore(p,features = list(gene6),name = "Endo");p<-AddModuleScore(p,features = list(gene7),name = "CAF")
colnames(p@meta.data)[85:91]<-c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score')
SpatialFeaturePlot(p,features = c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score'),alpha = c(1,0.6),ncol=4,images = df[df$file=='KT2',]$image,pt.size.factor = 2)&scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(-1, 1.5))
# ICA-TLS
p<-subset(sp,file=='KT1_2')
p$use <- ifelse(p$TLS %in% "ICA-TLS",1,0)
col = c("white", "#E6AB02")
SpatialFeaturePlot(p,features = c('use'),alpha = c(0,1),ncol=1,images = df[df$file=='KT1_2',]$image )+
  scale_fill_gradientn(colours = col)+NoLegend()
DefaultAssay(p)<-'spatial'
p<-subset(p,TLS=='ICA-TLS')
p<-AddModuleScore(p,features = list(gene1),name = "B");p<-AddModuleScore(p,features = list(gene2),name = "Plasma");p<-AddModuleScore(p,features = list(gene3),name = "fDC")
p<-AddModuleScore(p,features = list(gene4),name = "Tfh");p<-AddModuleScore(p,features = list(gene5),name = "T")
p<-AddModuleScore(p,features = list(gene6),name = "Endo");p<-AddModuleScore(p,features = list(gene7),name = "CAF")
colnames(p@meta.data)[85:91]<-c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score')
SpatialFeaturePlot(p,features = c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score'),alpha = c(1,0.6),ncol=4,images = df[df$file=='KT1_2',]$image,pt.size.factor = 2.5)&scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(-1, 1.5))


#Extended Data Fig. 9 ------

#Extended Data Fig. 9a ------


#Extended Data Fig. 9b ------
load('sp_correct.RData')
df<-read.table("sp_table.txt", sep='\t',header = T)
p<-subset(sp,ME=='unknown',invert=T);Idents(p)<-p$ME
hall<-as.data.frame(table(human_H$gs_name))
ex<-c("INFLAMMATORY_RESPONSE",
      "INTERFERON_GAMMA_RESPONSE","EPITHELIAL_MESENCHYMAL_TRANSITION","E2F_TARGETS","G2M_CHECKPOINT");ex<-paste("HALLMARK",ex,sep = "_")
hall<-hall[hall$Var1%in%ex,]
p<-subset(sp,development=='unknown',invert=T);Idents(p)<-p$development;DefaultAssay(p)<-'spatial'
for (i in hall$Var1){
  gene<-human_H[human_H$gs_name==i,]$gene_symbol
  p<-AddModuleScore(p,features = list(gene),name = i)
} 
colnames(p@meta.data)[84:91]<-as.character(hall$Var1)
myd=as.data.frame(scale(p@meta.data[,as.character(hall$Var1)]));colnames(myd)<-gsub('HALLMARK_','',colnames(myd))
myd=cbind(myd,p@meta.data[,c('development','file')])
s1<-myd[myd$development=='Nor',];s1<-aggregate(s1[,1:8],by=list(Group=s1$file),mean);colnames(s1)[1]<-'file';s1$development<-'Nor'
s2<-myd[myd$development=='Hyp',];s2<-aggregate(s2[,1:8],by=list(Group=s2$file),mean);colnames(s2)[1]<-'file';s2$development<-'Hyp'
s3<-myd[myd$development=='MiD',];s3<-aggregate(s3[,1:8],by=list(Group=s3$file),mean);colnames(s3)[1]<-'file';s3$development<-'MiD'
s4<-myd[myd$development=='MoD',];s4<-aggregate(s4[,1:8],by=list(Group=s4$file),mean);colnames(s4)[1]<-'file';s4$development<-'MoD'
s5<-myd[myd$development=='SD&CA',];s5<-aggregate(s5[,1:8],by=list(Group=s5$file),mean);colnames(s5)[1]<-'file';s5$development<-'SD&CA'
s6<-myd[myd$development=='ICA',];s6<-aggregate(s6[,1:8],by=list(Group=s6$file),mean);colnames(s6)[1]<-'file';s6$development<-'ICA'
s7<-myd[myd$development=='MCA',];s7<-aggregate(s7[,1:8],by=list(Group=s7$file),mean);colnames(s7)[1]<-'file';s7$development<-'MCA'
new<-rbind(s1,s2);new<-rbind(new,s3);new<-rbind(new,s4);new<-rbind(new,s5);new<-rbind(new,s6);new<-rbind(new,s7)
mycolor<-c('#1B9E77','#D95F02','#666666','#E7298A','#66A61E','#E6AB02','#A6761D')
my_comparison<-list(c('MoD','SD&CA'),c('SD&CA','ICA'))
ggboxplot(new,x="development",y="EPITHELIAL_MESENCHYMAL_TRANSITION",color = "development",palette = mycolor,add = "jitter",add.params = list(size=0.1),width = 0.7)+stat_compare_means(comparisons = my_comparison)+NoLegend()
my_comparison<-list(c('MoD','SD&CA'))
ggboxplot(new,x="development",y="E2F_TARGETS",color = "development",palette = mycolor,add = "jitter",add.params = list(size=0.1),width = 0.7)+stat_compare_means(comparisons = my_comparison)+NoLegend()
my_comparison<-list(c('MoD','SD&CA'))
ggboxplot(new,x="development",y="G2M_CHECKPOINT",color = "development",palette = mycolor,add = "jitter",add.params = list(size=0.1),width = 0.7)+stat_compare_means(comparisons = my_comparison)+NoLegend()
ggboxplot(new,x="development",y="INFLAMMATORY_RESPONSE",color = "development",palette = mycolor,add = "jitter",add.params = list(size=0.1),width = 0.7)+NoLegend()

#Extended Data Fig. 9c ------
Idents(my)<-my$ME;my$ME<-Idents(my)
features<-c('STAT5A','MAPK3','RACK1','TRAF1','CD6','GPSM3','CTSH','ULBP2','CYBA','CORO1A',
            'IL6','JAK3','IGFBP3','FSTL3','MMP1','TGFB1I1','COL1A1','COL1A2','SULF1','THBS2',
            'ISG15','IFI6','EGFR','CALM1','AKT1','CALML5','S100P')
exp <- my@assays$SCT@data[features,] %>% as.data.frame() %>% t() %>% as.data.frame()
exp$subC <- my@meta.data[rownames(exp),'ME']

trunc_z <- function(x) MinMax((x-mean(x))/sd(x),-2,2)

exp %>% reshape2::melt() %>% group_by(subC,variable) %>% #summarise(N = n()) %>% as.data.frame()
  summarise(mean = mean(value),frac = sum(value > 0)/length(value)) %>% 
  group_by(variable) %>% mutate(mean_z = trunc_z(mean)) -> plot_df
celltype <- c('IL2 Pathway','TNF Pathway','Inflammatory Response','T_Cytotoxicity Cell_Killing','Phagocytosis','IL6/IL10 Pathway',
              'Fibronectin Binding','TGFB Pathway','ECM Organization','Vasculature Development','Interferon Signaling','EGF Pathway','VEGF Signaling','Myeloid Immunity')
ctidx <- c(2,2,2,2,2,2,2,2,2,2,2,1,2,2)
gType <- rep(celltype,ctidx)
names(gType) <- features

plot_df$gType <- gType[plot_df$variable]
plot_df$gType <- factor(plot_df$gType,
                        levels = celltype)
color<-rev(getPalette(10))
ggplot(plot_df,aes(x = variable,y = subC, fill = mean_z,size = frac)) +
  geom_point(shape=21,color='black') +
  scale_fill_gradientn(colors = color) +
  scale_y_discrete(limits=rev(levels(plot_df$subC))) +
  scale_size_continuous(range = c(1,8),name = 'Proportion') +
  facet_grid(~gType,space = 'free_x',scale = 'free_x',labeller = label_wrap_gen(width=5)) +
  theme_bw() + xlab('') + ylab('') +
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_text(size = 10))



#Extended Data Fig. 10-----
#Extended Data Fig. 10 b -----
## summarize DEG in TLS + GSEA
p<-subset(sp,TLS=='unknown',invert=T);Idents(p)<-p$TLS;DefaultAssay(p)<-'spatial';p$TLS<-Idents(p)
my_levels<-c('Hyp-TLS','MiD-TLS','SD&CA-TLS','ICA-TLS');Idents(p)<-factor(Idents(p),levels = my_levels)
Idents(p)<-p$TLS
hyp<-FindMarkers(p,assay='spatial',ident.1 = 'Hyp-TLS',only.pos = T,test.use = 'roc');hyp$Gene<-rownames(hyp);hyp<-arrange(hyp,desc(myAUC))
mid<-FindMarkers(p,assay='spatial',ident.1 = 'MiD-TLS',only.pos = T,test.use = 'roc');mid$Gene<-rownames(mid);mid<-arrange(mid,desc(myAUC))
sa<-FindMarkers(p,assay='spatial',ident.1 = 'SD&CA-TLS',only.pos = T,test.use = 'roc');sa$Gene<-rownames(sa);sa<-arrange(sa,desc(myAUC))
ica<-FindMarkers(p,assay='spatial',ident.1 = 'ICA-TLS',only.pos = T,test.use = 'roc');ica$Gene<-rownames(ica);ica<-arrange(ica,desc(myAUC))
gene<-c('CXCL10','CXCL9','CXCR4','CXCL11','PTPRC','KRT4','RPS14','RPL28','SPRR3','TYROBP',
        'PFN1','FTH1','TSTD1','TXN2','SEM1','RHOA','IL4I1','MYO9B','IDO1','S100A7',
        'IGHG4','IGHG1','MIF','IFI27','IFI6','BIRC5','CCND1','CSF3','CXCL1','TNFRSF18',
        'COL1A1','COL3A1','SPARC','ITGB1','SERPINH1','SDC1','THBS1','FN1','HSPG2','DCN')

dotplot <- function(mca,features,celltype,ctidx,color){
  exp <- mca@assays$spatial@data[features,] %>% as.data.frame() %>% t() %>% as.data.frame()
  exp$subC <- mca@meta.data[rownames(exp),'TLS']
  
  trunc_z <- function(x) MinMax((x-mean(x))/sd(x),-2,2)
  
  exp %>% reshape2::melt() %>% group_by(subC,variable) %>% #summarise(N = n()) %>% as.data.frame()
    summarise(mean = mean(value),frac = sum(value > 0)/length(value)) %>% 
    group_by(variable) %>% mutate(mean_z = trunc_z(mean)) -> plot_df
  
  gType <- rep(celltype,ctidx)
  names(gType) <- features
  
  plot_df$gType <- gType[plot_df$variable]
  plot_df$gType <- factor(plot_df$gType,
                          levels = celltype)
  
  ggplot(plot_df,aes(x = variable,y = subC, fill = mean_z,size = frac)) +
    geom_point(shape=21,color='black') +
    scale_fill_gradientn(colors = color) +
    scale_y_discrete(limits=rev(levels(plot_df$subC))) +
    scale_size_continuous(range = c(1,8),name = 'Proportion',limits = c(0,1)) +
    facet_grid(~gType,space = 'free_x',scale = 'free_x') +
    theme_bw() + xlab('') + ylab('') +
    theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank())
}
T_plotG<-gene
area <- c('Hyp-TLS','MiD-TLS','SD&CA-TLS','ICA-TLS')
ctidx <- c(10,10,10,10)
dotplot(p,T_plotG,area,ctidx,color = rev(getPalette(10)))


#Extended Data Fig. 10 c -----
# GSEA
rank<-hyp[,c('Gene','myAUC')]
geneList<-rank$myAUC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
human_H<-msigdbr(species = "human", category = "H")
human_C2<-msigdbr(species = "human", category = "C2")
human_C2<-human_C2[human_C2$gs_subcat%in%c('CP:BIOCARTA','CP:KEGG','CP:REACTOME'),]
human<-rbind(human_C2,human_H);gmt<-human %>% dplyr::select(gs_name,gene_symbol)
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000,
                pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
library(GseaVis)
library(clusterProfiler)
geneSetID = c('REACTOME_CELLULAR_RESPONSES_TO_EXTERNAL_STIMULI','REACTOME_DEVELOPMENTAL_BIOLOGY','REACTOME_METABOLISM_OF_RNA')
gseaNb(object = gseaRes,geneSetID = geneSetID,addPoint = F, addPval = T,
       pCol = 'black', pHjust = 0,subPlot = 2,pvalX = 1,pvalY = 1.1,curveCol = brewer.pal(n = 3, name = "Set1"))

rank<-mid[,c('Gene','myAUC')]
geneList<-rank$myAUC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000,
                pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('KEGG_RIBOSOME','REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES','REACTOME_SIGNALING_BY_ROBO_RECEPTORS')
gseaNb(object = gseaRes,geneSetID = geneSetID,addPoint = F, addPval = T,
       pCol = 'black', pHjust = 0,subPlot = 2,pvalX = 1,pvalY = 1.1,curveCol = brewer.pal(n = 3, name = "Set1"))

rank<-sa[,c('Gene','myAUC')]
geneList<-rank$myAUC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000,
                pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING','REACTOME_SIGNALING_BY_INTERLEUKINS','REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM')
gseaNb(object = gseaRes,geneSetID = geneSetID,addPoint = F, addPval = T,
       pCol = 'black', pHjust = 0,subPlot = 2,pvalX = 1,pvalY = 1.1,curveCol = brewer.pal(n = 3, name = "Set1"))

rank<-ica[,c('Gene','myAUC')]
geneList<-rank$myAUC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000,
                pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('REACTOME_ECM_PROTEOGLYCANS','REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES','KEGG_ECM_RECEPTOR_INTERACTION')
gseaNb(object = gseaRes,geneSetID = geneSetID,addPoint = F, addPval = T,
       pCol = 'black', pHjust = 0,subPlot = 2,pvalX = 1,pvalY = 1.1,curveCol = brewer.pal(n = 3, name = "Set1"))




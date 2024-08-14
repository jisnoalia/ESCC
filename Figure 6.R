#Figure 6 ------
#Figure 6 a ------
my = subset(sc, L3_C == 'Fibroblast' & tissue == 'Tu');Idents(my)<-my$metastasis
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Metastasis','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
mycolor<-cell_color[cell_color$ct%in%levels(t$Cell_Subcluster),]$color
ggplot(t, aes(Metastasis,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Metastasis')+ylab('Proportion')+scale_fill_manual(values =mycolor)

#Figure 6 b ------
load("scRNA.RData")
tf<-subset(sc,tissue=='Tu');Idents(tf)<-tf$L3_C
my<-subset(tf,ident=c('Fibroblast'))
y<-subset(my,metastasis=='Y');n<-subset(my,metastasis=='N')
t<-as.data.frame(table(y@meta.data$patient,y@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(y@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Metastasis<-'Y'
com<-t
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Metastasis<-'N'
com<-rbind(com,t)
my_comparisons<-list(c('Y','N')) 
sub<-com[com$Subcluster=='FB_C3_COL1A1',]
colnames(sub)[5]<-"%POSTN+CAF in Fibroblast";sub$Metastasis<-factor(sub$Metastasis,levels = c('N','Y'))
ggboxplot(sub,x="Metastasis",y="%POSTN+CAF in Fibroblast",color = "Metastasis",palette = c("#377EB8","#E41A1C"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

load('sp_correct.RData')
gene5<-c('COL1A1','COL3A1','COL1A2','SPARC','FN1','POSTN','CST1','MMP11','CTHRC1','COL6A3')
p<-subset(sp,development=='ICA');Idents(p)<-p$development
p<-AddModuleScore(p,features = list(gene5),name = 'POSTN.CAF')
fea<-c('POSTN.CAF1')
myd=as.data.frame(scale(p@meta.data[,fea]))
myd=cbind(myd,p@meta.data[,c('file','development','tissue','metastasis')])
colnames(myd)<-c('POSTN+CAF_score','file','development','tissue','metastasis')
s1<-myd[myd$development=='ICA'&myd$metastasis=='Y',];s1<-aggregate(s1[,1],by=list(Group=s1$file),mean);colnames(s1)[1]<-'file';s1$development<-'ICA';s1$metastasis<-'Y'
s2<-myd[myd$development=='ICA'&myd$metastasis=='N',];s2<-aggregate(s2[,1],by=list(Group=s2$file),mean);colnames(s2)[1]<-'file';s2$development<-'ICA';s2$metastasis<-'N'
new<-rbind(s1,s2);colnames(new)[2]<-'POSTN+CAF_score'
new$metastasis<-factor(new$metastasis,levels = c('N','Y'))
ggboxplot(new,x="metastasis",y="POSTN+CAF_score",color = "metastasis",palette = c("#377EB8","#E41A1C"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

#Figure 6 c ------
Idents(sc)<-sc$L3_C;my = subset(sc, ident=c('Macrophage','Monocyte'));my<-subset(my,tissue == 'Tu');Idents(my)<-my$metastasis
t<-as.matrix(table(Idents(my),my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(my)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Metastasis','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
mycolor<-cell_color[cell_color$ct%in%levels(t$Cell_Subcluster),]$color
ggplot(t, aes(Metastasis,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Metastasis')+ylab('Proportion')+scale_fill_manual(values =mycolor)

#Figure 6 d ------
load("scRNA.RData")
tf<-subset(sc,tissue=='Tu');Idents(tf)<-tf$L3_C
my<-subset(tf,ident=c('Macrophage','Monocyte'))
y<-subset(my,metastasis=='Y');n<-subset(my,metastasis=='N')
t<-as.data.frame(table(y@meta.data$patient,y@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(y@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Metastasis<-'Y'
com<-t
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Metastasis<-'N'
com<-rbind(com,t)
my_comparisons<-list(c('Y','N')) 
sub<-com[com$Subcluster=='Mac_C4_LYVE1',]
colnames(sub)[5]<-"%LYVE1+MRTM in Macrophage+Monocyte";sub$Metastasis<-factor(sub$Metastasis,levels = c('N','Y'))
ggboxplot(sub,x="Metastasis",y="%LYVE1+MRTM in Macrophage+Monocyte",color = "Metastasis",palette = c("#377EB8","#E41A1C"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

load('sp_correct.RData')
gene6<-c('RNASE1','CCL18','C1QA','C1QB','C1QC','SELENOP','F13A1','PLTP','LGMN','LYVE1','CD68') # LYVE1+Mac
p<-subset(sp,development=='SD&CA');Idents(p)<-p$development
p<-AddModuleScore(p,features = list(gene6),name = 'LYVE1.Mac')
fea<-c('LYVE1.Mac1')
myd=as.data.frame(scale(p@meta.data[,fea]))
myd=cbind(myd,p@meta.data[,c('file','development','tissue','metastasis')])
colnames(myd)<-c('LYVE1+Mac_score','file','development','tissue','metastasis')
s1<-myd[myd$development=='SD&CA'&myd$metastasis=='Y',];s1<-aggregate(s1[,1],by=list(Group=s1$file),mean);colnames(s1)[1]<-'file';s1$development<-'SD&CA';s1$metastasis<-'Y'
s2<-myd[myd$development=='SD&CA'&myd$metastasis=='N',];s2<-aggregate(s2[,1],by=list(Group=s2$file),mean);colnames(s2)[1]<-'file';s2$development<-'SD&CA';s2$metastasis<-'N'
new<-rbind(s1,s2);colnames(new)[2]<-'LYVE1+MRTM_score'
new$metastasis<-factor(new$metastasis,levels = c('N','Y'))
ggboxplot(new,x="metastasis",y="LYVE1+MRTM_score",color = "metastasis",palette = c("#377EB8","#E41A1C"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

#Figure 6 e ------
my = subset(sc, L3_C == 'Macrophage')
gene<-c('MRC1','LYVE1','IL10','HLA-DRB5','HLA-DPA1')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 8)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())

#Figure 6 f ------
library(plotrix)
my = subset(sc, tissue == 'Tu' & L3_C == 'Macrophage')
t<-as.matrix(table(my@meta.data$tissue,my@meta.data$L4_C))
t<-as.data.frame(t)
total<-as.data.frame(table(my@meta.data$tissue))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
c4<-t[t$Cell_Subcluster=='Mac_C4_LYVE1',]
mysub <- as.matrix(table(my@meta.data$tissue,my@meta.data$L4_C))
roe.sub<-ROIE(mysub);roe.sub<-as.data.frame(roe.sub);roe.sub$Tissue<-rownames(roe.sub);roe.sub<-roe.sub[,c('Mac_C4_LYVE1','Tissue')];colnames(roe.sub)[1]<-'roe'
c4<-merge(c4,roe.sub,all=TRUE)
c4$Tissue<-factor(c4$Tissue,levels = c('PBMC','nLN','pLN','Nor','Adj','Tu'))
c4<-c4[order(c4$Tissue),]
c4[["Tissue"]] <- factor(c4[["Tissue"]], levels = as.character(c4[["Tissue"]]))
twoord.plot(lx=1:6,ly=c4$ratio,rx=1:6,ry=c4$roe,type=c('bar','line'),
            lcol = '#F06061', rcol = 'steelblue', 
            ylab = 'Proportion of cell', 
            rylab = 'RO/E', xtickpos=1:6, xticklab = c4$Tissue,lwd=2)

#Figure 6 h ------
library(ggsignif)
library(ggpubr)
library(ggprism)
mrtm <- data.frame(num = c(30,32,9,27),
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

#Figure 6 i ------
my<-subset(sc,tissue=='Tu');Idents(my)<-my$L3_C
fb<-subset(my,L3_C=='Fibroblast')
sub<-subset(fb,L4_C=='FB_C3_COL1A1')
p<-as.data.frame(table(sub$patient));colnames(p)<-c('Patient','sub')
pt<-as.data.frame(table(fb$patient));colnames(pt)<-c('Patient','total')
p<-merge(p,pt,by='Patient');p$ratio<-p$sub/p$total*100;p<-p[!is.nan(p$ratio),]
p<-p[,c(1,4)];colnames(p)<-c('Patient','FBRatio');com<-p
mac<-subset(my,ident=c('Macrophage','Monocyte'))
sub<-subset(mac,L4_C=='Mac_C4_LYVE1')
p<-as.data.frame(table(sub$patient));colnames(p)<-c('Patient','sub')
pt<-as.data.frame(table(mac$patient));colnames(pt)<-c('Patient','total')
p<-merge(p,pt,by='Patient');p$ratio<-p$sub/p$total*100;p<-p[!is.nan(p$ratio),]
p<-p[,c(1,4)];colnames(p)<-c('Patient','MacRatio')
com<-merge(com,p,by='Patient')
dat.lm <- lm(FBRatio ~ MacRatio, com)
formula <- sprintf("italic(y) == %.2f %+.2f * italic(x)", round(coef(dat.lm)[1],2), round(coef(dat.lm)[2],2))
r2 <- sprintf("italic(R^2) == %.3f", summary(dat.lm)$r.squared)
pvalue<-sprintf("italic(P-value) == %.6f", Regressionp(dat.lm))
labels <- data.frame(formula=formula, r2=r2, pvalue=pvalue,stringsAsFactors = FALSE)
ggplot(com, aes(x=FBRatio, y=MacRatio)) + geom_point(size=1)+ stat_smooth(method='lm', color="red", fill="#69b3a2", se=TRUE)+
  labs(x="POSTN+CAF/Fibroblast Ratio", y="LYVE1+MRTM/Macrophage+Monocyte Ratio") + theme_classic()+
  geom_text(data=labels, mapping=aes(x=50,y=45,label=formula), parse = TRUE, inherit.aes = FALSE,size=5) + 
  geom_text(data=labels, mapping=aes(x=50,y=40,label=r2),parse = TRUE,inherit.aes = FALSE, size=5)+
  geom_text(data=labels, mapping=aes(x=50,y=35,label=pvalue),parse = TRUE,inherit.aes = FALSE, size=5)+
  theme(axis.text = element_text(size=12,colour = "black"),axis.title = element_text(size=15))


#Figure 6 k ------
gene5<-c('COL1A1','COL3A1','COL1A2','SPARC','FN1','POSTN','CST1','MMP11','CTHRC1','COL6A3')
gene6<-c('RNASE1','CCL18','C1QA','C1QB','C1QC','SELENOP','F13A1','PLTP','LGMN','LYVE1','CD68')
p<-subset(sp,file=='QP');Idents(p)<-p$development;p<-subset(p,development=='unknown',invert=T);p$development<-Idents(p)
p<-AddModuleScore(p,features = list(gene5),name = 'POSTN.CAF')
p<-AddModuleScore(p,features = list(gene6),name = 'LYVE1.Mac')
colnames(p@meta.data)[84:85]<-c('POSTN+CAF_score','LYVE1+MRTM_score')
SpatialPlot(p,alpha = c(1,0.1),ncol=1,images = df[df$file=='QP',]$image,group.by = 'development',cols = "#66A61E")
SpatialFeaturePlot(p,features = c('POSTN+CAF_score'),alpha = c(1,0.6),ncol=1,images = df[df$file=='QP',]$image)&scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(-0.5, 2.5))
SpatialFeaturePlot(p,features = c('LYVE1+MRTM_score'),alpha = c(1,0.6),ncol=1,images = df[df$file=='QP',]$image)&scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(-0.5, 0.5))

#Figure 6 l ------
my = subset(sc, L3_C == 'Cancer' & tissue == 'Tu');Idents(my)<-my$metastasis
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
gene<-c('HLA-B','HLA-C','CD74','HLA-DRA','MDK','HLA-DRB1','VIM','CD24','KRT13','KRT6C','PI3','SPRR1B','SPRR1A','S100A9',
        'S100A8','SPRR3','SPRR2A','SPRR2D','KRT14','KRT16','SPRR2E','S100A7')
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

#Figure 6 m ------
deg<-FindMarkers(my,assay='RNA',logfc.threshold = 0,ident.1 = 'Y',only.pos = F,max.cells.per.ident = 1000);deg$Gene<-rownames(deg)
rank<-deg[,c('Gene','avg_log2FC')]
geneList<-rank$avg_log2FC
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
geneSetID = c('KEGG_CELL_ADHESION_MOLECULES_CAMS','REACTOME_KERATINIZATION')
gseaNb(object = gseaRes,geneSetID = geneSetID,addPoint = F, addPval = T,
       pCol = 'black', pHjust = 0,subPlot = 2,pvalX = 1,pvalY = 1.1,curveCol = brewer.pal(n = 2, name = "Set1"))


#Figure 6 n ------
load('scRNA.RData'); sc = subset(sc, L2_C=='Cancer'); sc = subset(sc, tissue=='pLN' | tissue=='Tu')
sc$class = paste(sc$tissue, sc$metastasis, sep='_'); Idents(sc) = 'class'
myg = rowSums(as.matrix(sc[['RNA']]@counts) > 0) >= 10; sc = subset(sc, features=rownames(sc[['RNA']]@counts[myg,]))

GS <- escape::getGeneSets(species="Homo sapiens", library="H")
colors = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(20))

res = escape::enrichIt(sc, gene.sets=GS, groups=1000, cores=8)
colnames(res) = str_remove(colnames(res), "HALLMARK_"); sc = AddMetaData(sc, res)


dittoHeatmap(sc, genes=NULL, metas=names(res), annot.by="class", fontsize=7, cluster_cols=TRUE, heatmap.colors=colors(50))

VlnPlot(sc, pt.size=0, features='DNA_REPAIR') + NoLegend()
VlnPlot(sc, pt.size=0, sort='decreasing', split.plot=T, split.by='class', features='DNA_REPAIR')

DotPlot(sc, features=colnames(res), cols=c("navy","firebrick3")) + coord_flip() +
  theme(axis.text.x=element_text(size=14,colour="black",angle=45,hjust=1), axis.title=element_blank(),
        axis.text.y=element_text(size=12,colour="black"), axis.line=element_line(linewidth=0.7),legend.text=element_text(size=12),
        legend.title=element_text(size=12))

ES2 <- sc@meta.data[,c(names(res), 'class')]; colnames(ES2)[ncol(ES2)] <- "cluster"
output <- getSignificance(ES2, group="cluster", fit="linear.model")

myd = output[, names(table(Idents(sc)))]
myp = c("ALLOGRAFT_REJECTION","APICAL_JUNCTION","APICAL_SURFACE","COAGULATION","HEDGEHOG_SIGNALING","HEME_METABOLISM",
        "KRAS_SIGNALING_DN","MYC_TARGETS_V2","MYOGENESIS","PANCREAS_BETA_CELLS","PEROXISOME","SPERMATOGENESIS","UV_RESPONSE_DN",
        "UV_RESPONSE_UP","XENOBIOTIC_METABOLISM","MITOTIC_SPINDLE","BILE_ACID_METABOLISM","IL2_STAT5_SIGNALING",
        "ESTROGEN_RESPONSE_EARLY","CHOLESTEROL_HOMEOSTASIS")
myd = myd[ !(rownames(myd) %in% myp),]
pheatmap(myd, cellwidth=18, cellheight=14, cluster_rows=T, cluster_cols=T, clustering_method='ward.D2', kmeans_k=NA, name='control',
         border_color="white", scale="row", drop_levels=T, show_rownames=TRUE, show_colnames=TRUE, color=colors, fontsize_col=13,
         fontsize_row=11, legend=T)

load('spatial.RData'); 
res = escape::enrichIt(sp, gene.sets=GS, groups=1000, cores=8)
colnames(res) = str_remove(colnames(res), "HALLMARK_"); sc = AddMetaData(sc, res)
myp = c('TGF_BETA_SIGNALING','NOTCH_SIGNALING','DNA_REPAIR','G2M_CHECKPOINT','MYC_TARGETS_V1','P53_PATHWAY','MTORC1_SIGNALING',
        'REACTIVE_OXYGEN_SPECIES_PATHWAY','OXIDATIVE_PHOSPHORYLATION','HYPOXIA','APOPTOSIS','GLYCOLYSIS',
        'UNFOLDED_PROTEIN_RESPONSE')
myd = myd[ rownames(myd) %in% myp,]
pheatmap(myd, cellwidth=18, cellheight=14, cluster_rows=T, cluster_cols=T, clustering_method='ward.D2', kmeans_k=NA, name='control',
         border_color="white", scale="row", drop_levels=T, show_rownames=TRUE, show_colnames=TRUE, color=colors, fontsize_col=13,
         fontsize_row=11, legend=T)


#Figure 6 o ------
s<-read.table("MacC4_CAF_Cancer_YN_summary.txt",header = T,sep = '\t')
target<-s[s$interacting_pair=='IGF2_IGF2R'&s$celltype_pairs=='FB_C3_COL1A1.Cancer_C4'&s$group=='N>Y',];com<-target
target<-s[s$interacting_pair=='LGALS9_SLC1A5'&s$celltype_pairs=='FB_C3_COL1A1.Cancer_C4'&s$group=='Y_Only',];com<-rbind(com,target)
s<-read.table("FBC3_YN_summary.txt",header = T,sep = '\t')
target<-s[s$interacting_pair=='DPP4_CCL11'&s$celltype_pairs=='FB_C3_COL1A1.FB_C3_COL1A1'&s$group=='N>Y',];com<-rbind(com,target)
com$celltype_pairs<-gsub('Cancer_C4','Cancer',com$celltype_pairs)
y<-com[,1:4];colnames(y)<-c('interacting_pair','celltype_pairs','mean','pvalue');y$group<-'Y'
n<-com[,c(1,2,5,6)];colnames(n)<-c('interacting_pair','celltype_pairs','mean','pvalue');n$group<-'N'
com<-rbind(y,n);com[is.na(com$mean),'mean']=0;com[is.na(com$pvalue),'pvalue']=1
com$celltype_pairs<-gsub('FB_C3_COL1A1','POSTN+CAF',com$celltype_pairs);com$celltype_pairs<-gsub('Mac_C4_LYVE1','LYVE1+MRTM',com$celltype_pairs)
com %>% 
  ggplot(aes(celltype_pairs,interacting_pair) )+ 
  geom_point(aes(fill=mean,size=-log10(pvalue+0.0001)),shape=21,color='black') +
  scale_fill_gradientn(colours = rev(getPalette(10)),name = 'Mean') +
  theme_bw() + coord_flip()+
  theme(axis.text.x = element_text(size=14,angle = 60,hjust=1,colour = "black"),panel.background = element_blank(),strip.background = element_rect(fill = c("#377EB8","#E41A1C")),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_text(size=10),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+xlab('Interacting cell type')+ylab('Paired Gene')+facet_wrap(~group)

#Figure 6 p ------
p<-subset(sp,file=='KT1_1')
p$colocalization<-'Other';s<-subset(p,ME=='ICA-ME');s<-subset(s,LGALS9>0&SLC1A5>0&COL1A1>0&KRT14>0);p@meta.data[rownames(p@meta.data)%in%Cells(s),'colocalization']<-'LGALS9/COL1A1_SLC1A5/KRT14'
s<-subset(p,ME=='ICA-ME');s<-subset(s,DPP4>0&CCL3L1>0&COL1A1>0&CD68>0);p@meta.data[rownames(p@meta.data)%in%Cells(s),'colocalization']<-'DPP4/COL1A1_CCL3L1/CD68'
s<-subset(p,ME=='ICA-ME');s<-subset(s,IGF2>0&IGF2R>0&COL1A1>0&KRT14>0);p@meta.data[rownames(p@meta.data)%in%Cells(s),'colocalization']<-'IGF2/COL1A1_IGF2R/KRT14'
SpatialPlot(p,group.by = 'colocalization',alpha = c(0.6,0),ncol=1,images = df[df$file=='KT1_1',]$image,cols = c("red","white"))&guides(fill = guide_legend(override.aes = list(size=5)))


#Extended Data Fig. 12 ------
#Extended Data Fig. 12 a ------
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

#Extended Data Fig. 12 b ------
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

#Extended Data Fig. 12 d ------
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

#Extended Data Fig. 12 e ------
my = subset(sc, L3_C == 'Macrophage')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
mycolor<-cell_color[cell_color$ct%in%levels(Idents(my)),]$color
VlnPlot(my, feature='CXCR4',pt.size=0, cols = mycolor)+NoLegend()
VlnPlot(my, feature='IL1B',pt.size=0, cols = mycolor)+NoLegend()

#Extended Data Fig. 12 f ------
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

#Extended Data Fig. 12 g ------
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


#Extended Data Fig. 12 h ------
my = subset(sc, L3_C == 'Cancer');my$group<-paste(my$metastasis,my$tissue,sep = '_')
Idents(my)<-my$group;my<-subset(my,ident=c('N_Tu','Y_pLN','Y_Tu'))
mhc<-c('HLA-A','HLA-B','HLA-C','TAP1','TAP2','B2M')
my<-AddModuleScore(my,features = list(mhc),name = 'MHCI_score');colnames(my@meta.data)[35]<-'MHCI_score'
my_levels<-c('N_Tu','Y_pLN','Y_Tu');Idents(my)<-factor(Idents(my),levels = my_levels)
my_comparison<-list(c('N_Tu','Y_pLN'),c('Y_pLN','Y_Tu'),c('N_Tu','Y_Tu'))
VlnPlot(my, feature='MHCI_score',pt.size=0, cols = c("#a44a3f","#a44a3f","#a44a3f"),y.max = 6)+NoLegend()+
  stat_compare_means(comparisons = my_comparison)+theme(axis.text.x = element_text(angle = 0,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
                                                        panel.grid.major=element_blank(),panel.grid.minor = element_blank())+ ylab('MHCI Score')+ggtitle('')

#Extended Data Fig. 12 i ------                                                  
my = subset(sc, L3_C == 'Cancer' & tissue == 'Tu');Idents(my)<-my$metastasis
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
VlnPlot(my, feature=c('IL1RN','HLA-DRA'),pt.size=0, cols = c("#377EB8","#E41A1C"))+NoLegend()&theme(axis.text.x = element_text(angle = 0,hjust = 1))

#Extended Data Fig. 12 j ------            
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

#Extended Data Fig. 12 k ------            
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
     

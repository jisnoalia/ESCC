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

#Figure 6 j ------
#left
top_acts_mat_sel = top_acts_mat[rownames(top_acts_mat) %in% c("CC13","CC8","CC2","CC7","CC4","CC14"),
  ]
celltypes = sort(rownames(sp)[grepl("^Mac|^FB",rownames(sp))])
celltypes = c("Epithelium","Cancer","FB.C1.CFD","FB.C2.IGF1" ,"FB.C3.COL1A1" ,"FB.C4.APOE","FB.C5.PDGFRB","FB.C6.ACTA2",
              "Endo.C1.ACKR1", "Endo.C2.FBLN5", "Endo.C3.RGCC","Endo.C4.CCL21",
              "Mac.C1.NLRP3","Mac.C2.SPP1" ,"Mac.C3.C1QC" ,"Mac.C4.LYVE1","Mac.C5.CXCL10","Mac.C6.MKI67",
              "Neutrophil",'CD8.C6.CD39','Platelet','T.C2.MKI67','DC.C3.CD1A',"CD4.C7.OX40"
  )
top_acts_mat_sel = top_acts_mat_sel[,colnames(top_acts_mat_sel) %in% celltypes]
  top_acts_mat_sel = top_acts_mat_sel[,match(celltypes,colnames(top_acts_mat_sel))]
  top_acts_mat_sel
  top_acts_mat_sel = top_acts_mat_sel[match(c("CC13","CC8","CC2","CC7","CC4","CC14"),rownames(top_acts_mat_sel)),]
  top_acts_mat_sel
  palette_length <- 100
  my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
  my_breaks <- c(seq(-1, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 1, length.out=floor(palette_length/2)))
  
  anno_cols = data.frame(ct = rownames(t(top_acts_mat_sel)),
                         Group = 'grey')
  anno_cols$Group = ifelse(anno_cols$ct %in% c("CD4.C7.OX40","CD8.C6.CD39","Endo.C3.RGCC","FB.C3.COL1A1","Mac.C2.SPP1"),"MIMER","other")
  rownames(anno_cols) = anno_cols$ct
  anno_cols$ct = NULL
  annotation_colors = list(
    Group = c('MIMER' = "brown","other" = 'grey')
  )
  pdf("./RCTD/ISCHIA/ISCHIA/Spatial_mapping_of_the_cellular_compostion_cc_15_select_mac_fib.pdf", 
      width = 5, height = 5)
  pheatmap(
    t(top_acts_mat_sel),
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
    cluster_rows = F,cluster_col = T,
    show_colnames = TRUE,
    show_rownames = TRUE,angle_col = 90
  )
  
  dev.off()

#right
celltypes = sort(rownames(sp)[grepl("^Mac|^FB",rownames(sp))])
FetchData(sp,c('cc_15',celltypes),slot = 'data') -> pre_cor_df
head(pre_cor_df)
pre_cor_df_development = pre_cor_df %>% filter(cc_15 %in% c("CC13","CC8","CC2","CC7","CC4","CC14"))

df_result <- pre_cor_df_development %>%
  mutate(
    Mac_rowsum = rowSums(across(starts_with("Mac"))),           
    across(starts_with("Mac"), ~ .x / Mac_rowsum, .names = "{.col}_ratio") , 
    FB_rowsum = rowSums(across(starts_with("FB"))),            
    across(starts_with("FB"), ~ .x / FB_rowsum, .names = "{.col}_ratio") 
  )
head(df_result)

df_result %>%
  ggplot(.,aes(x=FB.C3.COL1A1,y=Mac.C4.LYVE1)) +
  geom_point(size = 2)+
  stat_cor(show.legend = FALSE)+
  geom_smooth(method = "gam", formula = y~x)

library(ggrastr)
ggplot(df_result, aes(x=FB.C3.COL1A1, y=Mac.C4.LYVE1)) +
  geom_point_rast(color="#F37F4F",size=1)+ 
  stat_cor(method = "pearson",size=6,vjust=0.5,digits = 5)+
  stat_smooth(method='gam', 
              color="black", formula = y ~ s(x, k = 10))+ #
  labs(title='')+guides(size=FALSE)+
  theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  theme_classic()+
  theme(axis.text = element_text(size=12,colour = "black"),
        axis.title = element_text(size=15))

#Figure 6 l ------
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
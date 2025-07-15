#--------------------------------
#------Sup-Fig5
#------------------------------

#~~~~~~~~~~~Supplementary Fig.5~~~~~~~~~~~~~
library(ClusterGVis)
library(Seurat)
library(Mfuzz)
library(dplyr)
library(msigdbr)
library(fgsea)

pdf(file="TLS_mfuzz_delnontls.pdf", width=20, height=14)
load('/lustre/home/hycui/duxiao/spatial/spatial_pathology/sp_correct.RData')
df<-read.table("/lustre/home/hycui/duxiao/spatial/CellTrek/sp_table.txt", sep='\t',header = T)
group <- read.csv("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/TLS_cluster_mark.csv", row.names = 1, check.names=F, header = T)

p <- subset(sp, cells = rownames(group))
p <- AddMetaData(p, metadata = group)

dat<-p@assays$spatial@data
d<-NULL
tissue <- c("TLS_G1", "TLS_G2", "TLS_G3", "peri_TLS_G1", "peri_TLS_G2", "peri_TLS_G3")
for(i in tissue){
  cell<-rownames(p@meta.data[p$TLS_Cluster_mark==i,])
  new <- apply(dat[,cell], 1, median)
  d <-cbind(d,new)
}
colnames(d) <-tissue
dat<-Biobase::ExpressionSet(assayData = as.matrix(d))
dat <- Mfuzz::filter.NA(dat)
dat <- Mfuzz::fill.NA(dat, mode = 'mean')
dat <- Mfuzz::filter.std(dat, min.std = 0)
dat <- standardise(dat)
#getClusters(exp = d)
n <- 9
m <- mestimate(dat)
set.seed(123)
cl <- mfuzz(dat, c = n, m = m)
mfuzz.plot2(dat, cl = cl, mfrow = c(4, 3),colo="fancy",time.labels = colnames(cl$centers),xlab = 'TLS_Cluster_mark',ylab = 'Median expression of genes',x11=FALSE)
gene_cluster <- cbind(cl$cluster, cl$membership)
colnames(gene_cluster)[1] <- 'cluster'
gene_cluster<-as.data.frame(gene_cluster);gene_cluster$gene<-rownames(gene_cluster) 

all_gsea <- data.frame()
all_gsea1 <- data.frame()
markGenes <- c("ISG15","STAT1","CXCL14","HLA-B","PSMB8","TAP1","CTSC","B2M","PFN1","IFNGR2","TRAF5","IL1R1","JUN","NFKBIZ","CXCL1","TNFRSF14","TAB2","AURKAIP1","ENO1","GLUL","SAA2","TAPBPL","STMN1","TNFRSF1B","ITGAE","C1QC","UQCRH","CXCL16","BCR","ATP5PB","IRF6","LCK","IGFBP5","JAK1","ACKR1","COL5A2","NEAT1","PTPRC","A2M","FN1","ZFP36L2","ZAP70","CXCR4","PFN2","UQCC2","IGFBP6","IDH2","JUNB","JUND","ZFP36","SPRR3","S100A9","S100A8","IL1RN","SELL","BHLHE40","CSTA","KRT5","EIF4A2","CXCL9","CXCL13","SEC31A","NFKB1","GZMK","PPP2CA","MYC","HK1","BUB3","CD6","DDX6","DDX23","TUBA1A","KRT4","KRT13","KRT15","CDK4","LYZ","PSMD9","DNAJC15","GZMB","PSMA3","BATF","CCL2","IKZF3","STAT3","NME1","DDX42","ICAM2","EXOC7","NKG7")
human_H<-msigdbr(species = "human", category = "H")
human_C2<-msigdbr(species = "human", category = "C2")
human_C2<-human_C2[human_C2$gs_subcat%in%c('CP:BIOCARTA','CP:KEGG','CP:REACTOME'),]
human<-rbind(human_C2,human_H)
msigdbr_list = split(x = human$gene_symbol, f = human$gs_name)

for ( c in c(1:n) ) {
  mm<-gene_cluster[gene_cluster$cluster==c, c('gene','cluster',c)]
  colnames(mm)[3]<-'score'
  mm<-arrange(mm,desc(score))
  tg <- head(mm$gene, n=5)
  rank<-mm[,c('gene','score')]
  geneList<-rank$score
  names(geneList)=rank$gene
  geneList=sort(geneList,decreasing = TRUE)
  fg<-fgsea(msigdbr_list, geneList)
  fg<-fg[fg$leadingEdge%in%fg$leadingEdge[lapply(fg$leadingEdge,length)>3],]
  fg<-arrange(fg,desc(NES))
  fg<-fg[fg$NES>0,]
  gs<-fg[1:5,]
  gs <- gs[!(duplicated(gs$pathway)),]
  r <- gs[,c('pathway', 'pval', 'NES')]
  r <- as.data.frame(r)
  r <- na.omit(r)
  r$group <- paste("C", c, sep="")
  rownames(r) <- r$pathway
  result <- r[,c('group', 'pathway', "pval", "NES")]
  colnames(result)[2] <- "Description"
  colnames(result)[3] <- "pvalue"
  result$Description<-gsub('_',' ',result$Description)
  all_gsea <- rbind(all_gsea, result)
  fg <- fg[!(duplicated(fg$pathway)),]; r1 <- fg[,c('pathway', 'pval', 'NES')]; r1 <- as.data.frame(r1); r1 <- na.omit(r1); r1$group <- paste("C", c, sep=""); rownames(r1) <- r1$pathway; result1 <- r1[,c('group', 'pathway', "pval", "NES")]; colnames(result1)[2] <- "Description"; colnames(result1)[3] <- "pvalue"; result1$Description<-gsub('_',' ',result1$Description); all_gsea1 <- rbind(all_gsea1, result1)
}
cm <- clusterData(exp = d, cluster.method = "mfuzz", cluster.num = n, seed=123)
markGenes <- unique(markGenes)

col1 <- c("#D51F26","#272E6A", "#208A42","#89288F","#F47D2B","#9983BD","#8A9FD1","#C06CAB", "#D24B27")
times <- c(10,9,5,8,6,3,8,1,8)
col <-rep(col1,times)
need_gsea <- read.table("TLS_mfuzz_markgsva.txt", header=T, sep="\t")
visCluster(object = cm, plot.type = "both", column_names_rot = 45, show_row_dend = F, markGenes = markGenes, markGenes.side = "left", annoTerm.data = need_gsea, line.side = "left", go.col = col, add.bar = T, textbar.pos = c(0.85,0.15), go.size=10, ctAnno.col=col1)
dev.off()
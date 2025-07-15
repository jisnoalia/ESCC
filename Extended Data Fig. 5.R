#Extended Data Fig. 5 -----
#Extended Data Fig. 5 b -----
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('Macrophage','Monocyte'));Idents(my)<-my$L4_C
is <- c('CD274','PDCD1LG2','ICOSLG','CD276','VTCN1','VSIR','IDO1','TGFB1','IL10','CCL5','CCL17','CCL22','CXCL8',
        'CCL16','CCL18','IL17A','CCL24','IL1R1','IL4','IL13','ARG1','CD40LG','TPH1','XBP1','SOCS1')# Immunosuppression
ag <- c('VEGFA','CXCL8','CXCL1','CXCL2','TEK','ANGPT2','CXCL12','FLT1','CXCR4','CTSB','CTSS','SEMA4D','PROK2',
        'CXCR1','S1PR1','PDGFA','FGF2','MMP12','VEGFC','DNTTIP2','PGF','HIF1A','NFE2L2','CCND2','CCNE1','CD44',
        'E2F3','EDN1','EZH2','FGF18','FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2','NOTCH1','PTK2','SPP1','STC1',
        'TNFAIP6','TYMP','VAV2')# Angiogenesis
im <- c('MMP7','MMP2','MMP3','MMP9','CCL2','CCL3','CCL22','TGFB1','EGF','CCR2','FLT1','CTSS','CTSB','CXCR3',
        'WNT7B','ALOX5','CDH1','CCL18','CCL24','S1PR1','CXCL16','MAPK7','CSF1','SPARC','TLR4','VCAM1','CCL20')# Invasion and metastasis
ts <- c('IL1B','IL6','IL12A','IL23A','TNF','CXCL9','CXCL10','TLR2','TLR4','CXCL11','IFNG','CD40','FCGR2A',
        'ITGAX','IFNGR1','HLA-DPB1','HLA-DPA1','HLA-DRA','HLA-DRB1','HLA-DQA1','CD74','HLA-DRB5','IRF5')# Tumor suppression
pc <- c('MRC1','CD163','MERTK','C1QB','FCRLA','CD5L','CD81','GPNMB','CD36')# Phagocytosis
my = AddModuleScore(my,features=list(is),name='is'); my = AddModuleScore(my,features=list(ag),name='ag')
my = AddModuleScore(my,features=list(im),name='im'); my = AddModuleScore(my,features=list(ts),name='ts'); my = AddModuleScore(my,features=list(pc),name='pc')
my@meta.data = dplyr::rename(my@meta.data,ImmunoSuppression=is1, Angiogenesis=ag1, Metastasis=im1,
                            TumorSuppression=ts1, Phagocytosis=pc1)
fea <- c('ImmunoSuppression','Angiogenesis','Metastasis','TumorSuppression','Phagocytosis')
myd=as.data.frame(scale(my@meta.data[,fea])); myd=cbind(myd,my@meta.data[,c('tissue','L4_C')]); head(myd)
myd=myd[,-6]; myd=dplyr::group_by(myd, L4_C) %>% dplyr::summarise(across(all_of(fea), mean))
myd=as.data.frame(myd);myd=tibble::column_to_rownames(myd,var='L4_C'); myd=t(as.matrix(myd))
pheatmap(myd, cellwidth=30, cellheight=30, cluster_cols=TRUE, clustering_method='ward.D',
         cluster_rows=T, treeheight_row=0, kmeans_k=NA, border_color="white", scale="none", drop_levels=T,
         show_rownames=TRUE, show_colnames=TRUE, color=rev(getPalette(10)),
         fontsize=12, fontsize_row=15, legend=T, fontsize_col=15, display_numbers=ifelse(myd > 0.37,"*",""))

ma <- c('IL23A','TNF','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','CD74','HLA-DPA1','HLA-DPB1',
        'HLA-DRB5','HLA-DRB1','CD83','CD68','CD80','S100A8','S100A9')# M1
mb <- c('CCL20','CCL18','CCL22','LYVE1','VEGFA','CTSB','CTSD','TGFB1','MMP19','MMP9','CLEC7A','WNT7B','FASLG',
        'TNFSF12','TNFSF8','CD276','MSR1','FN1','IRF4','MRC1','CD163','ARG1','IL10','MARCO')# M2
my = AddModuleScore(my,features=list(ma),name='ma'); my = AddModuleScore(my,features=list(mb),name='mb')
my@meta.data = dplyr::rename(my@meta.data,Inflammatory=ma1, 'Anti-inflammatory'=mb1)
fea <- c('Anti-inflammatory','Inflammatory'); myd = as.data.frame(scale(my@meta.data[,fea]))
myd = cbind(myd,my@meta.data[,c('tissue','L4_C')]); myd = myd[,-3]
myd = dplyr::group_by(myd, L4_C) %>% dplyr::summarise(across(all_of(fea), mean))
myd = as.data.frame(myd); myd = tibble::column_to_rownames(myd,var='L4_C'); myd = t(as.matrix(myd))
pheatmap(myd, cellwidth=30, cellheight=30, cluster_cols=TRUE, clustering_method='ward.D',
         cluster_rows=T, treeheight_row=0, kmeans_k=NA, border_color="white", scale="none", drop_levels=T,
         show_rownames=TRUE, show_colnames=TRUE, color=rev(getPalette(10)),
         fontsize=12, fontsize_row=15, legend=T, fontsize_col=15)

#Extended Data Fig. 5 d-----
#top
s1 = subset(sp,development =='unknown',invert=T)
  ST.exp = s1@assays$spatial@data %>% as.data.frame()
  spot_abundance_category = s1@meta.data %>% as.matrix()
  
  ligand_receptor_df1 =  ligand_receptor_df
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
    # i =1 
    ct = ligand_receptor_df1[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  
  ligand_receptor_df1$lrpair <- paste0(ligand_receptor_df1$ligand, "_", ligand_receptor_df1$receptor)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df1, seq(nrow(ligand_receptor_df1)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) # step 2
  rownames(LigandR_mean.m) -> turn1
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("cc_15")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "cc_15"
  data_long1 <- reshape2::melt(data.plot, id.vars = c("cellID", "cc_15"), measure.vars = ligand_receptor_df1$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long1[, "cc_15"] <- factor(data_long1[, "cc_15"])
  
  
  data_long = data_long1
  turn = gsub("\\_$","",
              paste0(ligand_receptor_df$ligand,"_",ligand_receptor_df$receptor))
  data_long$LigandReceptor = factor(data_long$LigandReceptor,
                                    levels = turn)
  # data_long$cc_15 = ifelse(!data_long$cc_15 %in% c("CC3","CC5","CC9","CC10"),'other',data_long$cc_15)
  
  data_long_fil = data_long %>%  filter(cc_15 %in% c("CC3","CC5","CC9","CC10"))
  data_long_fil$cc_15 = factor(data_long_fil$cc_15,
                               levels = c("CC3","CC5","CC9","CC10",'other'))
  all_ccs <- unique(sp$cc_15)
  color_mapping15 <- setNames(dittoSeq::dittoColors()[1:15][1:length(all_ccs)], sort(all_ccs))
  
  my_comparsion = list(
    c("CC3","CC5"),
    c("CC5","CC9"),
    c("CC9","CC10")
  )
  g <- ggplot(data_long_fil,aes(x=LigandReceptor,y= LRmean,fill=cc_15))+ 
    # geom_point()+
    geom_boxplot(width = 0.4, outliers  = F, position = position_dodge(width = 0.7)) +
    # facet_grid(~LigandReceptor) + 
    scale_fill_manual(values = color_mapping15)+
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    stat_compare_means(mapping = aes(group = cc_15),label.y = 2,show.legend = F,label = "p.signif",size=2,step.increase = 0.05)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 8,color = "black"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  g

##bottom
 s1 = subset(sp,ME =='unknown',invert=T)
  ST.exp = s1@assays$spatial@data %>% as.data.frame()
  spot_abundance_category = s1@meta.data %>% as.matrix()
  
  ligand_receptor_df1 =  ligand_receptor_df
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
    # i =1 
    ct = ligand_receptor_df1[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  
  ligand_receptor_df1$lrpair <- paste0(ligand_receptor_df1$ligand, "_", ligand_receptor_df1$receptor)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df1, seq(nrow(ligand_receptor_df1)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) # step 2
  rownames(LigandR_mean.m) -> turn1
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("cc_15")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "cc_15"
  data_long1 <- reshape2::melt(data.plot, id.vars = c("cellID", "cc_15"), measure.vars = ligand_receptor_df1$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long1[, "cc_15"] <- factor(data_long1[, "cc_15"])
  
  
  data_long = data_long1
  turn = gsub("\\_$","",
              paste0(ligand_receptor_df$ligand,"_",ligand_receptor_df$receptor))
  data_long$LigandReceptor = factor(data_long$LigandReceptor,
                                    levels = turn)
  # data_long$cc_15 = ifelse(!data_long$cc_15 %in% c("CC3","CC5","CC9","CC10"),'other',data_long$cc_15)
  
  data_long_fil = data_long %>%  filter(cc_15 %in% c("CC3","CC5","CC9","CC10"))
  data_long_fil$cc_15 = factor(data_long_fil$cc_15,
                               levels = c("CC3","CC5","CC9","CC10",'other'))
  all_ccs <- unique(sp$cc_15)
  color_mapping15 <- setNames(dittoSeq::dittoColors()[1:15][1:length(all_ccs)], sort(all_ccs))
  
  my_comparsion = list(
    c("CC3","CC5"),
    c("CC5","CC9"),
    c("CC9","CC10")
  )
  g <- ggplot(data_long_fil,aes(x=LigandReceptor,y= LRmean,fill=cc_15))+ 
    # geom_point()+
    geom_boxplot(width = 0.4, outliers  = F, position = position_dodge(width = 0.7)) +
    # facet_grid(~LigandReceptor) + 
    scale_fill_manual(values = color_mapping15)+
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    stat_compare_means(mapping = aes(group = cc_15),label.y = 2,show.legend = F,
                       # label = "p.signif",
                       size=2,step.increase = 0.05
                       ,method = 'wilcox.test'
                       )+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 8,color = "black"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  g
#Extended Data Fig. 5 e-----
files = c('D_JT2','IT1','IT2','IT3','KT2','PT1','ST1','TT2',"QN","QP","QT1",'MLN_P','PT3LN_P_2','ZLN_P','YLN_3LN_P','QLn_PLn_N',
            'YN','YP','YT1',"YT2")
  LR_list = list()
  for (file  in files ) {
    p0 <- SpatialDimPlot(sp,group.by = 'new_mac_cc',images = df[df$file==file,'image'],stroke = NA,pt.size.factor = 1.4,image.alpha = 0)+
      scale_fill_manual(values = c(color_mapping15[ c("CC3","CC5","CC9","CC10")],'other'='lightgrey'),name='')+
      theme_minimal() +
      coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
      scale_x_continuous(breaks = common_xbreaks) +
      scale_y_continuous(breaks = common_xbreaks) +
      theme(
        aspect.ratio = 1,legend.position = 'top',legend.direction = 'horizontal')
    
    p1 <- SpatialFeaturePlot(sp, features =  "HLA_A_LILRB1", images = df[df$file==file,"image"],stroke = NA,
                             # alpha = c(1,0.6),
                             image.alpha = 0,pt.size.factor = 1.4) +
      theme_minimal() +
      scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),na.value='lightgrey',name = 'HLA-A_LILRB1'
      )+
      coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
      scale_x_continuous(breaks = common_xbreaks) +
      scale_y_continuous(breaks = common_xbreaks) +
      theme(
        aspect.ratio = 1,
        legend.title.position = "top" ,
        legend.justification = "center",
        legend.position = "top")
    
    p2 <- SpatialFeaturePlot(sp, features =  "CD47_SIRPA", images = df[df$file==file,"image"],stroke = NA,
                             # alpha = c(1,0.6),
                             image.alpha = 0,pt.size.factor = 1.4) +
      theme_minimal() +
      scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu")),na.value = "lightgrey"
      )+
      coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
      scale_x_continuous(breaks = common_xbreaks) +
      scale_y_continuous(breaks = common_xbreaks) +
      theme(
        aspect.ratio = 1,
        legend.title.position = "top" ,
        legend.justification = "center",
        legend.position = "top")
   p0 = p0 +labs(title = file)
   p =  p0|p1|p2 
   LR_list[[file]] = p 
   }
  
  for (file  in files ) {
  p_all = LR_list[[file]]
  ggsave(plot = p_all,filename = paste0(file,"_spatial_LR.pdf"),width = 12,height = 4)
  
  }
#Extended Data Fig. 5 f-----
load('sp_correct.RData')
df<-read.table("sp_table.txt", sep='\t',header = T)
my<-subset(sp,development=='unknown',invert=T)
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
myd<-myd[,c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')]
pheatmap(myd,color=rev(getPalette(10)),border_color = "grey60", fontsize = 12, angle_col = "90",clustering_method = "ward.D2",cluster_cols = F,cluster_rows = T,scale = "row")
#Extended Data Fig. 5 g-----
library(ggplot2)
library(cowplot)
pdf(file="HALLMARK_ME_DEV.pdf", width=18, height=8)
meta <- read.table("HALLMARK_meta.data.txt", header=T, sep="\t")

hallmark <- as.character(colnames(meta[84:133]))
m1 <- margin(t=0.5, r=0.1, b=0.5, l=0.1, unit="cm")
theme3 <- theme(plot.title = element_text(size = 16, face = "bold"), axis.title.y = element_text(size=15), axis.title.x = element_text(size=15), axis.text.x = element_text(size = 12), axis.text.y = element_text(size =14), legend.title = element_blank(), legend.text = element_text(size = 12), legend.position="right", plot.margin=m1, panel.background = element_rect(fill='white', color = 'black'))

plots <- vector("list", 50)
plots1 <- vector("list", 50)
c <- 0
for (h in hallmark) {
    h1 <- gsub("HALLMARK_", "", h)
    c = c + 1
    meta1 <- meta[meta$development != "unknown",]
    meta1$development <- factor(meta1$development, levels = c('Nor', 'Hyp', 'MiD', 'MoD', 'SD&CA', 'ICA', 'MCA'))
    de <- aggregate(meta1[,h], by=list(type=meta1$development),FUN=function(x) c(mean=mean(x), sd=sd(x), n <- length(x), sq=sqrt(n))); de$dev <- de$type; de$group = "DEV"; de$mean <- de$x[,1]; de$sd <- de$x[,2]; de$sem <- de$x[,2]/de$x[,4]

    meta2 <- meta[meta$ME != "unknown",]
    meta2$ME <-factor(meta2$ME,levels = c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA'))
    me <- aggregate(meta2[,h], by=list(type=meta2$ME),FUN=function(x) c(mean=mean(x), sd=sd(x), n <- length(x), sq=sqrt(n))); me$dev <- gsub("", "", me$type); me$group = "ME"; me$mean <- me$x[,1]; me$sd <- me$x[,2]; me$sem <- me$x[,2]/me$x[,4]

    merge <- rbind(de, me)
    merge$dev <- factor(merge$dev, levels = c('Nor', 'Hyp', 'MiD', 'MoD', 'SD&CA', 'ICA', 'MCA'))

    plots[[c]] <- ggplot(data=merge, aes(x=dev, y=mean, group=group)) + geom_point(aes(color=group)) + geom_line(aes(color=group)) + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, alpha=0.7) + labs(title=h1, x="Dev stages", y="Score (mean+-sem)") + theme3

    de1 <-  aggregate(meta1[,h], by=list(type=meta1$development,file=meta1$file),mean); de1$lab <- paste(de1$type, de1$file, sep="-")

    me1 <- aggregate(meta2[,h], by=list(type=meta2$ME, file=meta2$file),mean); me1$dev <- gsub("", "", me1$type); me1$lab <- paste(me1$dev, me1$file, sep="-")

    de2 <- de1[de1$lab %in% me1$lab,]
    merge1 <- cbind(de2, me1)
    merge1$ratio <- merge1[,3] / merge1[,7]
    merge1$dev <- factor(merge1$dev, levels = c('Nor', 'Hyp', 'MiD', 'MoD', 'SD&CA', 'ICA', 'MCA'))
    de2me <- aggregate(merge1$ratio, by=list(dev1=merge1$dev), FUN=function(x) c(mean=mean(x), sd=sd(x), n <- length(x), sq=sqrt(n)));
    de2me$mean <- de2me$x[,1]; de2me$sd <- de2me$x[,2]; de2me$sem <- de2me$x[,2]/de2me$x[,4]
    plots1[[c]] <- ggplot(data=de2me, aes(x=dev1, y=mean, group=1)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, alpha = 0.7) + labs(title=h1, x="Dev stages", y="Ratio+-sem") + theme3
}

for ( i in c(1,7,13,19,25,31,37,43) ) {
    j = i + 5
    p <- plot_grid(plotlist = plots[c(i:j)], ncol = 3, nrow = 2)
    print(p)
}
plot_grid(plotlist = plots[c(49,50)], ncol = 3, nrow = 2)

for ( i in c(1,7,13,19,25,31,37,43) ) {
    j = i + 5
    p <- plot_grid(plotlist = plots1[c(i:j)], ncol = 3, nrow = 2)
    print(p)
}
plot_grid(plotlist = plots1[c(49,50)], ncol = 3, nrow = 2)
dev.off()


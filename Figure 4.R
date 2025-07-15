#Figure 4a -----
load("./permutation_All.RData")
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
  geom_boxplot(width=0.5, outlier.shape = NA, alpha = 0.8)+ 
  scale_fill_manual(values = mycolor)+ stat_compare_means(paired = FALSE, comparisons = my_comparison) + 
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(size = 14)) + 
  theme_classic() + 
  NoLegend()

p+coord_fixed(ratio = 0.5)

# Figure 4b -----
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


# Figure 4c -----
load("cpdb_result.Rdata")
ggplot(cpdb,aes(x = LR,y = gene_pair, fill = mean,size = logPval)) +
  geom_point(shape=21,color='black') +
  scale_fill_gradientn(colors = color) +
  scale_y_discrete(limits=rev(levels(plot_df$subC))) +
  scale_size_continuous(range = c(1,8),name = 'Proportion') +
  facet_grid(~type,space = 'free_x',scale = 'free_x',labeller = label_wrap_gen(width=5)) +
  theme_bw() + xlab('') + ylab('') +
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_text(size = 10))
ligand_receptor_df = data.frame(ligand  = c("WNT5A","VEGFA","VEGFA","TNFSF4","TIMP1","TGFB1","SPP1","SPP1","SPP1","SPP1","SPP1","PVR","PDCD1","NECTIN2","MDK","MDK","LGALS9","IGF2","FN1","FN1","CXCL8","COL5A1",'COL4A2','COL4A1',"COL3A1","COL1A2",'COL1A1','COL1A1','COL1A1',"ANXA1"),
                                receptor = c("ANTXR1","KDR",'FLT1','TNFRSF4','FGFR2','TGFBR1',"PTGER4","CD44",'CCR8','ITGA9','ITGA4',"TIGIT","FAM3C",'TIGIT','SORL1','LRP1','HAVCR2','IGF2R','ITGA5','ITGA4','ACKR1','ITGA2','ITGA2','ITGA2','ITGA2','ITGA2','ITGA2','ITGA1','ITGA10',"FPR1"),
                                receptor2 = c("","","","","","TGFBR2","","","","ITGB1","ITGB1",'','','','',"","","","ITGB1","ITGB1",'',"ITGB1","ITGB1","ITGB1","ITGB1","ITGB1","ITGB1","ITGB1","ITGB1",""))
s1 = subset(sp,development =='unknown',invert=T)
  # slide = df[df$file == p,"image"]
  # # sample = df[df$file == p,"image"]
  # s1@images = s1@images[slide]
  
  ST.exp = s1@assays$spatial@data %>% as.data.frame()
  spot_abundance_category = s1@meta.data %>% as.matrix()
  
  ligand_receptor_df1 =  ligand_receptor_df %>% filter(receptor2 == "")
  ligand_receptor_df1$receptor2 = NULL
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
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
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("mimer_score_spot")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "mimer_score_spot"
  data_long1 <- reshape2::melt(data.plot, id.vars = c("cellID", "mimer_score_spot"), measure.vars = ligand_receptor_df1$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long1[, "mimer_score_spot"] <- factor(data_long1[, "mimer_score_spot"])
  
  
  
  ligand_receptor_df2 =  ligand_receptor_df %>% filter(receptor2 != "")
  # ligand_receptor_df2$receptor2 = NULL
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df2), function(i){
    ct = ligand_receptor_df2[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    
    return(tmp)
  })
  ligand_receptor_df2$lrpair <- paste0(ligand_receptor_df2$ligand, "_", ligand_receptor_df2$receptor,"_",ligand_receptor_df2$receptor2)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df2, seq(nrow(ligand_receptor_df2)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2],"_",x[3]))) # step 2
  rownames(LigandR_mean.m) -> turn2
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("mimer_score_spot")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "mimer_score_spot"
  data_long2 <- reshape2::melt(data.plot, id.vars = c("cellID", "mimer_score_spot"), measure.vars = ligand_receptor_df2$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long2[, "mimer_score_spot"] <- factor(data_long2[, "mimer_score_spot"])
  
  data_long = rbind(data_long1,data_long2)
  turn = gsub("\\_$","",
       paste0(ligand_receptor_df$ligand,"_",ligand_receptor_df$receptor,"_",ligand_receptor_df$receptor2))
  data_long$LigandReceptor = factor(data_long$LigandReceptor,
                                    levels = turn)
  g <- ggplot(data_long,aes(x=LigandReceptor,y= LRmean,fill=mimer_score_spot))+ 
    geom_boxplot(width = 0.4, outlier.shape = NA, lwd = 0.3, position = position_dodge(width = 0.7)) +
    # facet_grid(~LigandReceptor) + 
    scale_fill_manual(values = c("mimer"= 'brown','other' = 'grey'))+
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    stat_compare_means(aes(group=mimer_score_spot),label.y = 3,show.legend = F,label = "p.signif",size=2)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 8),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8))+ylim(0,3.5)
  g
s1 = subset(sp,ME =='unknown',invert=T)
  slide = df[df$file == p,"image"]
  # # sample = df[df$file == p,"image"]
  # s1@images = s1@images[slide]
  
  ST.exp = s1@assays$spatial@data %>% as.data.frame()
  spot_abundance_category = s1@meta.data %>% as.matrix()
  
  ligand_receptor_df1 =  ligand_receptor_df %>% filter(receptor2 == "")
  ligand_receptor_df1$receptor2 = NULL
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
    ct = ligand_receptor_df1[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  ligand_receptor_df1$lrpair <- paste0(ligand_receptor_df1$ligand, "_", ligand_receptor_df1$receptor)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df1, seq(nrow(ligand_receptor_df1)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) # step 2
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("mimer_score_spot")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "mimer_score_spot"
  data_long1 <- reshape2::melt(data.plot, id.vars = c("cellID", "mimer_score_spot"), measure.vars = ligand_receptor_df1$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long1[, "mimer_score_spot"] <- factor(data_long1[, "mimer_score_spot"])
  
  
  
  ligand_receptor_df2 =  ligand_receptor_df %>% filter(receptor2 != "")
  # ligand_receptor_df2$receptor2 = NULL
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df2), function(i){
    ct = ligand_receptor_df2[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  ligand_receptor_df2$lrpair <- paste0(ligand_receptor_df2$ligand, "_", ligand_receptor_df2$receptor,"_",ligand_receptor_df2$receptor2)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df2, seq(nrow(ligand_receptor_df2)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2],"_",x[3]))) # step 2
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("mimer_score_spot")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "mimer_score_spot"
  data_long2 <- reshape2::melt(data.plot, id.vars = c("cellID", "mimer_score_spot"), measure.vars = ligand_receptor_df2$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long2[, "mimer_score_spot"] <- factor(data_long2[, "mimer_score_spot"])
  
  data_long = rbind(data_long1,data_long2)
  
  turn = gsub("\\_$","",
              paste0(ligand_receptor_df$ligand,"_",ligand_receptor_df$receptor,"_",ligand_receptor_df$receptor2))
  data_long$LigandReceptor = factor(data_long$LigandReceptor,
                                    levels = turn)
  
  g <- ggplot(data_long,aes(x=LigandReceptor,y= LRmean,fill=mimer_score_spot))+ 
    geom_boxplot(width = 0.4, outlier.shape = NA, lwd = 0.3, position = position_dodge(width = 0.7)) +
    # facet_grid(~LigandReceptor) + 
    scale_fill_manual(values = c("mimer"= 'brown','other' = 'grey'))+
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    stat_compare_means(aes(group=mimer_score_spot),label.y = 3,show.legend = F,label = "p.signif",size=2)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 8),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8))+ylim(0,3.5)
  g
# Figure 4d -----
library(spatstat)
library(png)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)
library(Seurat)
library(ggpubr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggraph)
library(tidygraph)
codex2 = codex
names(codex2@images) = unique(codex2$Sample)

unique_types_cols = rep('grey',length(unique_types))
names(unique_types_cols) = unique_types

Mimer = c("SPP1+Mac",'ISG15+Mac','Ki-67+Mac',
          'CD4T-aTreg','Tcells.cc',
          'PDCD1+Tex','PDCD1+Tex2','LAG3+Tex',
          'COL1A1+CAF','POSTN+CAF',
          'PLVAP+Endo','Vimentin+Endo')
unique_types_cols[names(unique_types_cols) %in% Mimer] = 'red'
CN_ratio_4 %>%
  distinct(Sample,neighborhood10) -> tt
names(table(tt$Sample))[table(tt$Sample)==2] -> use_2_dist

table(codex2$neighborhood10,codex$Sample) %>% as.data.frame() %>% 
  rename('neighborhood10'='Var1','Sample'='Var2','number'='Freq') %>%
  filter(neighborhood10 %in% c("CN4","CN12")) %>% filter()
unique_types <- unique(Mimer) 
CNS = c('CN4',"CN12")
type_dist_df_CN_sample_list <- list()
Samples = use_2_dist
for (sa in Samples) {
  location = GetTissueCoordinates(codex2@images[[sa]])
  location %>%
    column_to_rownames("cell") -> pos
  
  dist_mat <- as.matrix(dist(pos))
  dist_mat <- (dist_mat + t(dist_mat))/2
  diag(dist_mat) <- 0
  colnames(dist_mat) <- rownames(pos)
  rownames(dist_mat) <- rownames(pos)
  distance_matrix_df = dist_mat
  dist_mat[1:4,1:4]
  
  type_combinations <- expand.grid(type1 = unique_types, type2 = unique_types) %>%
    filter(as.numeric(type1) <= as.numeric(type2)) # 确保只处理一次 (A-B 和 B-A 算一次)
  cell_ids = codex2@meta.data %>%
    mutate(cell_id = rownames(.),
           cell_type = subCelltype) %>% 
    filter(Sample==sa) %>%
    dplyr::select(cell_id,cell_type,neighborhood10)
  # 用于存储结果的列表
  type_dist_df_CN_list <- list()
  # CNS = c('CN4',"CN12")
  for (CN in CNS) {
    cell_id_to_type = cell_ids%>%
      filter(neighborhood10 == CN)
    type_dist_results_list <- list()
    
    # 遍历所有细胞类型组合
    for (i in 1:nrow(type_combinations)) {
      type1 <- type_combinations$type1[i] %>% as.character()
      type2 <- type_combinations$type2[i] %>% as.character()
      # 获取属于当前类型的细胞ID
      cells_type1_ids <- cell_id_to_type %>%
        filter(cell_type == type1) %>%
        pull(cell_id)
      cells_type2_ids <- cell_id_to_type %>%
        filter(cell_type == type2) %>%
        pull(cell_id)
      # 确保有细胞属于这些类型
      if (length(cells_type1_ids) == 0 || length(cells_type2_ids) == 0) {
        next # 跳过没有细胞的类型
      }
      
      # 从距离矩阵中提取所有 Type1 和 Type2 之间两两距离的子集
      # 使用行名和列名进行子集选择
      sub_matrix <- distance_matrix_df[cells_type1_ids, cells_type2_ids, drop = FALSE]
      
      # 将子矩阵展平为向量
      flat_distances <- as.vector(sub_matrix)
      
      # 如果是同一种类型 (Type1 == Type2)，我们需要处理对角线 (自距离为0)
      if (type1 == type2) {
        # 移除值为0的自连接距离
        flat_distances <- flat_distances[flat_distances != 0]
      }
      
      # 计算各种距离统计量
      if (length(flat_distances) > 0) {
        mean_dist <- mean(flat_distances)
        min_dist <- min(flat_distances)
        max_dist <- max(flat_distances)
      } else {
        mean_dist <- NA
        min_dist <- NA
        max_dist <- NA
      }
      
      # 存储结果
      type_dist_results_list[[i]] <- data.frame(
        Type1 = type1,
        Type2 = type2,
        Mean_Distance = mean_dist,
        Min_Distance = min_dist,
        Max_Distance = max_dist
      )
    }
    
    # 将结果列表合并为数据框
    type_dist_df <- do.call(rbind, type_dist_results_list)
    
    type_dist_df$CN = CN
    
    type_dist_df_CN_list[[CN]] = type_dist_df
  }
  type_dist_df_CN <- do.call(rbind,type_dist_df_CN_list)
  type_dist_df_CN$Sample = sa
  type_dist_df_CN_sample_list[[sa]] = type_dist_df_CN
  print(paste0(i,CN,sa,"done!!!"))
}

type_dist_df_CN_sample_list_df = do.call(rbind,type_dist_df_CN_sample_list)

colnames(codex@meta.data)
codex@meta.data %>%
  group_by(neighborhood10,Sample,subCelltype) %>% 
  summarise(number = n()) %>%
  filter(subCelltype%in%'PDCD1+Tex')

type_dist_df_CN_sample_list_df %>%
  group_by(CN,Type1,Type2) %>%
  filter(!is.na(Min_Distance) )%>% #Min_Distance
  summarise(mean_Mean_Distance = mean(Min_Distance)) %>%  #Min_Distance
  # filter(!Type1%in%'PDCD1+Tex') %>%
  # filter(!Type2%in%'PDCD1+Tex') %>%
  filter(CN=='CN4') -> patient_data


codex@meta.data %>% 
  group_by(neighborhood10,subCelltype) %>%
  summarise(number = n()) %>% 
  filter(neighborhood10=='CN4') %>%
  filter(subCelltype %in% Mimer) -> CN_number_count


edges_df <- patient_data %>%
  select(from = Type1, to = Type2, weight = mean_Mean_Distance) %>%
  filter(from != to)


    edges_df <- patient_data %>%
      select(from = Type1, to = Type2, weight = mean_Mean_Distance) %>%
      # 确保没有自环，除非你想表示自己与自己的距离
      filter(from != to)
    
    # 创建所有可能出现的细胞类型节点列表
    nodes_df <- data.frame(name = unique(c(edges_df$from, edges_df$to)))
    nodes_df <- nodes_df %>% mutate(subCelltype = name) %>%left_join(CN_number_count,by = 'subCelltype')
    
    # 创建 igraph 对象
    # 如果你的细胞类型是分类变量，确保factor levels一致
    g_patient <- tbl_graph(nodes = nodes_df, edges = edges_df)
    # g_patient <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)
    
    # 为节点添加一些属性，例如，可以为所有图统一分配颜色
    V(g_patient)$type_group <- V(g_patient)$name # 这里节点名称就是细胞类型
    
    # 计算 'proximity' 来控制边的视觉属性
    max_dist_val <- max(E(g_patient)$weight)
    E(g_patient)$proximity <- max_dist_val - E(g_patient)$weight # 距离越小，proximity越大
    
    # 绘制 ggraph
    
    gr1_layout <- create_layout(g_patient, layout = "linear",circular = T)
    p_patient <- ggraph(gr1_layout) +
      geom_edge_link(aes(alpha = proximity, width = proximity),
                     edge_colour = "grey", lineend = "round") +
      scale_edge_width_continuous(range = c(0.5, 1)) +
      scale_edge_alpha_continuous(range = c(0.4, 1)) +
      geom_node_point(aes(color=name,size=log10(number))) + # 节点大小可以放大一些 #aes(color = type_group), 
      geom_node_text(aes(label = name), repel = TRUE, size = 4) + # 显示细胞类型名称
      # scale_color_brewer(palette = "Set2") + # 为细胞类型选择调色板
      # scale_size_continuous(range = c(0.5,10))+
      scale_color_manual(values = cell_color)+
      labs(title = paste0("Cell Type Network for CN4"),
           edge_alpha = "Proximity", edge_width = "Proximity", color = "Cell Type",size="log10(Cell Count)") +
      theme_graph() +
      theme(legend.position = "right") + guides(color=guide_legend(ncol=2, override.aes = list(size=4)),
                                                size=guide_legend(ncol=2),
                                                )
    
    p_patient

    
    {
      type_dist_df_CN_sample_list_df %>%
        group_by(CN,Type1,Type2) %>%
        filter(!is.na(Min_Distance) )%>%
        summarise(mean_Mean_Distance = mean(Min_Distance)) %>%
        filter(!Type1%in%'PDCD1+Tex') %>%
        filter(!Type2%in%'PDCD1+Tex') %>%
        filter(CN=='CN12') -> patient_data
      
      patient_data %>%
        filter(Type1 %in% Mimer & Type2 %in% Mimer) -> filter_CN12

      codex@meta.data %>% 
        group_by(neighborhood10,subCelltype) %>%
        summarise(number = n()) %>% 
        filter(neighborhood10=='CN12') %>%
        filter(subCelltype %in% Mimer) -> CN_number_count

      edges_df <- patient_data %>%
        select(from = Type1, to = Type2, weight = mean_Mean_Distance) %>%
        filter(from != to)
      library(ggraph)
      library(tidygraph)
      
      # 准备 igraph 需要的边列表 (edge list)
      # 这里我们使用 Mean_Distance 作为权重
      patient_data$mimer=ifelse(patient_data$Type1 %in% Mimer & patient_data$Type1 %in% Mimer ,'Mimer','other')
      edges_df <- patient_data %>%
        select(from = Type1, to = Type2, weight = mean_Mean_Distance) %>%
        # 确保没有自环，除非你想表示自己与自己的距离
        filter(from != to)
      
      # 创建所有可能出现的细胞类型节点列表
      nodes_df <- data.frame(name = unique(c(edges_df$from, edges_df$to)))
      nodes_df <- nodes_df %>% mutate(subCelltype = name) %>%left_join(CN_number_count,by = 'subCelltype')
      
      # 创建 igraph 对象
      # 如果你的细胞类型是分类变量，确保factor levels一致
      g_patient <- tbl_graph(nodes = nodes_df, edges = edges_df)
      # g_patient <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)
      
      # 为节点添加一些属性，例如，可以为所有图统一分配颜色
      V(g_patient)$type_group <- V(g_patient)$name # 这里节点名称就是细胞类型
      
      # 计算 'proximity' 来控制边的视觉属性
      max_dist_val <- max(E(g_patient)$weight)
      E(g_patient)$proximity <- max_dist_val - E(g_patient)$weight # 距离越小，proximity越大
      
      # 绘制 ggraph
      
      gr1_layout <- create_layout(g_patient, layout = "linear",circular = T)
      p_patient12 <- ggraph(gr1_layout) +
        geom_edge_link(aes(alpha = proximity, width = proximity),
                       edge_colour = "grey", lineend = "round") +
        scale_edge_width_continuous(range = c(0.5, 1)) +
        scale_edge_alpha_continuous(range = c(0.4, 1)) +
        geom_node_point(aes(color=name,size=log10(number)) )+ # 节点大小可以放大一些 #aes(color = type_group), 
        geom_node_text(aes(label = name), repel = TRUE, size = 4) + # 显示细胞类型名称
        # scale_color_brewer(palette = "Set2") + # 为细胞类型选择调色板
        # scale_size_continuous(range = c(0,10),limits = c(0,10))+
        scale_color_manual(values = cell_color)+
        labs(title = paste0("Cell Type Network for CN12"),
             edge_alpha = "Proximity", edge_width = "Proximity", color = "Cell Type",size="log10(Cell Count)") +
        theme_graph() +
        theme(legend.position = "right") + guides(color=guide_legend(ncol=4, override.aes = list(size=4)),
                                                  size=guide_legend(ncol=2),
        )
     
}
p_patient|p_patient12
    
# Figure 4f -----
df<-data.frame('Score Type'=c('Immunosuppression','Angiogenesis','Invasion and metastasis','Tumor suppression','Phagocytosis','Inflammatory macrophage','Anti-inflammatory macrophage','Exhaustion','Toxicity','MHCI score'))
df$gene<-'';df[1,]$gene<-paste('CD274','PDCD1LG2','ICOSLG','CD276','VTCN1','VSIR','IDO1','TGFB1','IL10','CCL5','CCL17','CCL22','CXCL8',
                               'CCL16','CCL18','IL17A','CCL24','IL1R1','IL4','IL13','ARG1','CD40LG','TPH1','XBP1','SOCS1',sep = ",")
df[2,]$gene<-paste('VEGFA','CXCL8','CXCL1','CXCL2','TEK','ANGPT2','CXCL12','FLT1','CXCR4','CTSB','CTSS','SEMA4D','PROK2',
                   'CXCR1','S1PR1','PDGFA','FGF2','MMP12','VEGFC','DNTTIP2','PGF','HIF1A','NFE2L2','CCND2','CCNE1','CD44',
                   'E2F3','EDN1','EZH2','FGF18','FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2','NOTCH1','PTK2','SPP1','STC1',
                   'TNFAIP6','TYMP','VAV2',sep = ",")
df[3,]$gene<-paste('MMP7','MMP2','MMP3','MMP9','CCL2','CCL3','CCL22','TGFB1','EGF','CCR2','FLT1','CTSS','CTSB','CXCR3',
                   'WNT7B','ALOX5','CDH1','CCL18','CCL24','S1PR1','CXCL16','MAPK7','CSF1','SPARC','TLR4','VCAM1','CCL20',sep = ",")
df[4,]$gene<-paste('IL1B','IL6','IL12A','IL23A','TNF','CXCL9','CXCL10','TLR2','TLR4','CXCL11','IFNG','CD40','FCGR2A',
                   'ITGAX','IFNGR1','HLA-DPB1','HLA-DPA1','HLA-DRA','HLA-DRB1','HLA-DQA1','CD74','HLA-DRB5','IRF5',sep = ",")
df[5,]$gene<-paste('MRC1','CD163','MERTK','C1QB','FCRLA','CD5L','CD81','GPNMB','CD36',sep = ",")
df[6,]$gene<-paste('IL23A','TNF','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','CD74','HLA-DPA1','HLA-DPB1',
                   'HLA-DRB5','HLA-DRB1','CD83','CD68','CD80','S100A8','S100A9',sep = ",")
df[7,]$gene<-paste('CCL20','CCL18','CCL22','LYVE1','VEGFA','CTSB','CTSD','TGFB1','MMP19','MMP9','CLEC7A','WNT7B','FASLG',
                   'TNFSF12','TNFSF8','CD276','MSR1','FN1','IRF4','MRC1','CD163','ARG1','IL10','MARCO',sep = ",")
df[8,]$gene<-paste("LAYN","TIGIT","CTLA4","PDCD1","HAVCR2","ENTPD1","EOMES","TOX",sep = ",")
df[9,]$gene<-paste("PRF1","IFNG","GNLY","NKG7","GZMB","GZMH","FGFBP2","CCL4",sep = ",")
df[10,]$gene<-paste("HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M",sep = ",")

head(df)

marker_list = list(
  exhaustion = c("LAYN","TIGIT","CTLA4","PDCD1","HAVCR2","ENTPD1","EOMES","TOX"),
  pro_fibrotic_signature = c( "CCL22", "CSF1", "CHIT1", "FCMR", "SIGLEC15", "CTSK", "COL22A1", "CHI3L1",
                              "SPP1", "SPARC", "MMP9", "MMP7", "GPC4", "TIMP3", "MATK", "LIPA", "PALLD", "FDX1", "LPL"),
  ECM = ECM ,
  Cytolytics_effector = c("PRF1","IFNG","GNLY","NKG7","GZMB","GZMH","FGFBP2","CCL4"),
  
  Phagocytosis = c('MRC1','CD163','MERTK','C1QB','FCRLA','CD5L','CD81','GPNMB','CD36'),
  Anti_inflammatory = c('CCL20','CCL18','CCL22','LYVE1','VEGFA','CTSB','CTSD','TGFB1','MMP19','MMP9','CLEC7A','WNT7B','FASLG',
                        'TNFSF12','TNFSF8','CD276','MSR1','FN1','IRF4','MRC1','CD163','ARG1','IL10','MARCO'),
  
  Immunosuppression = c('CD274','PDCD1LG2','ICOSLG','CD276','VTCN1','VSIR','IDO1','TGFB1','IL10','CCL5','CCL17','CCL22','CXCL8',
                        'CCL16','CCL18','IL17A','CCL24','IL1R1','IL4','IL13','ARG1','CD40LG','TPH1','XBP1','SOCS1'),
  Angiogenesis = c('VEGFA','CXCL8','CXCL1','CXCL2','TEK','ANGPT2','CXCL12','FLT1','CXCR4','CTSB','CTSS','SEMA4D','PROK2',
                   'CXCR1','S1PR1','PDGFA','FGF2','MMP12','VEGFC','DNTTIP2','PGF','HIF1A','NFE2L2','CCND2','CCNE1','CD44',
                   'E2F3','EDN1','EZH2','FGF18','FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2','NOTCH1','PTK2','SPP1','STC1',
                   'TNFAIP6','TYMP','VAV2'),
  Inflammatory = c('IL23A','TNF','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','CD74','HLA-DPA1','HLA-DPB1',
                   'HLA-DRB5','HLA-DRB1','CD83','CD68','CD80','S100A8','S100A9'),
  Tumorsuppression = c('IL1B','IL6','IL12A','IL23A','TNF','CXCL9','CXCL10','TLR2','TLR4','CXCL11','IFNG','CD40','FCGR2A',
                       'ITGAX','IFNGR1','HLA-DPB1','HLA-DPA1','HLA-DRA','HLA-DRB1','HLA-DQA1','CD74','HLA-DRB5','IRF5'),
  pro_metastasis = c('MMP7','MMP2','MMP3','MMP9','CCL2','CCL3','CCL22','TGFB1','EGF','CCR2','FLT1','CTSS','CTSB','CXCR3',
                     'WNT7B','ALOX5','CDH1','CCL18','CCL24','S1PR1','CXCL16','MAPK7','CSF1','SPARC','TLR4','VCAM1','CCL20')
)
sp = AddModuleScore(sp,features = marker_list)
colnames(sp@meta.data)

sp@meta.data = dplyr::rename(sp@meta.data,exhaustion=Cluster1, pro_fibrotic_signature=Cluster2, ECM=Cluster3,
                             Cytolytics_effector=Cluster4, Phagocytosis=Cluster5,Anti_inflammatory=Cluster6 ,
                             Immunosuppression = Cluster7,Angiogenesis = Cluster8,Inflammatory = Cluster9,
                             Tumorsuppression = Cluster10,pro_metastasis = Cluster11)
fea <-names(marker_list)

fea = c("exhaustion", "pro_fibrotic_signature" ,"ECM" , "Immunosuppression" , "Angiogenesis", "Tumorsuppression", "pro_metastasis" )
dominated = subset(sp,development%in%"unknown",invert=T)
dominated = subset(dominated,new_mimer=="other",invert=T)
myd=as.data.frame(scale(dominated@meta.data[,fea]));
myd=cbind(myd,dominated@meta.data[,c('tissue','new_mimer')]);
head(myd)
# myd=myd[,-6];
myd=dplyr::group_by(myd, new_mimer) %>% dplyr::summarise(across(all_of(fea), mean))
myd=as.data.frame(myd);myd=tibble::column_to_rownames(myd,var='new_mimer'); myd=t(as.matrix(myd))
library(pheatmap)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
colors = c(colorRampPalette(colors = rev(getPalette(10))[1:5])(length(bk)/2),
           colorRampPalette(colors = rev(getPalette(10))[6:10])(length(bk)/2))
pdf("mimer_function_in_dominated.pdf",width = 8,height = 8)
# myd[myd< -0.5] = -2
pheatmap(myd, cellwidth=30, cellheight=30, cluster_cols=TRUE, clustering_method='ward.D',
         cluster_rows=F, treeheight_row=0, kmeans_k=NA, border_color="white", 
         scale="row",
         drop_levels=T,
         show_rownames=TRUE, show_colnames=TRUE, 
         color = colors,
         # color=rev(getPalette(10)),
         legend_breaks=seq(-2,2,0.5),breaks = bk,
         fontsize=12, fontsize_row=15, legend=T, fontsize_col=15 #, display_numbers=ifelse(myd > 0.37,"*","")
         )
dev.off()

ME = subset(sp,ME%in%"unknown",invert=T)
ME = subset(ME,new_mimer=="other",invert=T)
myd=as.data.frame(scale(ME@meta.data[,fea]));
myd=cbind(myd,ME@meta.data[,c('tissue','new_mimer')]);
head(myd)
# myd=myd[,-6];
myd=dplyr::group_by(myd, new_mimer) %>% dplyr::summarise(across(all_of(fea), mean))
myd=as.data.frame(myd);myd=tibble::column_to_rownames(myd,var='new_mimer'); myd=t(as.matrix(myd))
library(pheatmap)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
pdf("mimer_function_in_ME.pdf",width = 8,height = 8)
pheatmap(myd, cellwidth=30, cellheight=30, cluster_cols=TRUE, 
         clustering_method='ward.D',
         cluster_rows=F, treeheight_row=0, kmeans_k=NA, border_color="white",
         scale="row",
         drop_levels=T,
         show_rownames=TRUE, show_colnames=TRUE, 
         color = colors,
         legend_breaks=seq(-2,2,0.5),breaks = bk,
         # color=rev(getPalette(10)),
         fontsize=12, fontsize_row=15, legend=T, fontsize_col=15 #, display_numbers=ifelse(myd > 0.37,"*","")
)
dev.off()
# Figure 4g -----
sub.dist.ME_only = list()
# "COL1A1_only" "SPP1_only"   "mimer"       "OX40_only"   "PDCD1_only"  "RGCC_only"  
only = c("COL1A1_only","SPP1_only","OX40_only","PDCD1_only","RGCC_only")
for (ct in only) {
  for (i in unique(sp$file)) {
    
    ima = df$image[df$file == i]
    print(paste0(i," begin!!! "))
    location = GetTissueCoordinates(sp@images[[ima]])
    spot_diameter_fullres = sp@images[[ima]]@scale.factors$fiducial
    location$imagerow <- location$imagerow*(65/spot_diameter_fullres)# converts pixels to microns
    location$imagecol <- location$imagecol*(65/spot_diameter_fullres)
    sub = subset(sp,file == i)
    sub@images = sub@images[ima]
    # SpatialDimPlot(sub,group.by = "new_mimer",images = ima,stroke = NA)
    
    plot.dat = location %>% 
      mutate(cellID = rownames(.)) %>%
      left_join(sub@meta.data,by = "cellID")
    
    
    plot.dat %>%
      dplyr::select(c("imagerow","imagecol","cellID")) %>%
      column_to_rownames("cellID") -> pos
    
    dist_mat <- as.matrix(dist(pos))
    dist_mat <- dist_mat + t(dist_mat)
    diag(dist_mat) <- 0
    
    colnames(dist_mat) <- rownames(pos)
    rownames(dist_mat) <- rownames(pos)
    if (sum(plot.dat$new_mimer %in% ct) >1) {
      dist_filter = dist_mat[rownames(dist_mat) %in% plot.dat[plot.dat$ME != "unknown",]$cellID,
                             colnames(dist_mat) %in% plot.dat[plot.dat$new_mimer %in% ct,]$cellID]
      
      min_df <- apply(dist_filter, 1, function(x) min(x[!(x==0)])) 
      # min_ct = apply(dist_filter, 1, function(x) which.min(x[!(x==0)])) 
      min_df <- as.data.frame(min_df)
      colnames(min_df) <- c("min")
      min_df$id = rownames(min_df)
      sub.dist.ME_only[[paste0(ct,'_',i)]] = min_df
      print(paste0(ct,'_',i," done!!! "))
    }
    
  }
  
}
load("./spatial_pathology/sp_meta_info.Rdata")
sub.dist.all = do.call(rbind,sub.dist.ME)
sub.dist.all$type = 'Mimer'

# "COL1A1_only" "SPP1"    "OX40"   "PDCD1"  "RGCC"  
sub.dist.ME_only_COL1A1 = sub.dist.ME_only[grepl("^COL1A1",names(sub.dist.ME_only))]
sub.dist.all_COL1A1 = do.call(rbind,sub.dist.ME_only_COL1A1)
sub.dist.all_COL1A1$type = 'COL1A1_only'

sub.dist.ME_only_SPP1 = sub.dist.ME_only[grepl("^SPP1",names(sub.dist.ME_only))]
sub.dist.all_SPP1 = do.call(rbind,sub.dist.ME_only_SPP1)
sub.dist.all_SPP1$type = 'SPP1_only'

sub.dist.ME_only_OX40 = sub.dist.ME_only[grepl("^OX40",names(sub.dist.ME_only))]
sub.dist.all_OX40 = do.call(rbind,sub.dist.ME_only_OX40)
sub.dist.all_OX40$type = 'OX40_only'

sub.dist.ME_only_PDCD1 = sub.dist.ME_only[grepl("^PDCD1",names(sub.dist.ME_only))]
sub.dist.all_PDCD1 = do.call(rbind,sub.dist.ME_only_PDCD1)
sub.dist.all_PDCD1$type = 'PDCD1_only'

sub.dist.ME_only_RGCC = sub.dist.ME_only[grepl("^RGCC",names(sub.dist.ME_only))]
sub.dist.all_RGCC = do.call(rbind,sub.dist.ME_only_RGCC)
sub.dist.all_RGCC$type = 'RGCC_only'

only_distinct = rbind(sub.dist.all_COL1A1,sub.dist.all_SPP1,sub.dist.all_OX40,sub.dist.all_PDCD1,sub.dist.all_RGCC,sub.dist.all)
only_distinct$cellID = only_distinct$id
table(table(only_distinct$id))
only_distinct$id = NULL
only_distinct %>% left_join(meta,by = "cellID") -> only_plot_df
library(ggrastr)
library(ggcorrplot)
library(ggstar)
library(ggpubr)
tmp=subset(only_plot_df,min < 20)
table(is.na(tmp$ME))
tmp = tmp[!is.na(tmp$ME),]
setwd("./spatial_pathology/distinct/with_only")
cols_mimer = c("Mimer" = 'brown','other'='grey','OX40_only' = "#1fb1aa",
               "SPP1_only" = "#c2bc8d","PDCD1_only" = "#569c9d",
               "RGCC_only" = "#775095","COL1A1_only" = "#b96c35")
               immune_cols  = c(173:183)
pdf(file  = "The_Distance_20_to_type_pathway_score_with_only_ME.pdf",width = 6,height = 4)
for (i in immune_cols) {
  ggplot(tmp,aes(x=min,y= tmp[[i]] ,color=type) )+
    # geom_point_rast() +
    scale_color_manual(values = cols_mimer )+
    stat_cor(show.legend = FALSE)+
    geom_smooth(method = "gam", formula = y~x,level = 0.99)+
    # facet_grid(~ME)+
    labs(y=colnames(tmp)[i],x="Distance to Spot")+
    theme_classic() -> p0
  print(p0)
}
dev.off()
# Figure 4h -----
#up
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
#bottom
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
# Figure 4i & 4j-----
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
# Figure 4l -----
codex2 = codex
names(codex2@images) = unique(codex2$Sample)
unique_types_cols = rep('grey',length(unique_types))
names(unique_types_cols) = unique_types
Mimer = c("SPP1+Mac",'ISG15+Mac','Ki-67+Mac',
          'CD4T-aTreg','Tcells.cc',
          'PDCD1+Tex','PDCD1+Tex2','LAG3+Tex',
          'COL1A1+CAF','POSTN+CAF',
          'PLVAP+Endo','Vimentin+Endo')
unique_types_cols[names(unique_types_cols) %in% Mimer] = 'red'


unique_types <- unique(Mimer) 

CNS = unique(codex2$neighborhood10_new)[unique(codex2$neighborhood10_new)!='other']

type_dist_results_list <- list()
type_dist_df_CN_list <- list()
type_dist_df_CN_sample_list <- list()
Samples = names(codex2@images) 

for (sa in Samples) {
  location = GetTissueCoordinates(codex2@images[[sa]])
  location %>%
    column_to_rownames("cell") -> pos
  dist_mat <- as.matrix(dist(pos))
  dist_mat <- (dist_mat + t(dist_mat))/2
  diag(dist_mat) <- 0
  colnames(dist_mat) <- rownames(pos)
  rownames(dist_mat) <- rownames(pos)
  distance_matrix_df = dist_mat
  dist_mat[1:4,1:4]
  
  type_combinations <- expand.grid(type1 = unique_types, type2 = unique_types) %>%
    filter(as.numeric(type1) <= as.numeric(type2)) # 确保只处理一次 (A-B 和 B-A 算一次)
  cell_ids = codex2@meta.data %>%
    mutate(cell_id = rownames(.),
           cell_type = subCelltype) %>% 
    filter(Sample==sa) %>%
    dplyr::select(cell_id,cell_type,neighborhood10)
  # 用于存储结果的列表
    cell_id_to_type = cell_ids

    # 遍历所有细胞类型组合
    for (i in 1:nrow(type_combinations)) {
      type1 <- type_combinations$type1[i] %>% as.character()
      type2 <- type_combinations$type2[i] %>% as.character()
      
      
      # 获取属于当前类型的细胞ID
      cells_type1_ids <- cell_id_to_type %>%
        filter(cell_type == type1) %>%
        pull(cell_id)
      
      cells_type2_ids <- cell_id_to_type %>%
        filter(cell_type == type2) %>%
        pull(cell_id)
      
      # 确保有细胞属于这些类型
      if (length(cells_type1_ids) == 0 || length(cells_type2_ids) == 0) {
        next # 跳过没有细胞的类型
      }
      
      # 从距离矩阵中提取所有 Type1 和 Type2 之间两两距离的子集
      # 使用行名和列名进行子集选择
      sub_matrix <- distance_matrix_df[cells_type1_ids, cells_type2_ids, drop = FALSE]
      
      # 将子矩阵展平为向量
      flat_distances <- as.vector(sub_matrix)
      
      # 如果是同一种类型 (Type1 == Type2)，我们需要处理对角线 (自距离为0)
      if (type1 == type2) {
        # 移除值为0的自连接距离
        flat_distances <- flat_distances[flat_distances != 0]
      }
      
      # 计算各种距离统计量
      if (length(flat_distances) > 0) {
        mean_dist <- mean(flat_distances)
        min_dist <- min(flat_distances)
        q25 = quantile(flat_distances,0.25)
        q20 = quantile(flat_distances,0.20)
        max_dist <- max(flat_distances)
      } else {
        mean_dist <- NA
        min_dist <- NA
        max_dist <- NA
        q25 <- NA
        q20 <- NA
      }
      
      # 存储结果
      type_dist_results_list[[i]] <- data.frame(
        Type1 = type1,
        Type2 = type2,
        Mean_Distance = mean_dist,
        Min_Distance = min_dist,
        Max_Distance = max_dist
      )
    }
    
    # 将结果列表合并为数据框
    type_dist_df <- do.call(rbind, type_dist_results_list)
    
    type_dist_df$Sample = sa
  type_dist_df_CN_sample_list[[sa]] = type_dist_df
  print(paste0(i,CN,sa,"done!!!"))
}

type_dist_df_CN_sample_list_df_all  = do.call(rbind,type_dist_df_CN_sample_list)
save(type_dist_df_CN_sample_list_df_all,file = "all_region_distinct.Rdata")
head(type_dist_df_CN_sample_list_df_all)
codex@meta.data %>%
  distinct(Sample,.keep_all = T) %>%
  dplyr::select(Sample,Tissue,Patient_ID) -> meta_tissue
theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   # 图例键为空
    legend.text = element_text(color="black",size=10), # 定义图例文本
    legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
    legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
    legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
    legend.background=element_blank()) 
}
tissue_cols <- c("Adj" = "#377EB8", "Tumor" = "#E41A1C") 
type_dist_df_CN_sample_list_df_all %>%
  group_by(Sample,) %>%
  summarise(mean_distance = mean(Min_Distance)) %>%
  left_join(meta_tissue,by = "Sample") %>%
  # filter(mean_distance <25) %>%
  ggplot(.,aes(x=Tissue,y=mean_distance)) +
  geom_boxplot(aes(fill = Tissue),lwd = 0.3,width=0.5,outliers = F, position = position_dodge(width = 0.2))+
  geom_jitter(aes(color=Tissue),
             width = 0.2)+
  stat_compare_means(comparisons = list(c("Adj","Tumor")),method = 'wilcox.test')+
  scale_fill_manual(values = tissue_cols)+
  scale_color_manual(values = alpha(tissue_cols,0.5))+labs(y='The average minimum distance\n between Mimer members')+
  theme_bw()
ggsave(filename = 'Average_minimum_distance.pdf',width = 3,height = 4)

type_dist_df_CN_sample_list_df_all %>%
  group_by(Sample,) %>%
  summarise(mean_distance = mean(Mean_Distance)) %>%
  left_join(meta_tissue,by = "Sample") %>%
  ggplot(.,aes(x=Tissue,y=mean_distance)) +
  geom_boxplot(aes(fill = Tissue),lwd = 0.3,width=0.5,outliers = F, position = position_dodge(width = 0.2))+
  geom_jitter(aes(color=Tissue),
              width = 0.2)+
  stat_compare_means(comparisons = list(c("Adj","Tumor")),method = 'wilcox.test')+
  scale_fill_manual(values = tissue_cols)+
  scale_color_manual(values = alpha(tissue_cols,0.5))+labs(y='The average mean distance\n between Mimer members')+
  theme_bw()


type_dist_df_CN_sample_list_df_all %>%
  filter(!is.na(Mean_Distance)) %>%
  group_by(Sample) %>%
  summarise(mean_distance = mean(Mean_Distance)) -> mean_min_distance 

mean_min_distance = mean_min_distance %>%
  left_join(meta_tissue,by='Sample') %>%
  filter(Tissue =='Tumor')
survival_info <- readxl::read_xlsx("~/ESCC_codex/CN/data.xlsx",sheet = 2)
survival_info <- survival_info %>%
  filter(Patient_ID %in% mean_min_distance$Patient_ID)%>%
  left_join(mean_min_distance,by='Patient_ID')

cn_var= "mean_distance"
group_var = "mean_distance_gorup"
threshold <- mean(survival_info[[cn_var]], na.rm = TRUE)
survival_info$time = as.numeric(survival_info$`Survival_time (months)`)
survival_info$status = survival_info$`Survival_status (1, dead; 0, alive)`
value <- surv_cutpoint(survival_info, time = "time", event = "status", variables = "mean_distance",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])

survival_info[[group_var]] <- ifelse(survival_info[[cn_var]] > cut_off, "Far", "Near")
survival_info[[group_var]] <- factor(survival_info[[group_var]], levels = c("Near", "Far"))

library(survival)
# 生存分析
head(survival_info)

fit <- survfit(as.formula(paste("Surv(time, status) ~", group_var)), data = survival_info)  
df <- data.frame(
  time = fit$time,
  surv = fit$surv,
  lower = fit$lower,
  upper = fit$upper,
  strata = rep(names(fit$strata), fit$strata)
)  

cox_formula <- as.formula(paste("Surv(time, status) ~", group_var))
cox_fit <- coxph(cox_formula, data = survival_info)
cox_summary <- summary(cox_fit)

hr <- round(cox_summary$coefficients[1, "exp(coef)"], 3)
p_hr <- signif(cox_summary$coefficients[1, "Pr(>|z|)"], 3)
ci_lower <- round(cox_summary$conf.int[1, "lower .95"], 3)
ci_upper <- round(cox_summary$conf.int[1, "upper .95"], 3)

surv_diff <- survdiff(as.formula(paste("Surv(time, status) ~", group_var)), data = survival_info)
p_logrank <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p_logrank_label <- signif(p_logrank, 3)

# 合并标签
hr_label <- paste0(
  "n(Far) = 50","\n",
  "n(Near) = 20","\n",
  "Logrank p = ", p_logrank_label, "\n",
  "p(HR) = ", p_hr , "\n",
  "HR (high/low) = ", hr, " (", ci_lower, "-", ci_upper, ")\n"
)
library(survival)
library(survminer)
p <- ggsurvplot(
  fit, 
  data = survival_info, 
  conf.int.style = 'step',
  pval = F, 
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("red", "#5887ec"),
  legend.title = '',
  risk.table.y.text = TRUE ,
  risk.table.height = 0.25,
  risk.table.title = "", 
  title = paste0("KM Curve for ", cn_var))
p
max_time <- max(survival_info$time, na.rm = TRUE)
p$plot <- p$plot + annotate("text", x = max_time, y = 0.8, label = hr_label, size = 4, hjust=1) + 
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
p

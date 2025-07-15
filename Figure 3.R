
#--------------------------------
# --------Fig3.a
#--------------------------------

# 加载必要的库
library(Seurat)
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(pheatmap)
library(viridis)

setwd('D:/Postdoc/Projects/Single_cell/ZWM')
# ----------------加载数据。
load(file = 'E:/Data_download/Baiao_cluster/spt_TANF_subset.Rdata')

# 创建TLS核心区域子集对象，用于提取高变基因。

core_tls_cells <- rownames(TAN_subset@meta.data)[
  TAN_subset@meta.data$TLS_mark == 1 & 
  !grepl("Periphery|Scattered", TAN_subset@meta.data$TLS_cluster_label)
]
core_tls_subset <- subset(TAN_subset, cells = core_tls_cells)

#----------1. 提取TLS的基因表达矩阵
topn = 200

# 激活包含目标基因的 assay
DefaultAssay(core_tls_subset) <- "spatial"

# 步骤2：筛选高变基因（仅基于核心区域）
core_tls_subset <- FindVariableFeatures(core_tls_subset,
                                      selection.method = "vst",
                                      nfeatures = 2000)

# 提取高变基因列表
core_hvg <- VariableFeatures(core_tls_subset)[1:topn]

# 提取hvg矩阵
hvg_matrix <- GetAssayData(core_tls_subset, slot = "data")[core_hvg, ]

# 1.2 按TLS聚合
tls_labels <- core_tls_subset@meta.data$TLS_cluster_label
tls_feature_matrix <- t(apply(hvg_matrix, 1, function(x) tapply(x, tls_labels, mean)))

# 1.3. Z-score标准化
tls_zscore <- t(scale(t(tls_feature_matrix)))

# 1.4. 聚类
gene_cluster <- hclust(dist(tls_zscore, method = "euclidean"), method = "ward.D2")
tls_cluster <- hclust(dist(t(tls_zscore), method = "euclidean"), method = "ward.D2")

# -------------2. 创建注释数据框
# 2.1 读取包含TLS周边区域的标记的信息
meta.data.new <- fread(file = './SP_Hallmarks.csv') %>% as.data.frame()
meta.data <- meta.data.new[,c(1,grep('HALLMARK',colnames(meta.data.new)))]
row.names(meta.data) <- as.character(meta.data$cellID)
TAN_subset <- AddMetaData(TAN_subset, metadata = meta.data)

# 2.1.2 读取BCR的突变信息
# 读取原始数据并标准化cell_id格式
db_spt_result <- fread('./SPT_spot_vjccdr3_indices_add_mu_freq.csv') %>% 
  as.data.frame() %>%
  mutate(cell_id = gsub("\\.", "_", cell_id))

# 创建包含所有三个基因座的宽表格
spot_indices_wide <- db_spt_result %>%
  filter(locus %in% c("IGH", "IGK", "IGL")) %>%  # 包含所有目标基因座
  pivot_wider(
    names_from = locus,
    values_from = c(unique_seqs, sum_consensus, sorted_consensus, shannon, gini,sorted_mu_freq, avg_mu_freq, max_mu_freq),
    names_glue = "{locus}.{.value}"  # 自动生成列名如 IGH.unique_clone
  ) %>% column_to_rownames(var = "cell_id")

TAN_subset <- AddMetaData(TAN_subset, metadata = spot_indices_wide)

# 2.1.3 读取自定义通路打分
# 读取自定义通路文件并标准化行名
meta_custompath <- fread(file = './HALLMARK_custompath_TLS_meta.txt') %>% 
  as.data.frame() %>% 
  column_to_rownames(var = 'V1') 

# 选择11个关键通路
custom_paths <- c('exhaustion','Cytolytics_effector',
                  'Phagocytosis','Anti_inflammatory',
                  'Angiogenesis','Inflammatory','Tumorsuppression')

# 添加到TAN_subset的metadata
TAN_subset <- AddMetaData(TAN_subset, metadata = meta_custompath[,custom_paths])

meta_df <- TAN_subset@meta.data 
# 将peripheral_label_formatted转为字符型
meta_df$peripheral_label_formatted <- as.character(meta_df$peripheral_label_formatted)

# 2.2 获取热图中的TLS标签（排除Non-TLS）
tls_labels_in_heatmap <- colnames(tls_zscore)

# 创建注释数据框（新增自定义通路列）
annotation_df <- data.frame(
  row.names = tls_labels_in_heatmap,
  patient = character(length(tls_labels_in_heatmap)),
  tissue = character(length(tls_labels_in_heatmap)),
  cluster_size = integer(length(tls_labels_in_heatmap)),
  TLS_9_avg = numeric(length(tls_labels_in_heatmap)),
  TLS_12_avg = numeric(length(tls_labels_in_heatmap)),
  TLS_50_avg = numeric(length(tls_labels_in_heatmap)),

  # BCR指标
  IGL.unique_seqs_avg = numeric(length(tls_labels_in_heatmap)),
  IGK.unique_seqs_avg = numeric(length(tls_labels_in_heatmap)),
  IGH.unique_seqs_avg = numeric(length(tls_labels_in_heatmap)),
  IGH.avg_mu_freq_avg = numeric(length(tls_labels_in_heatmap)),
  IGK.avg_mu_freq_avg = numeric(length(tls_labels_in_heatmap)),
  IGL.avg_mu_freq_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_IGH.unique_seqs_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_IGK.unique_seqs_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_IGL.unique_seqs_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_IGH.avg_mu_freq_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_IGK.avg_mu_freq_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_IGL.avg_mu_freq_avg = numeric(length(tls_labels_in_heatmap)),

  # TLS核心区域指标
  TLS_ALLOGRAFT_REJECTION_avg = numeric(length(tls_labels_in_heatmap)), 

  # === 新增核心区域自定义通路（11个）===
  TLS_exhaustion_avg = numeric(length(tls_labels_in_heatmap)),
  TLS_Cytolytics_effector_avg = numeric(length(tls_labels_in_heatmap)),
  TLS_Phagocytosis_avg = numeric(length(tls_labels_in_heatmap)),
  TLS_Anti_inflammatory_avg = numeric(length(tls_labels_in_heatmap)),
  TLS_Angiogenesis_avg = numeric(length(tls_labels_in_heatmap)),
  TLS_Inflammatory_avg = numeric(length(tls_labels_in_heatmap)),
  TLS_Tumorsuppression_avg = numeric(length(tls_labels_in_heatmap)),
  
  # TLS周边区域指标
  peripheral_ALLOGRAFT_REJECTION_avg = numeric(length(tls_labels_in_heatmap)),
  
  # === 新增周边区域自定义通路（11个）===
  peripheral_exhaustion_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_Cytolytics_effector_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_Phagocytosis_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_Anti_inflammatory_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_Angiogenesis_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_Inflammatory_avg = numeric(length(tls_labels_in_heatmap)),
  peripheral_Tumorsuppression_avg = numeric(length(tls_labels_in_heatmap)),
  stringsAsFactors = FALSE
)

# 2.3 填充注释数据（添加自定义通路计算）
for (tls_label in tls_labels_in_heatmap) {
  # 提取核心区域信息
  core_cells <- meta_df[meta_df$TLS_cluster_label == tls_label, ]

  if (nrow(core_cells) > 0) {
    annotation_df[tls_label, "patient"] <- core_cells$patient[1]
    annotation_df[tls_label, "tissue"] <- core_cells$tissue[1]
    annotation_df[tls_label, "cluster_size"] <- core_cells$cluster_size[1]
    
    # 计算核心区域的三个TLS特征平均值
    annotation_df[tls_label, "TLS_9_avg"] <- mean(core_cells$TLS_9, na.rm = TRUE)
    annotation_df[tls_label, "TLS_12_avg"] <- mean(core_cells$TLS_12, na.rm = TRUE)
    annotation_df[tls_label, "TLS_50_avg"] <- mean(core_cells$TLS_50, na.rm = TRUE)

    # BCR指标（含NA处理）
    annotation_df[tls_label, "IGH.avg_mu_freq_avg"] <- mean(core_cells$IGH.avg_mu_freq, na.rm = TRUE)
    annotation_df[tls_label, "IGK.avg_mu_freq_avg"] <- mean(core_cells$IGK.avg_mu_freq, na.rm = TRUE)
    annotation_df[tls_label, "IGL.avg_mu_freq_avg"] <- mean(core_cells$IGL.avg_mu_freq, na.rm = TRUE)
    annotation_df[tls_label, "TLS_ALLOGRAFT_REJECTION_avg"] <- mean(core_cells$HALLMARK_ALLOGRAFT_REJECTION, na.rm = TRUE)

    # 计算核心区域的自定义通路平均值
    for (path in custom_paths) {
      col_name <- paste0("TLS_", path, "_avg")
      annotation_df[tls_label, col_name] <- mean(core_cells[[path]], na.rm = TRUE)
    }

    # 将BCR指标中的NA替换为0后再计算平均值
    bcr_cells <- core_cells[, c("IGL.unique_seqs", "IGK.unique_seqs", "IGH.unique_seqs")]
    bcr_cells[is.na(bcr_cells)] <- 0  # 将所有NA替换为0
    
    annotation_df[tls_label, "IGL.unique_seqs_avg"] <- mean(bcr_cells$IGL.unique_seqs)
    annotation_df[tls_label, "IGK.unique_seqs_avg"] <- mean(bcr_cells$IGK.unique_seqs)
    annotation_df[tls_label, "IGH.unique_seqs_avg"] <- mean(bcr_cells$IGH.unique_seqs)
  }

  # 计算周边区域的指标
  peripheral_cells <- meta_df[
    meta_df$TLS_cluster_label == "Non-TLS" & 
    meta_df$peripheral_label_formatted == tls_label,
  ]

  if (nrow(peripheral_cells) > 0) {
    # 周边区域ALLOGRAFT_REJECTION
    annotation_df[tls_label, "peripheral_ALLOGRAFT_REJECTION_avg"] <- mean(
      peripheral_cells$HALLMARK_ALLOGRAFT_REJECTION,
      na.rm = TRUE
    )

    # ============ 计算周边区域的BCR指标 ============
    # 对于unique_seqs指标：将NA替换为0后再平均
    peripheral_uniqueseqs <- peripheral_cells[, c("IGH.unique_seqs", "IGK.unique_seqs", "IGL.unique_seqs")]
    peripheral_uniqueseqs[is.na(peripheral_uniqueseqs)] <- 0
    annotation_df[tls_label, "peripheral_IGH.unique_seqs_avg"] <- mean(peripheral_uniqueseqs$IGH.unique_seqs)
    annotation_df[tls_label, "peripheral_IGK.unique_seqs_avg"] <- mean(peripheral_uniqueseqs$IGK.unique_seqs)
    annotation_df[tls_label, "peripheral_IGL.unique_seqs_avg"] <- mean(peripheral_uniqueseqs$IGL.unique_seqs)

    # 对于avg_mu_freq指标：忽略NA直接平均
    annotation_df[tls_label, "peripheral_IGH.avg_mu_freq_avg"] <- mean(peripheral_cells$IGH.avg_mu_freq, na.rm = TRUE)
    annotation_df[tls_label, "peripheral_IGK.avg_mu_freq_avg"] <- mean(peripheral_cells$IGK.avg_mu_freq, na.rm = TRUE)
    annotation_df[tls_label, "peripheral_IGL.avg_mu_freq_avg"] <- mean(peripheral_cells$IGL.avg_mu_freq, na.rm = TRUE)

    # ============ 新增：计算周边区域的自定义通路平均值 ============
    for (path in custom_paths) {
      col_name <- paste0("peripheral_", path, "_avg")
      annotation_df[tls_label, col_name] <- mean(peripheral_cells[[path]], na.rm = TRUE)
    }

  } else {
    annotation_df[tls_label, "peripheral_ALLOGRAFT_REJECTION_avg"] <- NA

    # ============ 新增：周边区域BCR指标为空时设为NA ============
    annotation_df[tls_label, "peripheral_IGH.unique_seqs_avg"] <- NA
    annotation_df[tls_label, "peripheral_IGK.unique_seqs_avg"] <- NA
    annotation_df[tls_label, "peripheral_IGL.unique_seqs_avg"] <- NA
    annotation_df[tls_label, "peripheral_IGH.avg_mu_freq_avg"] <- NA
    annotation_df[tls_label, "peripheral_IGK.avg_mu_freq_avg"] <- NA
    annotation_df[tls_label, "peripheral_IGL.avg_mu_freq_avg"] <- NA

    # ============ 新增：周边区域自定义通路设为NA ============
    for (path in custom_paths) {
      col_name <- paste0("peripheral_", path, "_avg")
      annotation_df[tls_label, col_name] <- NA
    }
  }
}

# 2.4 添加TLS聚类类别
# 对TLS进行ward.D2聚类并切分为3类
tls_cluster_3 <- cutree(tls_cluster, k = 3)
# 将聚类结果添加到注释数据框
annotation_df$TLS_cluster <- factor(tls_cluster_3[rownames(annotation_df)],
                                   levels = 3:1, #顺序按照热图从左到右为G1-G3，所以level设定反过来，3:1,保持与ZXH的标注一致
                                   labels = paste0("TLS_G", 1:3))

# 2.5 保存TLS聚类与列名对应关系表
tls_cluster_mapping <- data.frame(
  TLS_label = rownames(annotation_df),
  TLS_cluster = annotation_df$TLS_cluster
)
# 保存到CSV文件
write.csv(tls_cluster_mapping, 
          file = "TLS_cluster_mapping.csv",
          row.names = FALSE)

#------------------3. 绘图设置（添加自定义通路颜色映射）
# 创建调色板
blue_palette <- colorRampPalette(c("grey90", "#0066CC"))(100)
red_pallette <- colorRampPalette(c("grey90", "darkred"))(100)
lightblue_palette <- colorRampPalette(c("grey90", "#33A1DE"))(100)
lightred_palette <- colorRampPalette(c("grey90", "#FF6666"))(100)
darkgreen_palette <- colorRampPalette(c("grey90", "#006d2c"))(100)
lightgreen_palette <- colorRampPalette(c("grey90", "#66C2A5"))(100)  # 新增浅绿色

# 设置TLS聚类类别的颜色
tls_cluster_colors <- setNames(brewer.pal(3, "Set1"),
                              levels(annotation_df$TLS_cluster))

annotation_colors <- list(
  patient = setNames(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set3"))[seq_len(nlevels(factor(annotation_df$patient)))],
                    levels(factor(annotation_df$patient))),
  tissue = setNames(brewer.pal(8, "Set2")[seq_len(nlevels(factor(annotation_df$tissue)))],
                   levels(factor(annotation_df$tissue))),
  TLS_9_avg = viridis(100),
  TLS_12_avg = viridis(100),
  TLS_50_avg = viridis(100),

  # TLS聚类类别颜色设置
  TLS_cluster = tls_cluster_colors,

  # 6种BCR相关注释
  IGL.unique_seqs_avg = blue_palette,
  IGK.unique_seqs_avg = blue_palette,
  IGH.unique_seqs_avg = blue_palette,
  IGH.avg_mu_freq_avg = red_pallette,
  IGK.avg_mu_freq_avg = red_pallette,
  IGL.avg_mu_freq_avg = red_pallette,

  # 周边区域BCR指标
  peripheral_IGL.unique_seqs_avg = lightblue_palette,
  peripheral_IGK.unique_seqs_avg = lightblue_palette,
  peripheral_IGH.unique_seqs_avg = lightblue_palette,
  peripheral_IGH.avg_mu_freq_avg = lightred_palette,
  peripheral_IGK.avg_mu_freq_avg = lightred_palette,
  peripheral_IGL.avg_mu_freq_avg = lightred_palette,

  # ===== 自定义通路颜色设置 =====
  # 核心区域使用深绿渐变（与ALLOGRAFT一致）
  TLS_ALLOGRAFT_REJECTION_avg = darkgreen_palette,
  TLS_exhaustion_avg = darkgreen_palette,
  TLS_Cytolytics_effector_avg = darkgreen_palette,
  TLS_Phagocytosis_avg = darkgreen_palette,
  TLS_Anti_inflammatory_avg = darkgreen_palette,
  TLS_Angiogenesis_avg = darkgreen_palette,
  TLS_Inflammatory_avg = darkgreen_palette,
  TLS_Tumorsuppression_avg = darkgreen_palette,

  # 周边区域使用浅绿渐变
  peripheral_ALLOGRAFT_REJECTION_avg = lightgreen_palette,
  peripheral_exhaustion_avg = lightgreen_palette,
  peripheral_Cytolytics_effector_avg = lightgreen_palette,
  peripheral_Phagocytosis_avg = lightgreen_palette,
  peripheral_Anti_inflammatory_avg = lightgreen_palette,
  peripheral_Angiogenesis_avg = lightgreen_palette,
  peripheral_Inflammatory_avg = lightgreen_palette,
  peripheral_Tumorsuppression_avg = lightgreen_palette
)

# 调整热图尺寸
options(repr.plot.width = 18, repr.plot.height = topn/4)

#------------------4. 绘制热图
p <- pheatmap(
  mat = tls_zscore,
  cluster_rows = gene_cluster,
  cluster_cols = tls_cluster,
  show_rownames = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  main = "TLS Subtypes by HVG Expression",
  annotation_col = annotation_df,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  annotation_legend = TRUE,
  gaps_col = c(7,49), # 三组区分的位置
  cutree_col = 3,
  silent = TRUE,
)

print(p)

#--------------------------------
# ----------------Fig3.c-d
#--------------------------------
library(ggplot2)
library(ggsignif)
library(rstatix)
library(dplyr)
library(tidyr)

# 获取热图中列（TLS）的排序顺序（从左到右）
column_order <- colnames(tls_zscore)[p$tree_col$order]

allcells <- colnames(core_avg_grouped)[-c(1,ncol(core_avg_grouped))]
all_imm <- allcells[grep(pattern = '^B|^CD4|^CD8|^Mac|^Mast|Neutrophil|^Pla|^T.|^DC|^NK|MAIT',allcells)]

# 1. 核心区数据矩阵
core_mat <- core_avg_grouped %>% 
select(TLS_cluster_label, all_of(all_imm)) %>% 
slice(match(column_order, TLS_cluster_label)) %>% 
column_to_rownames("TLS_cluster_label") %>% 
t()  # 转置：行为细胞类型，列为TLS

# 2. 周边区数据矩阵
peri_mat <- peri_avg_grouped %>% 
select(TLS_id, all_of(all_imm)) %>% 
slice(match(column_order, TLS_id)) %>% 
column_to_rownames("TLS_id") %>% 
t()

# 3. 筛选目标细胞
core_sig_cells <- sig_diff %>% 
  filter(region == "Core", fdr < 0.05) %>% 
  pull(cell_type)

core_sig_cells <- c(core_sig_cells,'Plasma')

#去掉非免疫细胞
core_sig_cells <- core_sig_cells[!core_sig_cells %in% c('Neuron','FB.C4.APOE','FB.C3.COL1A1','Endo.C2.FBLN5','Endo.C3.RGCC','FB.C2.IGF1','FB.C1.CFD','Endo.C4.CCL21','FB.C6.ACTA2')]

peri_sig_cells <- sig_diff %>% 
  filter(region == "Peri", fdr < 0.05) %>% 
  pull(cell_type)
  
#去掉非免疫细胞
peri_sig_cells <- peri_sig_cells[!peri_sig_cells %in% c('Neuron','FB.C4.APOE','FB.C3.COL1A1','Endo.C2.FBLN5','Endo.C3.RGCC','FB.C2.IGF1','FB.C1.CFD','Endo.C4.CCL21','FB.C6.ACTA2')]

# 4. 准备核心区数据
core_long <- core_avg_grouped %>%
  select(TLS_cluster_label, all_of(core_sig_cells)) %>%
  # 添加分组信息
  left_join(rownames_to_column(annotation_df,var = 'TLS_cluster_label') %>% select(TLS_cluster_label, TLS_cluster), by = "TLS_cluster_label") %>%
  # 转换为长格式
  pivot_longer(
    cols = all_of(core_sig_cells),
    names_to = "CellType",
    values_to = "Abundance"
  )

# 5. 准备周边区数据
peri_long <- peri_avg_grouped %>%
  select(TLS_id, all_of(peri_sig_cells)) %>%
  # 添加分组信息
  left_join(rownames_to_column(annotation_df,var = 'TLS_id') %>% select(TLS_id, TLS_cluster), by = "TLS_id") %>%
  # 转换长格式
  pivot_longer(
    cols = all_of(peri_sig_cells),
    names_to = "CellType",
    values_to = "Abundance"
  )
colnames(core_long)[1] <- 'TLS_id'

# 6. 添加区域标识
core_long$Region <- "Core"
peri_long$Region <- "Peri"

# 7. 合并数据
combined_data <- bind_rows(core_long, peri_long)

# 8.差异细胞类型箱线图绘制
# 加载必要的包
library(tidyverse)
library(ggpubr)
library(ggsignif)

# 设置颜色方案和全局参数
cluster_colors <- c('#E41A1C', '#377EB8', '#4DAF4A')  # 红、蓝、绿
options(repr.plot.width = 16, repr.plot.height = 12)  # 增大画布尺寸

# 将数据分为Core和Peri区域
core_data <- combined_data %>% filter(Region == "Core")
peri_data <- combined_data %>% filter(Region == "Peri")

# 设置统计比较组
comparisons <- list(
  c("TLS_G1", "TLS_G2"),
  c("TLS_G1", "TLS_G3"),
  c("TLS_G2", "TLS_G3")
)

# 创建绘图函数（使用星号标记显著性）
create_abundance_plot <- function(data, region_name) {
  # 按CellType计算最大丰度值和标记位置（优化高度计算）
  y_positions <- data %>%
    group_by(CellType) %>%
    summarise(
      max_val = max(Abundance, na.rm = TRUE),
      q3 = quantile(Abundance, 0.75, na.rm = TRUE)  # 使用第三四分位数作为基准
    ) %>%
    mutate(
      base_y = max_val * 1.05,  # 基础高度
      step1 = base_y * 1.10,    # 第一比较组高度
      step2 = base_y * 1.20,    # 第二比较组高度
      step3 = base_y * 1.30     # 第三比较组高度
    )
  
  # 创建箱线图
  p <- ggplot(data, aes(x = TLS_cluster, y = Abundance, fill = TLS_cluster)) +
    geom_boxplot(width = 0.7, alpha = 0.8, outlier.size = 0, outlier.alpha = 0) +
    geom_jitter(width = 0.3, alpha = 0.5, size = 1.5, shape = 21, color = "black") +
    scale_fill_manual(values = cluster_colors) +
    facet_wrap(~ CellType, scales = "free", ncol = 5) +
    labs(
      title = paste(region_name, "Region: CellType Abundance"),
      subtitle = "Significance levels: *** p ≤ 0.001; ** p ≤ 0.01; * p ≤ 0.05",
      x = "TLS Groups",
      y = "Abundance"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold")
    )
  
  # 添加显著性标记（使用星号替代数字P值）
  valid_celltypes <- data %>%
    group_by(CellType) %>%
    summarise(groups = n_distinct(TLS_cluster)) %>%
    filter(groups == 3) %>%
    pull(CellType)
  
  if(length(valid_celltypes) > 0) {
    # 为每个细胞类型添加标记
    for(ct in valid_celltypes) {
      y_pos <- y_positions %>% filter(CellType == ct)
      
      p <- p +
        # G1 vs G2
        geom_signif(
          data = filter(data, CellType == ct),
          comparisons = list(comparisons[[1]]),
          test = "wilcox.test",
          map_signif_level = TRUE,
          textsize = 3,
          tip_length = 0.01,
          y_position = y_pos$step1,
          vjust = 0.5
        ) +
        # G1 vs G3
        geom_signif(
          data = filter(data, CellType == ct),
          comparisons = list(comparisons[[2]]),
          test = "wilcox.test",
          map_signif_level = TRUE,
          textsize = 3,
          tip_length = 0.01,
          y_position = y_pos$step2,
          vjust = 0.5
        ) +
        # G2 vs G3
        geom_signif(
          data = filter(data, CellType == ct),
          comparisons = list(comparisons[[3]]),
          test = "wilcox.test",
          map_signif_level = TRUE,
          textsize = 3,
          tip_length = 0.01,
          y_position = y_pos$step3,
          vjust = 0.5
        )
    }
  }

  return(p)
}

# 创建并显示图形
core_plot <- create_abundance_plot(core_data, "Core")
peri_plot <- create_abundance_plot(peri_data, "Peri")

# 显示图形
print(core_plot)
print(peri_plot)

# # 保存图形
# ggsave("./TLS_core_celltype_abundance_add_plasma.pdf", core_plot,
#        width = 12, height = 8)

# ggsave("./TLS_peri_celltype_abundance.pdf", peri_plot,
#        width = 12, height = 12)

#--------------------------------
# ----------------Fig3.e
#--------------------------------
# 1. 准备数据 ----------------------------------------------------------------
# 提取6种BCR指标数据和TLS聚类类别
bcr_data <- annotation_df[, c("TLS_cluster",
                             "IGL.unique_seqs_avg",
                             "IGK.unique_seqs_avg",
                             "IGH.unique_seqs_avg",

                             "IGL.avg_mu_freq_avg",
                             "IGK.avg_mu_freq_avg",
                             "IGH.avg_mu_freq_avg"),
                             ]

# 2. 转换数据为长格式（便于ggplot绘图）---------------------------------------
library(tidyr)
bcr_long <- pivot_longer(
  bcr_data,
  cols = -TLS_cluster,
  names_to = "BCR_metric",
  values_to = "Value"
)

# 3. 创建更友好的指标名称标签----------------------------------------
metric_labels <- c(
  "IGL.unique_seqs_avg" = "IGL Unique Seqs",
  "IGK.unique_seqs_avg" = "IGK Unique Seqs",
  "IGH.unique_seqs_avg" = "IGH Unique Seqs",

  "IGL.avg_mu_freq_avg" = "IGL Mutation Freq",
  "IGK.avg_mu_freq_avg" = "IGK Mutation Freq",
  "IGH.avg_mu_freq_avg" = "IGH Mutation Freq"
)

bcr_long$Metric_Label <- factor(bcr_long$BCR_metric, 
                               levels = names(metric_labels),
                               labels = metric_labels)

library(ggplot2)
library(ggsignif)
library(rstatix)

# 1. 准备数据（保持不变）
bcr_data <- annotation_df[, c("TLS_cluster", 
                             "IGL.unique_seqs_avg", 
                             "IGK.unique_seqs_avg", 
                             "IGH.unique_seqs_avg",
                             "IGL.avg_mu_freq_avg",
                             "IGK.avg_mu_freq_avg",
                             "IGH.avg_mu_freq_avg")]
bcr_long <- pivot_longer(
  bcr_data,
  cols = -TLS_cluster,
  names_to = "BCR_metric",
  values_to = "Value"
)


# 3. 创建更友好的指标名称标签----------------------------------------
metric_labels <- c(
  "IGL.unique_seqs_avg" = "IGL Unique Seqs",
  "IGK.unique_seqs_avg" = "IGK Unique Seqs",
  "IGH.unique_seqs_avg" = "IGH Unique Seqs",

  "IGL.avg_mu_freq_avg" = "IGL Mutation Freq",
  "IGK.avg_mu_freq_avg" = "IGK Mutation Freq",
  "IGH.avg_mu_freq_avg" = "IGH Mutation Freq"
)

bcr_long$Metric_Label <- factor(bcr_long$BCR_metric, 
                               levels = names(metric_labels),
                               labels = metric_labels)

# 2. 计算各组最大Y值（用于自动调整标记高度）
max_vals <- bcr_long %>% 
  group_by(BCR_metric) %>% 
  summarise(max_val = max(Value, na.rm = TRUE))

# 3. 设置比较组和标记高度
comparisons <- list(
  c("TLS_G1", "TLS_G2"),
  c("TLS_G1", "TLS_G3"),
  c("TLS_G2", "TLS_G3")
)

options(repr.plot.width = 12,repr.plot.height = 8)
# 4. 绘制箱线图并添加显著性标记
cluster_colors <- c('#E41A1C','#377EB8','#4DAF4A')
bcr_boxplot <- ggplot(bcr_long, aes(x = TLS_cluster, y = Value, fill = TLS_cluster)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.size = 0, outlier.alpha = 0) +
  geom_jitter(width = 0.3, alpha = 0.5) +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(~ Metric_Label, scales = "free_y", ncol = 3) +
  labs(title = "BCR Features", x = "TLS Groups", y = "Value") +
  theme_bw(base_size = 12) +
  
  # 关键修复：使用 ggsignif 添加标记 
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",             # 非参数检验
    map_signif_level = TRUE,          # 自动转星号
    tip_length = 0.01,                # 竖线长度
    textsize = 4,                     # 星号大小
    step_increase = 0.1,              # 避免重叠
  )

# 5. 保存结果
print(bcr_boxplot)
ggsave("./TLS_core_BCR_Comparison.pdf", width = 12, height = 8)

#--------------------------------
# ----------------Fig3.f
#--------------------------------
library(ggplot2)
library(ggsignif)  # 关键修复包
library(rstatix)

# 1. 准备数据 ----------------------------------------------------------------
# 已有包含TLS聚类结果的annotation_df
# 提取6种BCR指标数据和TLS聚类类别
bcr_data <- annotation_df[, c("TLS_cluster", 
                             "peripheral_IGL.unique_seqs_avg", 
                             "peripheral_IGK.unique_seqs_avg", 
                             "peripheral_IGH.unique_seqs_avg",

                             "peripheral_IGL.avg_mu_freq_avg",
                             "peripheral_IGK.avg_mu_freq_avg",
                             "peripheral_IGH.avg_mu_freq_avg"),
                             ]

# 2. 转换数据为长格式（便于ggplot绘图）---------------------------------------
library(tidyr)
bcr_long <- pivot_longer(
  bcr_data,
  cols = -TLS_cluster,
  names_to = "BCR_metric",
  values_to = "Value"
)

# 3. 创建更友好的指标名称标签----------------------------------------
metric_labels <- c(
  "peripheral_IGL.unique_seqs_avg" = "peripheral IGL Unique Seqs",
  "peripheral_IGK.unique_seqs_avg" = "peripheral IGK Unique Seqs",
  "peripheral_IGH.unique_seqs_avg" = "peripheral IGH Unique Seqs",

  "peripheral_IGL.avg_mu_freq_avg" = "peripheral IGL Mutation Freq",
  "peripheral_IGK.avg_mu_freq_avg" = "peripheral IGK Mutation Freq",
  "peripheral_IGH.avg_mu_freq_avg" = "peripheral IGH Mutation Freq"
)

bcr_long$Metric_Label <- factor(bcr_long$BCR_metric, 
                               levels = names(metric_labels),
                               labels = metric_labels)


# 2. 计算各组最大Y值（用于自动调整标记高度）
max_vals <- bcr_long %>% 
  group_by(BCR_metric) %>% 
  summarise(max_val = max(Value, na.rm = TRUE))

# 3. 设置比较组和标记高度
comparisons <- list(
  c("TLS_G1", "TLS_G2"),
  c("TLS_G1", "TLS_G3"),
  c("TLS_G2", "TLS_G3")
)

options(repr.plot.width = 12,repr.plot.height = 8)
# 4. 绘制箱线图并添加显著性标记
cluster_colors <- c('#E41A1C','#377EB8','#4DAF4A')
bcr_boxplot <- ggplot(bcr_long, aes(x = TLS_cluster, y = Value, fill = TLS_cluster)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.size = 0, outlier.alpha = 0) +
  geom_jitter(width = 0.3, alpha = 0.5) +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(~ Metric_Label, scales = "free_y", ncol = 3) +
  labs(title = "BCR Features", x = "TLS Groups", y = "Value") +
  theme_bw(base_size = 12) +
  
  # 关键修复：使用 ggsignif 添加标记
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",             # 非参数检验
    map_signif_level = TRUE,          # 自动转星号
    tip_length = 0.01,                # 竖线长度
    textsize = 4,                     # 星号大小
    step_increase = 0.1,              # 避免重叠
  )

# 5. 保存结果
print(bcr_boxplot)

ggsave("./TLS_peri_BCR_Comparison.pdf", width = 12, height = 8)

#--------------------------------
# ----------------Fig3.b
#--------------------------------

# 1. 加载描绘TLS区域所需的函数。
library(ggplot2)
library(dplyr)
library(purrr)
library(RANN)
library(igraph)
library(concaveman)


# 1. 数据预处理函数
prepare_spatial_data <- function(mca, img_name, marker_col = "TLS_50") {
  coords <- GetTissueCoordinates(mca, image = img_name)
  if(is.null(rownames(coords))) rownames(coords) <- colnames(mca)
  
  slice_metadata <- mca@meta.data %>% 
    filter(slice_name == img_name) %>%
    select(cellID = rownames(.), all_of(marker_col))
  
  thresh <- quantile(slice_metadata[[marker_col]], 0.9)
  
  spatial_df <- coords %>%
    as.data.frame() %>%
    mutate(cellID = rownames(.)) %>%
    inner_join(slice_metadata, by = "cellID") %>%
    mutate(TLS_mark = ifelse(.data[[marker_col]] >= thresh, 1, 0))
  
  return(spatial_df)
}

# 2. 计算六边形距离
calculate_hex_distance <- function(coords_df) {
  if (!all(c("imagerow", "imagecol") %in% colnames(coords_df))) {
    stop("数据框必须包含imagerow和imagecol列")
  }
  
  dist_matrix <- as.matrix(dist(coords_df[, c("imagerow", "imagecol")], method = "euclidean"))
  diag(dist_matrix) <- NA
  min_dists <- apply(dist_matrix, 1, min, na.rm = TRUE)
  
  return(median(min_dists))
}

# 3. 获取Visium邻居
get_visium_neighbors <- function(coords_df, hex_distance, radius_multiplier = 1.3) {
  # 参数校验：确保系数在合理范围内 (1.01, 1.7]
  if (radius_multiplier <= 1.01 || radius_multiplier > 1.7) {
    warning("radius_multiplier should be in (1.0, 1.7] for Visium hex grid. Using default 1.3.")
    radius_multiplier <- 1.3
  }
  
  coords <- as.matrix(coords_df[, c("imagerow", "imagecol")])
  nn <- nn2(coords, searchtype = "radius", radius = radius_multiplier * hex_distance)
  
  neighbors <- vector("list", nrow(coords_df))
  for (i in 1:nrow(coords_df)) {
    valid_idx <- nn$nn.idx[i, ]
    valid_dists <- nn$nn.dists[i, ]
    # 严格筛选：仅包含半径内且非自身点
    neighbors[[i]] <- valid_idx[valid_dists > 0 & valid_dists <= radius_multiplier * hex_distance]
  }
  return(neighbors)
}


# 4. 查找连通分量
find_connected_components <- function(spatial_df) {
  required_cols <- c("imagerow", "imagecol", "TLS_mark")
  if (!all(required_cols %in% colnames(spatial_df))) {
    stop("缺少必要列: ", paste(setdiff(required_cols, colnames(spatial_df)), collapse = ", "))
  }
  
  tls_points <- which(spatial_df$TLS_mark == 1)
  if (length(tls_points) == 0) {
    message("没有找到TLS标记点")
    return(list(components = list(), neighbors = list(), hex_dist = NA))
  }
  
  hex_dist <- calculate_hex_distance(spatial_df)
  neighbors <- get_visium_neighbors(spatial_df, hex_dist)
  
  adj_matrix <- matrix(0, nrow = nrow(spatial_df), ncol = nrow(spatial_df))
  for (i in 1:length(neighbors)) {
    if (spatial_df$TLS_mark[i] == 1) {
      valid_nbrs <- neighbors[[i]][spatial_df$TLS_mark[neighbors[[i]]] == 1]
      if(length(valid_nbrs) > 0) {
        adj_matrix[i, valid_nbrs] <- 1
      }
    }
  }
  
  g <- graph.adjacency(adj_matrix, mode = "undirected")
  comp <- components(g)
  tls_indices <- which(spatial_df$TLS_mark == 1)
  components <- split(tls_indices, comp$membership[tls_indices])
  
  return(list(components = components, neighbors = neighbors, hex_dist = hex_dist))
}

# 5. 添加聚类标签
add_cluster_labels <- function(spatial_df, components, neighbors, min_size = 5) {
  spatial_df$TLS_cluster <- 0
  spatial_df$cluster_size <- 0
  spatial_df$is_edge <- FALSE
  
  cluster_sizes <- map_int(components, length)
  large_clusters <- components[cluster_sizes >= min_size]
  small_clusters <- components[cluster_sizes < min_size]
  
  for (i in seq_along(large_clusters)) {
    cluster_points <- large_clusters[[i]]
    spatial_df$TLS_cluster[cluster_points] <- i
    spatial_df$cluster_size[cluster_points] <- length(cluster_points)
  }
  
  for (i in seq_along(small_clusters)) {
    cluster_points <- small_clusters[[i]]
    spatial_df$TLS_cluster[cluster_points] <- -i
    spatial_df$cluster_size[cluster_points] <- length(cluster_points)
  }
  
  for (i in which(spatial_df$TLS_cluster > 0)) {
    nbrs <- neighbors[[i]]
    if (any(spatial_df$TLS_cluster[nbrs] != spatial_df$TLS_cluster[i], na.rm = TRUE)) {
      spatial_df$is_edge[i] <- TRUE
    }
  }
  
  result_summary <- list(
    total_spots = nrow(spatial_df),
    tls_spots = sum(spatial_df$TLS_mark),
    large_clusters = length(large_clusters),
    small_clusters = length(small_clusters),
    min_size = min_size
  )
  
  if (length(large_clusters) > 0) {
    cluster_features <- spatial_df %>%
      filter(TLS_cluster > 0) %>%
      group_by(TLS_cluster) %>%
      summarise(
        size = n(),
        .groups = "drop"
      )
    result_summary$cluster_features <- cluster_features
  }
  
  return(list(spatial_df = spatial_df, summary = result_summary))
}
# 6：标记TLS周边区域
mark_peripheral_regions <- function(spatial_df, neighbors, hex_dist, max_distance = 5) {
  # 初始化周边区域标记列
  spatial_df$peripheral_label <- 0
  spatial_df$peripheral_distance <- NA_real_

  # 获取所有边缘点(edge points)
  edge_points <- which(spatial_df$is_edge & spatial_df$TLS_cluster > 0)

  # 如果没有边缘点，直接返回
  if (length(edge_points) == 0) {
    message("未找到边缘点，无法标记周边区域")
    return(spatial_df)
  }

  # 创建空间坐标矩阵
  coords <- as.matrix(spatial_df[, c("imagecol", "imagerow")])

  # 计算所有点到边缘点的距离
  nn_search <- nn2(
    data = coords[edge_points, , drop = FALSE],
    query = coords,
    k = min(50, length(edge_points)),
    searchtype = "radius",
    radius = max_distance * hex_dist
  )

  # 初始化周边区域数据框
  peripheral_candidates <- data.frame(
    point_idx = integer(),
    cluster = integer(),
    distance = numeric()
  )

  # 收集所有可能的周边点
  for (i in 1:nrow(nn_search$nn.idx)) {
    valid_nbrs <- nn_search$nn.idx[i, ]
    valid_dists <- nn_search$nn.dists[i, ]
    valid_idx <- which(valid_dists > 0 & valid_dists <= max_distance * hex_dist)
    
    if (length(valid_idx) > 0) {
      # 获取对应的边缘点索引和聚类标签
      edge_idx <- edge_points[valid_nbrs[valid_idx]]
      clusters <- spatial_df$TLS_cluster[edge_idx]
      distances <- valid_dists[valid_idx]
      
      # 排除已经是TLS的点
      if (spatial_df$TLS_mark[i] == 0) {
        peripheral_candidates <- rbind(
          peripheral_candidates,
          data.frame(
            point_idx = rep(i, length(valid_idx)),
            cluster = clusters,
            distance = distances
          )
        )
      }
    }
  }

  # 如果没有候选点，直接返回
  if (nrow(peripheral_candidates) == 0) {
    message("未找到符合条件的周边区域点")
    return(spatial_df)
  }

  # 按点和距离分组，保留最近邻的聚类
  peripheral_candidates <- peripheral_candidates %>%
    group_by(point_idx) %>%
    arrange(distance) %>%
    slice(1) %>%
    ungroup()

  # 获取聚类面积信息
  cluster_areas <- spatial_df %>%
    filter(TLS_cluster > 0) %>%
    group_by(TLS_cluster) %>%
    summarise(cluster_area = n(), .groups = "drop")

  # 标记周边区域
  for (i in 1:nrow(peripheral_candidates)) {
    idx <- peripheral_candidates$point_idx[i]
    cluster_id <- peripheral_candidates$cluster[i]

    # 如果点未被标记或当前聚类面积更大
    if (spatial_df$peripheral_label[idx] == 0) {
      spatial_df$peripheral_label[idx] <- cluster_id
      spatial_df$peripheral_distance[idx] <- peripheral_candidates$distance[i]
    } else {
      # 比较聚类面积
      current_area <- cluster_areas$cluster_area[cluster_areas$TLS_cluster == spatial_df$peripheral_label[idx]]
      new_area <- cluster_areas$cluster_area[cluster_areas$TLS_cluster == cluster_id]

      if (new_area > current_area) {
        spatial_df$peripheral_label[idx] <- cluster_id
        spatial_df$peripheral_distance[idx] <- peripheral_candidates$distance[i]
      }
    }
  }
  return(spatial_df)
}

# 7.  旋转坐标
rotate_visium_coords <- function(df, angle) {
  # 备份原始坐标
  df$imagerow_raw <- df$imagerow
  df$imagecol_raw <- df$imagecol

  # 计算组织中心点
  center_x <- (min(df$imagerow) + max(df$imagerow)) / 2
  center_y <- (min(df$imagecol) + max(df$imagecol)) / 2

  # 应用旋转（基于笛卡尔坐标系）
  if (angle == 0) {
    # 无旋转
    df$imagerow <- df$imagerow_raw
    df$imagecol <- df$imagecol_raw
  } else if (angle == 90) {
    # 逆时针90度
    df$imagerow <- center_x + (df$imagecol_raw - center_y)
    df$imagecol <- center_y - (df$imagerow_raw - center_x)
  } else if (angle == 180) {
    # 180度
    df$imagerow <- 2 * center_x - df$imagerow_raw 
    df$imagecol <- df$imagecol_raw - 2 * center_y
  } else if (angle == 270) {
    # 逆时针270度（顺时针90度）
    df$imagerow <- center_x - (df$imagecol_raw - center_y)
    df$imagecol <- center_y + (df$imagerow_raw - center_x)
  }
  return(df)
}

# 8. 修改可视化函数，增加周边区域显示
visualize_tls_clusters_with_periphery <- function(spatial_df, summary,
                                                 periphery_size = 0.8,
                                                 periphery_alpha = 0.4,
                                                 periphery_shape = 16) {
  # 创建新的聚类标签列，按聚类大小排序
  spatial_df <- spatial_df %>%
    mutate(
      temp_label = case_when(
        TLS_cluster > 0 ~ TLS_cluster,
        peripheral_label > 0 ~ peripheral_label + 1000,  # 周边区域使用偏移量
        TLS_cluster < 0 ~ -1,
        TRUE ~ 0
      )
    )

  # 获取大聚类的大小信息并按大小降序排序
  if (summary$large_clusters > 0) {
    # 获取主要聚类信息
    main_clusters <- spatial_df %>%
      filter(TLS_cluster > 0) %>%
      group_by(TLS_cluster) %>%
      summarise(size = mean(cluster_size), .groups = "drop") %>%
      arrange(desc(size)) %>%
      mutate(new_label = paste0("TLS", row_number()))

    # 创建映射关系
    label_mapping <- setNames(
      main_clusters$new_label,
      main_clusters$TLS_cluster
    )

    # 应用新标签
    spatial_df$cluster_class <- case_when(
      spatial_df$temp_label > 0 & spatial_df$temp_label < 1000 ~
        label_mapping[as.character(spatial_df$temp_label)],
      spatial_df$temp_label >= 1000 ~ 
        paste0("Periphery of ", label_mapping[as.character(spatial_df$peripheral_label)]),
      spatial_df$temp_label == -1 ~ "Small cluster",
      TRUE ~ "Not TLS"
    )

    # 设置因子水平
    cluster_levels <- c(
      main_clusters$new_label,
      paste0("Periphery of ", main_clusters$new_label),
      "Small cluster",
      "Not TLS"
    )
  } else {
    spatial_df$cluster_class <- case_when(
      spatial_df$temp_label == -1 ~ "Small cluster",
      TRUE ~ "Not TLS"
    )
    cluster_levels <- c("Small cluster", "Not TLS")
  }

  # 设置因子水平
  spatial_df$cluster_class <- factor(spatial_df$cluster_class, levels = cluster_levels)

  # 创建点类型列
  spatial_df$point_type <- case_when(
    spatial_df$TLS_cluster > 0 & !spatial_df$is_edge ~ "Core",
    spatial_df$TLS_cluster > 0 & spatial_df$is_edge ~ "Edge",
    spatial_df$peripheral_label > 0 ~ "Periphery",
    spatial_df$TLS_cluster < 0 ~ "Small cluster",
    TRUE ~ "Not TLS"
  )

  # 固定颜色设置
  fixed_colors <- c(
    "Not TLS" = "gray90", 
    "Small cluster" = "black"
  )

  # 为大聚类和周边区域生成颜色
  if (summary$large_clusters > 0) {
    # 主聚类颜色
    main_colors <- setNames(
      rainbow(nrow(main_clusters)), 
      main_clusters$new_label
    )

    # 周边区域颜色 - 使用相同颜色但不同形状
    periphery_colors <- setNames(
      main_colors,
      paste0("Periphery of ", main_clusters$new_label)
    )

    color_values <- c(main_colors, periphery_colors, fixed_colors)
  } else {
    color_values <- fixed_colors
  }

  # 创建基础图
  # 首先旋转坐标,和visium对应。
  spatial_df <- rotate_visium_coords(spatial_df,angle = 180)
  p <- ggplot(spatial_df, aes(x = imagecol, y = imagerow)) +
    # 先绘制非TLS点作为背景
    geom_point(
      data = filter(spatial_df, point_type == "Not TLS"),
      aes(color = cluster_class),
      size = 1, alpha = 0.6
    ) +
    # 绘制周边区域点 - 使用三角形形状区分
    geom_point(
      data = filter(spatial_df, point_type == "Periphery"),
      aes(color = cluster_class),
      size = periphery_size, 
      alpha = periphery_alpha, 
      shape = periphery_shape  # 使用三角形形状
    ) +
    # 绘制小聚类点
    geom_point(
      data = filter(spatial_df, point_type == "Small cluster"),
      aes(color = cluster_class),
      size = 1, alpha = 0.8
    ) +
    # 绘制TLS核心和边缘点 - 使用圆形
    geom_point(
      data = filter(spatial_df, point_type %in% c("Core", "Edge")),
      aes(color = cluster_class, size = point_type),
      alpha = 1, shape = 16  # 使用圆形
    ) +
    scale_color_manual(values = color_values) +
    scale_size_manual(values = c(Core = 2.5, Edge = 1.5)) +
    theme_minimal() +
    labs(
      title = paste("TLS Clusters and Periphery in", unique(spatial_df$slice_name)),
      subtitle = sprintf(
        "%d large clusters (≥%d spots), %d peripheral regions", 
        summary$large_clusters, summary$min_size, 
        length(unique(spatial_df$peripheral_label[spatial_df$peripheral_label > 0]))
      ),
      x = "X Position", y = "Y Position",
      color = "Cluster Type", size = "Point Type"
    ) +
    coord_fixed() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )

  # 添加聚类轮廓
  if (summary$large_clusters > 0) {
    for (cl in main_clusters$TLS_cluster) {
      cl_points <- spatial_df %>% filter(TLS_cluster == cl)
      hull <- concaveman(as.matrix(cl_points[, c("imagecol", "imagerow")]))
      p <- p + 
        geom_polygon(
          data = as.data.frame(hull), 
          aes(x = V1, y = V2), 
          fill = NA, color = "black", linetype = "solid", linewidth = 0.6, alpha = 0.7
        )
    }
  }
  return(p)
}

# 2. --------------------------绘制三个TLS_Group和周边区域的映射图
# 2.1-------------加载映射TLS_Group visualize的函数
#  可视化函数 ----
visualize_tls_by_group <- function(spatial_df, summary,
                                   periphery_size = 1,
                                   periphery_alpha = 0.6) {
  # 定义颜色映射

  group_colors <- c(
    "G1" = "#E41A1C",  # 红色
    "G2" = "#377EB8",  # 青色
    "G3" = "#4DAF4A",  # 蓝色
    "Small cluster" = "black",
    "Not TLS" = "gray90"
  )

  spatial_df <- rotate_visium_coords(spatial_df,angle = 180)

  # 创建基础图
  p <- ggplot(spatial_df, aes(x = imagecol, y = imagerow)) +
    # 非TLS背景点
    geom_point(
      data = filter(spatial_df, plot_group == "Not TLS"),
      color = group_colors["Not TLS"],
      size = 0.8, alpha = 0.6, shape = 16
    ) +
    # 小聚类点
    geom_point(
      data = filter(spatial_df, plot_group == "Small cluster"),
      color = group_colors["Small cluster"],
      size = 1, alpha = 0.8, shape = 16
    ) +
    # 周边区域点（圆形）
    geom_point(
      data = filter(spatial_df, point_type == "Periphery"),
      aes(color = plot_group),
      size = periphery_size, 
      alpha = periphery_alpha, 
      shape = 16
    ) +
    # TLS核心和边缘点（圆形）
    geom_point(
      data = filter(spatial_df, point_type %in% c("Core", "Edge")),
      aes(color = plot_group), size = 1,
      alpha = 1, shape = 16
    ) +
    # 设置颜色和大小
    scale_color_manual(values = group_colors) +
    scale_size_manual(values = c(Core = 1, Edge = 1)) +
    # 主题和标签
    theme_minimal() +
    labs(
      title = paste("TLS Clusters by Group in", unique(spatial_df$slice_name)),
      subtitle = sprintf("%d large clusters | %d peripheral regions",
                        summary$large_clusters, 
                        n_distinct(spatial_df$peripheral_label)),
      x = "X Position", y = "Y Position",
      color = "TLS Group", size = "Point Type"
    ) +
    coord_fixed() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # 添加TLS轮廓
  large_clusters <- spatial_df %>%
    filter(TLS_cluster > 0) %>%
    distinct(TLS_cluster, plot_group) %>%
    arrange(TLS_cluster)
  
  for (i in 1:nrow(large_clusters)) {
    cl <- large_clusters$TLS_cluster[i]
    group_val <- large_clusters$plot_group[i]
    
    cl_points <- spatial_df %>% 
      filter(TLS_cluster == cl, point_type %in% c("Core", "Edge"))
    
    if (nrow(cl_points) > 2) {
      hull <- concaveman(as.matrix(cl_points[, c("imagecol", "imagerow")]))
      p <- p + 
        geom_polygon(
          data = as.data.frame(hull), 
          aes(x = V1, y = V2), 
          fill = NA, 
          color = 'black',  # 使用分组颜色
          linetype = "solid", 
          linewidth = 0.5, 
          alpha = 0.8
        )
    }
  }

  return(p)
}

# --------- 2.2 绘制单张slice的TLS_Group
# 主程序
img_name <- "AF01ST000118"

coords <- GetTissueCoordinates(mca, image = img_name)
percentile_90 <- quantile(mca@meta.data[, "TLS_50"], probs = 0.9)[[1]]

slice_metadata <- mca@meta.data %>%
  filter(slice_name == img_name) %>%
  select(cellID, slice_name, TLS_50) %>%
  mutate(TLS_mark = ifelse(TLS_50 < percentile_90, 0, 1))

spatial_df <- slice_metadata %>%
  inner_join(coords %>% mutate(cellID = rownames(.)), by = "cellID")

# 1. 查找连通分量和计算邻居 ----
comp_result <- find_connected_components(spatial_df)

# 2. 添加聚类标签 ----
clust_result <- add_cluster_labels(
  spatial_df, 
  comp_result$components, 
  comp_result$neighbors, 
  min_size = 5  # 最小聚类尺寸
)

# 3. 标记周边区域 ----
spatial_df_with_periphery <- mark_peripheral_regions(
  clust_result$spatial_df,
  comp_result$neighbors,
  comp_result$hex_dist,
  max_distance = 5  # 周边区域最大距离
)

# 4. 添加TLS_ID并与外部数据合并 ----
# 创建TLS_ID列
spatial_df_final <- spatial_df_with_periphery %>%
  mutate(
    TLS_ID = case_when(
      TLS_cluster > 0 ~ paste0(slice_name, "-TLS", TLS_cluster),
      peripheral_label > 0 ~ paste0("peri_", slice_name, "-TLS", peripheral_label),
      TRUE ~ NA_character_
    )
  )

# 5. 与外部分类数据合并
# 5.1 读取TLS_ID和三种类型TLS聚类的对应关系信息。
cor_file_slice_final <- fread(file = 'D:/Postdoc/Projects/Single_cell/ZWM/cor_file_slice_TLS.csv') %>% as.data.frame()
cor_file_slice_final <- cor_file_slice_final[,c('TLS_ID','TLS_Cluster')]
cor_file_slice_final_peri <- cor_file_slice_final
cor_file_slice_final_peri$TLS_ID <- paste0('peri_',cor_file_slice_final_peri$TLS_ID)
cor_file_slice_final_peri$TLS_Cluster <- paste0('peri_',cor_file_slice_final_peri$TLS_Cluster)
cor_file_slice_final <- rbind(cor_file_slice_final,cor_file_slice_final_peri)

# 5.2 合并对应关系准备可视化TLS group
spatial_df_final <- spatial_df_final %>%
  left_join(cor_file_slice_final, by = "TLS_ID") %>%
  mutate(
    # 提取组别信息（G1/G2/G3）
    TLS_Group = str_extract(TLS_Cluster, "G\\d+"),
    # 创建统一的绘图分组
    plot_group = case_when(
      !is.na(TLS_Group) ~ TLS_Group,
      TLS_cluster < 0 ~ "Small cluster",
      TRUE ~ "Not TLS"
    ),
    # 定义点类型
    point_type = case_when(
      TLS_cluster > 0 & !is_edge ~ "Core",
      TLS_cluster > 0 & is_edge ~ "Edge",
      peripheral_label > 0 ~ "Periphery",
      TLS_cluster < 0 ~ "Small cluster",
      TRUE ~ "Not TLS"
    )
  )

# 6. 可视化 ----
final_plot <- visualize_tls_by_group(
  spatial_df_final,
  clust_result$summary
)

# 显示结果
options(repr.plot.width=5,repr.plot.height=5)
print(final_plot)


# 3. -------------------------绘制单张图片单个基因映射图
# 3.1-----加载映射基因的函数
library(RColorBrewer)
my_gradient <- colorRampPalette(c("#094baf",'#ffff31', "#d40a0a"))# "#00abb8", "#7fbc41","#7fbc41",'yellow', "#f39b7f",
colors_custom <- my_gradient(100)

visualize_custom_score <- function(
    spatial_df,
    summary,
    target_column,  # 必选：需展示的特征列名（如"TLS_12"）
    periphery_size = 1,
    periphery_alpha = 0.6,
    color_palette = gradient
) {
  # 校验目标列是否存在且为数值型
  if (!target_column %in% names(spatial_df)) {
    stop("错误：数据框中不存在列 '", target_column, "'")
  }
  if (!is.numeric(spatial_df[[target_column]])) {
    stop("错误：'", target_column, "' 必须是数值型列")
  }

  # 旋转坐标（保持与原始流程一致）
  spatial_df <- rotate_visium_coords(spatial_df, angle = 180)

  # 区域选择逻辑：仅核心/周边区域显示目标特征值（其他区域设为NA）
  spatial_df$display_value <- ifelse(
    spatial_df$TLS_cluster > 0 | spatial_df$peripheral_label > 0,  # 用TLS_50生成的区域标签
    spatial_df[[target_column]],  # 动态引用目标列
    NA
  )

  # 点类型分类（核心/边缘/周边/其他）
  spatial_df$point_type <- case_when(
    spatial_df$TLS_cluster > 0 & !spatial_df$is_edge ~ "Core",
    spatial_df$TLS_cluster > 0 & spatial_df$is_edge ~ "Edge",
    spatial_df$peripheral_label > 0 ~ "Periphery",
    TRUE ~ "Other"
  )

  # 基础绘图
  p <- ggplot(spatial_df, aes(x = imagecol, y = imagerow)) +
    # 绘制其他区域（灰色背景）
    geom_point(
      data = filter(spatial_df, point_type == "Other"),
      color = "gray90", size = periphery_size, alpha = periphery_alpha
    ) +
    # 绘制周边区域
    geom_point(
      data = filter(spatial_df, point_type %in% c("Periphery","Edge","Core")),
      aes(color = display_value),
      size = periphery_size, 
      alpha = periphery_alpha
    ) +
    # 设置颜色梯度（动态范围）
    scale_color_gradientn(
      colors = color_palette,
      na.value = "transparent",
      name = paste0(target_column, ""),
      limits = range(spatial_df[[target_column]], na.rm = TRUE)
    ) +
    theme_void() +
    coord_fixed() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  # 添加大聚类轮廓（基于TLS_50定义的核心区域）
  if (summary$large_clusters > 0) {
    main_clusters <- spatial_df %>%
      filter(TLS_cluster > 0) %>%
      group_by(TLS_cluster) %>%
      summarise(size = n(), .groups = "drop") %>%
      filter(size >= summary$min_size) %>%
      arrange(desc(size))
    for (cl in main_clusters$TLS_cluster) {
      cl_points <- spatial_df %>% filter(TLS_cluster == cl)
      hull <- concaveman::concaveman(as.matrix(cl_points[, c("imagecol", "imagerow")]))
      p <- p + 
        geom_polygon(
          data = as.data.frame(hull), 
          aes(x = V1, y = V2), 
          fill = NA, color = "black", 
          linetype = "solid", linewidth = 0.5, alpha = 0.7
        )
    }
  }
  return(p)
}

#---------------- 3.2 绘制单张基因slie的基因表达映射
img_name <- "TF02ST021022"
coords <- GetTissueCoordinates(mca, image = img_name)
percentile_90 <- quantile(mca@meta.data[, "TLS_50"], probs = 0.9)[[1]]
target_feature <- "TNFSF13"  # 可替换为任意基因名（如"CD3E"）
slice_metadata <- mca@meta.data %>%
  filter(slice_name == img_name) %>%
  select(cellID, slice_name, TLS_50, all_of(target_feature)) %>%  # 动态包含目标列
  mutate(TLS_mark = ifelse(TLS_50 < percentile_90, 0, 1))

spatial_df <- slice_metadata %>%
  inner_join(coords %>% mutate(cellID = rownames(.)), by = "cellID")

# 后续步骤保持不变
comp_result <- find_connected_components(spatial_df)
clust_result <- add_cluster_labels(spatial_df, comp_result$components, comp_result$neighbors, min_size = 5)
spatial_df_with_periphery <- mark_peripheral_regions(
  clust_result$spatial_df, comp_result$neighbors, comp_result$hex_dist, max_distance = 5
)

# 绘图
custom_plot <- visualize_custom_score(
  spatial_df_with_periphery,
  clust_result$summary,
  target_column = target_feature,  # 动态指定展示列
  color_palette =  colors_custom # 自定义色阶
)

# 输出图形
options(repr.plot.width = 8, repr.plot.height = 7)
print(custom_plot)


# 4.-------------- 绘制映射TLS-flare图
# 4.1------------加载数据：
load("E:/Data_download/Baiao_cluster/sp.RData")
meta.data.new<-fread(file = 'E:/Data_download/Cluster_data/SP_Hallmarks.csv') %>% as.data.frame()
TLS_meta<-fread(file = 'D:/Postdoc/Projects/Single_cell/ZWM/TLS_12_9_50_meta.txt') %>% as.data.frame()
meta.data.bak <- mca@meta.data
mca@meta.data <- TLS_meta
row.names(mca@meta.data) <- mca@meta.data$cellID

# 为 mca@meta.data 添加 slice_name 信息
slice_names <- c()

# 遍历 mca@images 并为每个 spot 添加对应的 slice_name
for (slice_index in seq_along(mca@images)) {
  slice_name <- names(mca@images)[slice_index]
  slice_spots <- rownames(mca@images[[slice_index]]@coordinates)
  
  # 为这些 spots 添加对应的 slice_name
  slice_names <- c(slice_names, setNames(rep(slice_name, length(slice_spots)), slice_spots))
}

# 将 slice_name 信息添加到 mca@meta.data
mca@meta.data$slice_name <- slice_names[rownames(mca@meta.data)]

mca@meta.data$TLS_mark <- ifelse(mca@meta.data$TLS == "unknown", 0, 1)
image_name = mca@meta.data$slice_name[mca$TLS_mark==1][1]


#4.2.------------------- 读取根据筛选差异基因，并且完成利用GSVA score比较了in-house治疗队列responder之后筛选的差异基因表
G3_marker_data <- fread(file = './Genelist/Genes_core_4groups_loop10_pct0.25_fc0.5_0.6perc_top500_GSVA_boxplot.csv') %>% as.data.frame()
peri_TLS_G3_genelist <- G3_marker_data$gene[G3_marker_data$cluster == 'TLS_G3']
DefaultAssay(mca) <- "spatial"
mca<-AddModuleScore(mca,features = list(peri_TLS_G3_genelist),name = "TLS_G3_score")
colnames(mca@meta.data)[ncol(mca@meta.data)] <- 'TLS_G3_score'

# 4.3------------------------加载visualize addmodule的映射函数
library(RColorBrewer)
my_gradient <- colorRampPalette(c("#094baf",'#ffff31', "#d40a0a"))# "#00abb8", "#7fbc41","#7fbc41",'yellow', "#f39b7f",
colors_custom <- my_gradient(100)

visualize_custom_score <- function(
    spatial_df, 
    summary,
    target_column,  # 必选：需展示的特征列名（如"TLS_12"）
    periphery_size =1,
    periphery_alpha = 1,
    color_palette = gradient
) {
  # 校验目标列是否存在且为数值型
  if (!target_column %in% names(spatial_df)) {
    stop("错误：数据框中不存在列 '", target_column, "'")
  }
  if (!is.numeric(spatial_df[[target_column]])) {
    stop("错误：'", target_column, "' 必须是数值型列")
  }
  
  # 旋转坐标（保持与原始流程一致）
  spatial_df <- rotate_visium_coords(spatial_df, angle = 180)
  
  # 区域选择逻辑：仅核心/周边区域显示目标特征值（其他区域设为NA）
  spatial_df$display_value <- ifelse(
    spatial_df$TLS_cluster > 0 | spatial_df$peripheral_label > 0,  # 用TLS_50生成的区域标签
    spatial_df[[target_column]],  # 动态引用目标列
    NA
  )
  
  # 点类型分类（核心/边缘/周边/其他）
  spatial_df$point_type <- case_when(
    spatial_df$TLS_cluster > 0 & !spatial_df$is_edge ~ "Core",
    spatial_df$TLS_cluster > 0 & spatial_df$is_edge ~ "Edge",
    spatial_df$peripheral_label > 0 ~ "Periphery",
    TRUE ~ "Other"
  )
  
  # 基础绘图
  p <- ggplot(spatial_df, aes(x = imagecol, y = imagerow)) +
    # 绘制其他区域（灰色背景）
    geom_point(
      data = filter(spatial_df, point_type == "Other"),
      color = "gray90", size = periphery_size, alpha = periphery_alpha
    ) +
    # 绘制周边区域
    geom_point(
      data = filter(spatial_df, point_type %in% c("Periphery","Edge","Core")),
      aes(color = display_value),
      size = periphery_size, 
      alpha = periphery_alpha
    ) +
    # 设置颜色梯度（动态范围）
    scale_color_gradientn(
      colors = color_palette,
      na.value = "transparent",
      name = paste0(target_column, ""),
      limits = range(spatial_df[[target_column]], na.rm = TRUE)
    ) +
    theme_void() +
    coord_fixed() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  # 添加大聚类轮廓（基于TLS_50定义的核心区域）
  if (summary$large_clusters > 0) {
    main_clusters <- spatial_df %>%
      filter(TLS_cluster > 0) %>%
      group_by(TLS_cluster) %>%
      summarise(size = n(), .groups = "drop") %>%
      filter(size >= summary$min_size) %>%
      arrange(desc(size))
    
    for (cl in main_clusters$TLS_cluster) {
      cl_points <- spatial_df %>% filter(TLS_cluster == cl)
      hull <- concaveman::concaveman(as.matrix(cl_points[, c("imagecol", "imagerow")]))
      p <- p + 
        geom_polygon(
          data = as.data.frame(hull), 
          aes(x = V1, y = V2), 
          fill = NA, color = "black", 
          linetype = "solid", linewidth = 0.5, alpha = 0.7
        )
    }
  }
  
  return(p)
}

# 4.4------------------------ 绘制单张slice的TLS_G3 (TLS-flare)的映射
# 主程序数据准备（确保包含目标列）
img_name <- "TF02ST021022"
coords <- GetTissueCoordinates(mca, image = img_name)
percentile_90 <- quantile(mca@meta.data[, "TLS_50"], probs = 0.9)[[1]]

# 关键修改：同时选择TLS_50和目标特征列（如TLS_12）
target_feature <- "TLS_G3_score" 
slice_metadata <- mca@meta.data %>%
  filter(slice_name == img_name) %>%
  select(cellID, slice_name, TLS_50, all_of(target_feature)) %>%  # 动态包含目标列
  mutate(TLS_mark = ifelse(TLS_50 < percentile_90, 0, 1))

spatial_df <- slice_metadata %>%
  inner_join(coords %>% mutate(cellID = rownames(.)), by = "cellID")

# 后续步骤保持不变
comp_result <- find_connected_components(spatial_df)
clust_result <- add_cluster_labels(spatial_df, comp_result$components, comp_result$neighbors, min_size = 5)
spatial_df_with_periphery <- mark_peripheral_regions(
  clust_result$spatial_df, comp_result$neighbors, comp_result$hex_dist, max_distance = 5
)

# 使用修改后的函数绘图
custom_plot <- visualize_custom_score(
  spatial_df_with_periphery,
  clust_result$summary,
  target_column = target_feature,  # 动态指定展示列
  color_palette =  colors_custom # 自定义色阶
)

# 输出图形
options(repr.plot.width = 8, repr.plot.height = 7)
print(custom_plot)


setwd('D:/Postdoc/Projects/Single_cell/ZWM')
#--------------------------------
#------Fig3.h
#------------------------------

# 加载必要包 ----------------------------------------------------------------
library(GSVA)
library(limma)
library(BiocParallel)
library(ggpubr)
library(patchwork)
library(data.table)
library(pROC)
library(DESeq2)

# 1. 数据准备与分组转换 ---------------------------------------------------
combined_data <- fread('./ImmTherapy_cohort/Processed_immcohort_therapy_response_survival_gene_expr.csv')

# 创建二分类分组变量：响应组(CR+PR) vs 非响应组(SD+PD)
combined_data$ResponseGroup <- ifelse(
  combined_data$Best_overall_response %in% c("CR", "PR"), "Responder",
  "Non-Responder"  # SD和PD合并为非响应组
)
combined_data <- combined_data[!is.na(ResponseGroup), ]
combined_data$ResponseGroup <- factor(
  combined_data$ResponseGroup,
  levels = c("Responder", "Non-Responder")
)

# 提取纯数值表达矩阵
gene_start_col <- which(colnames(combined_data) == "TSPAN6")
gene_end_col <- which(colnames(combined_data) == "ResponseGroup") - 1
expr_matrix <- as.matrix(combined_data[, gene_start_col:gene_end_col, with = FALSE])
mode(expr_matrix) <- "integer"
rownames(expr_matrix) <- combined_data$ID

# 分组向量
group <- combined_data$ResponseGroup

# 验证数据结构
stopifnot(
  is.numeric(expr_matrix[1, 1]),
  nrow(expr_matrix) == length(group),
  all(!is.na(group))
)
message(paste("有效样本数:", nrow(expr_matrix), "| 基因数:", ncol(expr_matrix)))

# # 2. 使用DESeq2进行差异分析（两组比较）---------------------------------------
# colData <- data.frame(
#   row.names = rownames(expr_matrix),
#   condition = group
# )

# dds <- DESeqDataSetFromMatrix(
#   countData = t(expr_matrix),
#   colData = colData,
#   design = ~ condition
# )

# # 执行Wald检验（两组比较）
# dds <- DESeq(dds)
# saveRDS(dds,file = './Genelist/Immcohort_deseq2_dds.RDS')
dds <- readRDS(file = './Genelist/Immcohort_deseq2_dds.RDS')
res <- results(dds, alpha = 0.1)  # 放宽FDR阈值
# diff_genes <- rownames(subset(res, padj < 0.1 & abs(log2FoldChange) > 0.5))
diff_genes <- rownames(subset(res, pvalue < 0.2))

# 3. GSVA分析 ---------------------------------------------------------------
#markers <- fread('./Genelist/TLS_Group_marker_7groups_0.1_fc0.5_all_non_TLS.csv')
markers <- top_markers

# markers <- combined_markers

# clusters <- c('Non-TLS', 'peri_TLS_G3', 'TLS_G3', 'peri_TLS_G2', 'TLS_G2', 'peri_TLS_G1', 'TLS_G1')
clusters <- c('Non-TLS', 'TLS_G3',  'TLS_G2', 'TLS_G1')

# 初始化存储结果
plot_list <- list()
auc_results <- data.frame()
validate_markers <- data.frame()

for (cluster in clusters) {
  gene_list <- unique(markers$gene[markers$cluster == cluster])
  valid_genes <- intersect(gene_list, diff_genes)
  
  # 保存经过Deseq筛选的marker基因
  validate_markers <- rbind(validate_markers, markers[markers$cluster == cluster & markers$gene %in% valid_genes,])
  
  if (length(valid_genes) < 3) {
    message(paste0("跳过 ", cluster, "：有效基因数不足 (", length(valid_genes), ")"))
    next
  }
  
  # 使用gsvaParam API
  gsva_param <- gsvaParam(
    exprData = t(expr_matrix),  # 全部基因而不是差异基因
    geneSets = list(cluster = valid_genes),
    kcdf = "Poisson"  # Count数据用Poisson
  )
  gsva_scores <- gsva(gsva_param, verbose = FALSE, BPPARAM = SerialParam())
  scores <- gsva_scores[1, ]
  
  # 创建绘图数据
  plot_data <- data.frame(
    Sample = names(scores),
    Score = scores,
    Group = group,
    Cluster = cluster
  )
  
  # 创建箱线图（显示有效基因数）
  p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = c("#00BA38", "#F8766D")) +  # 两组颜色
    labs(title = paste0(cluster, "\n(", length(valid_genes), " genes)"),
         y = "GSVA Score",
         x = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size=10),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1,size=10)) +
    stat_compare_means(
      method = "wilcox.test",  # 两组比较使用Wilcoxon检验
      label = "p.format",
      tip.length = 0.02,
      size = 5
    )
  
  plot_list[[cluster]] <- p
  
  # 计算AUC（二分类）
  binary_scores <- plot_data$Score
  binary_labels <- ifelse(plot_data$Group == "Responder", 1, 0)
  
  roc_obj <- roc(response = binary_labels, predictor = binary_scores, quiet = TRUE)
  auc_value <- auc(roc_obj)
  
  auc_results <- rbind(auc_results, data.frame(
    Cluster = cluster,
    AUC = round(auc_value, 3),
    Genes = length(valid_genes)
  ))
}

# 4. 图形拼接与输出 -------------------------------------------------------
if (length(plot_list) > 0) {
  # 设置拼接布局
  combined_plot <- wrap_plots(plot_list, ncol = 4) + 
    plot_annotation(title = "GSVA score(Responder vs Non-Responder)",
                    subtitle = paste("DESeq2 filtered genes (p<0.2), total", length(valid_genes), "genes"),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                                  plot.subtitle = element_text(hjust = 0.5, size = 12)))
  
  # 添加共享图例
  legend_plot <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#00BA38", "#F8766D"),
                      name = "Response Group",
                      labels = c("Responder (CR+PR)", "Non-Responder (SD+PD)")) +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold"))
  
  legend <- get_legend(legend_plot)
  
  # 最终拼接
  final_plot <- combined_plot / legend + 
    plot_layout(heights = c(10, 1))
  
  # 显示图形
  options(repr.plot.width = 14, repr.plot.height = 6)
  print(final_plot)
  
  # 保存结果
  #   ggsave("GSVA_Analysis_Plot_TwoGroups.pdf", final_plot, width = 14, height = 16)
  #   fwrite(auc_results, "GSVA_AUC_Results_TwoGroups.csv")
} else {
  warning("没有足够的有效基因集生成图表，请检查差异基因筛选和基因集匹配")
}

options(repr.plot.width = 3, repr.plot.height = 5)

ggplot(auc_results, aes(x='AUC',y=Cluster, fill=AUC)) +
  geom_tile() + 
  geom_text(
    aes(label = sprintf("%.2f", AUC)),  # 显示2位小数的AUC值
    color = "black", 
    size = 8
  ) +
  scale_fill_gradient2(low="white", high="red", midpoint=0.5)+
  theme(axis.text = element_text(size=12,),
        axis.text.x = element_text(angle = 45,hjust = 1)
  )



# fwrite(validate_markers,file=  "./Genelist/Genes_core_4groups_5loop_pct0.25_fc0.5_0.6perc_GSVA_boxplot.csv")

# for (group in clusters) {
#     writeLines(as.character(validate_markers$gene[validate_markers$cluster == group]), paste0("./Genelist/",group,"_peri4groups_5loop_pct0.25_fc0.5_0.6perc.txt"),sep = ',')
# }

#-----------------------
# ----- Fig.3 i
#-----------------------

# 加载必要包 ----------------------------------------------------------------
library(GSVA)
library(limma)
library(BiocParallel)
library(ggpubr)
library(patchwork)
library(data.table)
library(pROC)
library(DESeq2)
library(ggplot2)

# 1. 数据准备与分组转换 ---------------------------------------------------
combined_data <- fread('./ImmTherapy_cohort/Processed_immcohort_therapy_response_survival_gene_expr.csv')

# 创建二分类分组变量：响应组(CR+PR) vs 非响应组(SD+PD)
combined_data$ResponseGroup <- ifelse(
  combined_data$Best_overall_response %in% c("CR", "PR"), "Responder",
  "Non-Responder"  # SD和PD合并为非响应组
)
combined_data <- combined_data[!is.na(ResponseGroup), ]
combined_data$ResponseGroup <- factor(
  combined_data$ResponseGroup,
  levels = c("Responder", "Non-Responder")
)

# 提取纯数值表达矩阵
gene_start_col <- which(colnames(combined_data) == "TSPAN6")
gene_end_col <- which(colnames(combined_data) == "ResponseGroup") - 1
expr_matrix <- as.matrix(combined_data[, gene_start_col:gene_end_col, with = FALSE])
mode(expr_matrix) <- "integer"
rownames(expr_matrix) <- combined_data$ID

# 分组向量
group <- combined_data$ResponseGroup

# 验证数据结构
stopifnot(
  is.numeric(expr_matrix[1, 1]),
  nrow(expr_matrix) == length(group),
  all(!is.na(group))
)
message(paste("有效样本数:", nrow(expr_matrix), "| 基因数:", ncol(expr_matrix)))

# 2. 使用DESeq2进行差异分析（两组比较）---------------------------------------
# colData <- data.frame(
#   row.names = rownames(expr_matrix),
#   condition = group
# )

# dds <- DESeqDataSetFromMatrix(
#   countData = t(expr_matrix),
#   colData = colData,
#   design = ~ condition
# )

# # 执行Wald检验（两组比较）
# dds <- DESeq(dds)
dds <- readRDS(file = './Genelist/Immcohort_deseq2_dds.RDS')
res <- results(dds, alpha = 0.1)  # 放宽FDR阈值
diff_genes <- rownames(subset(res, pvalue < 0.2))

# 3. GSVA分析与ROC曲线绘制--------------------------------------------------
markers <- top_markers
# clusters <- c('Non-TLS', 'peri_TLS_G3', 'TLS_G3', 'peri_TLS_G2', 'TLS_G2', 'peri_TLS_G1', 'TLS_G1')
clusters <- c('Non-TLS', 'TLS_G3',  'TLS_G2', 'TLS_G1')

# 初始化存储结果
plot_list <- list()
auc_results <- data.frame()

for (cluster in clusters) {
  gene_list <- unique(markers$gene[markers$cluster == cluster])
  valid_genes <- intersect(gene_list, diff_genes)
  
  if (length(valid_genes) < 3) {
    message(paste0("跳过 ", cluster, "：有效基因数不足 (", length(valid_genes), ")"))
    next
  }
  
  # 使用gsvaParam API
  gsva_param <- gsvaParam(
    exprData = t(expr_matrix),  # 仅差异基因
    geneSets = list(cluster = valid_genes),
    kcdf = "Poisson"  # Count数据用Poisson
  )
  gsva_scores <- gsva(gsva_param, verbose = FALSE, BPPARAM = SerialParam())
  scores <- gsva_scores[1, ]
  
  # 准备ROC曲线数据
  binary_labels <- ifelse(group == "Responder", 1, 0)
  roc_obj <- roc(response = binary_labels, predictor = scores, quiet = TRUE)
  auc_value <- auc(roc_obj)
  
  # 创建ROC曲线数据框
  roc_data <- data.frame(
    Specificity = roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities
  )
  
  p <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_line(color = "steelblue", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(title = paste0(cluster, "\n(n = ", length(valid_genes), ")"),
         x = "1 - Specificity (False Positive Rate)",
         y = "Sensitivity (True Positive Rate)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    ) +
    # 修复：使用annotate替代geom_text避免惰性求值问题
    annotate("text", 
             x = 0.75, y = 0.25, 
             label = paste0("AUC = ", round(auc_value, 3)),
             size = 5, 
             color = "red", 
             fontface = "bold")
  
  plot_list[[cluster]] <- p
  
  auc_results <- rbind(auc_results, data.frame(
    Cluster = cluster,
    AUC = round(auc_value, 3),
    Genes = length(valid_genes)
  ))
}

# 4. 整合所有ROC曲线到一张大图 ---------------------------------------------
# 计算布局：4x2网格
combined_roc_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(
    title = "ROC Curves of TLS Cluster Signatures",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# 保存组合图
options(repr.plot.width = 16, repr.plot.height = 4)
pdf(file = "./Genelist/core_loop10_cobined_ROC_Curves.pdf", width = 16, height = 4)
print(combined_roc_plot)
dev.off()

# 5. 输出AUC结果
# fwrite(auc_results, "./Genelist/core_loop10_AUC_Results.csv")


#-----------------------
# ----- Fig.3 j
#-----------------------
# 加载必要包 ----------------------------------------------------------------
library(survival)
library(survminer)
library(patchwork)

validate_markers <- fread( "./Genelist/Genes_core_4groups_loop10_pct0.25_fc0.5_0.6perc_top500_GSVA_boxplot.csv")
# 设置生存数据类型 (PFS 或 OS) 和分组方法 --------------------------------
surv_type <- "OS"                  # 可更改为 "PFS"
# group_method <- "optimal"          # "median" 或 "optimal" (基于surv_cutpoint的最佳截断值)
group_method <- "median"          # "median" 或 "optimal" (基于surv_cutpoint的最佳截断值)

# 循环所有基因集生成生存曲线 -----------------------------------------------
survplot_list <- list()
survtable_list <- list()
for (cluster in clusters) {
  # 提取当前基因集
  gene_list <- unique(validate_markers$gene[validate_markers$cluster == cluster]) # combined_markers
  gene_sets <- list()
  gene_sets[[cluster]] <- gene_list
  if (length(gene_list) == 0) next
  
  # GSVA分析
  params <- gsvaParam(
    exprData = t(expr_matrix),
    geneSets = gene_sets,
    kcdf = "Poisson",
    minSize = 1
  )
  gsva_scores <- gsva(params, verbose = FALSE, BPPARAM = SerialParam())
  scores <- gsva_scores[1, ]
  
  # 准备生存数据
  if (surv_type == "OS") {
    surv_data <- data.frame(
      ID = combined_data$ID,
      GSVA_Score = scores,
      Status = combined_data$OS_months_state,
      Time = combined_data$OS
    )
    xlab_text <- "OS Days"
  } else {
    surv_data <- data.frame(
      ID = combined_data$ID,
      GSVA_Score = scores,
      Status = combined_data$OS_months_state,
      Time = combined_data$PFS
    )
    xlab_text <- "PFS Days"
  }
  
  # 分组逻辑 ------------------------------------------------------------
  if (group_method == "median") {
    median_score <- median(surv_data$GSVA_Score)
    surv_data$GSVA_Group <- ifelse(
      surv_data$GSVA_Score >= median_score,
      "High GSVA",
      "Low GSVA"
    )
    group_title <- paste(cluster, "\n(Median Split)")
    
  } else if (group_method == "optimal") {
    cutpoint <- surv_cutpoint(
      data = surv_data,
      time = "Time",
      event = "Status",
      variables = "GSVA_Score",
      minprop = 0.1
    )
    best_cutoff <- as.numeric(cutpoint$cutpoint[1])
    surv_data$GSVA_Group <- ifelse(
      surv_data$GSVA_Score >= best_cutoff,
      "High GSVA",
      "Low GSVA"
    )
    group_title <- paste(cluster, "\n(Optimal Cutoff)")
  }
  
  # 拟合生存曲线
  fit <- survfit(
    Surv(Time, Status) ~ GSVA_Group,
    data = surv_data
  )
  
  # 绘制生存曲线
  surv_plot <- ggsurvplot(
    fit,
    data = surv_data,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    palette = c("#E74C3C", "#3498DB"),
    legend.labs = c("High GSVA", "Low GSVA"),
    legend.title = "GSVA Group",
    xlab = xlab_text,
    ylab = paste0(surv_type),
    ggtheme = theme_bw(),
    title = group_title,  # 动态显示分组方法
    font.title = c(11, "bold"),
    censor = TRUE,               # 显示删失点
    censor.size = 3,             # 标记点大小（默认3）
    censor.shape = "|",          # 标记形状（默认为"+"）
    censor.colour = "black"      # 标记颜色（默认与曲线一致）
  )
  survplot_list[[cluster]] <- surv_plot$plot
  survtable_list[[cluster]] <- surv_plot$table
}

# 图形拼接与输出 ---------------------------------------------------------
# 拼接生存曲线 (3x3布局)
surv_combined <- wrap_plots(survplot_list, ncol = 4) + 
  plot_annotation(
    title = paste0("TLS GeneSets", 
                   ifelse(surv_type == "OS", "OS", "PFS"),
                   " (", ifelse(group_method == "median", "median", "bestcut"), ")"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# 添加共享图例
legend <- get_legend(
  survplot_list[[1]] + 
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1))
)

final_surv_plot <- surv_combined / legend + 
  plot_layout(heights = c(10, 1))

# 显示并保存
options(repr.plot.width = 16, repr.plot.height = 4)
print(final_surv_plot)
# ggsave(
#   paste0("TLS_Survival_", surv_type, "_", group_method, ".png"), 
#   final_surv_plot, width = 14, height = 10, dpi = 300
# )

table_combined <- wrap_plots(survtable_list, ncol = 4) + 
  plot_annotation(
    title = paste0("TLS GeneSets", 
                   ifelse(surv_type == "OS", "OS", "PFS"),
                   " (", ifelse(group_method == "median", "median", "bestcut"), ")"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

print(table_combined)

#--------------------------------
#------Fig3.g
#------------------------------
validate_markers_core <- fread(file=  "./Genelist/Genes_core_4groups_loop10_pct0.25_fc0.5_0.6perc_top500_GSVA_boxplot.csv")

clusters <- c('Non-TLS', 'TLS_G3', 'TLS_G2',  'TLS_G1')

DefaultAssay(tls_balanced) <- "spatial"

# 为每个cluster类型计算模块分数
for (cluster_type in clusters) {
  # 获取当前cluster的基因列表
  cluster_genelist <- validate_markers_core$gene[validate_markers_core$cluster == cluster_type] # 使用Deseq筛选之后的基因列表
  
  # 确保基因列表不为空
  if (length(cluster_genelist) > 0) {
    # 创建临时名称避免命名冲突
    temp_name <- paste0("temp_", cluster_type)
    
    # 计算模块分数
    tls_balanced <- AddModuleScore(
      object = tls_balanced,
      features = list(cluster_genelist),
      name = temp_name
    )
    
    # 构建最终列名
    new_col_name <- paste0(cluster_type, "_score")
    
    # 获取最后添加的列（AddModuleScore在最后添加新列）
    n <- ncol(tls_balanced@meta.data)
    temp_col <- tls_balanced@meta.data[, n, drop = FALSE]
    
    # 重命名后替换
    colnames(temp_col) <- new_col_name
    tls_balanced@meta.data <- cbind(tls_balanced@meta.data[, 1:(n-1)], temp_col)
  } else {
    warning(paste("No genes found for cluster:", cluster_type))
  }
}

# 加载必要包
library(ggplot2)
library(patchwork)
library(dplyr)

# 获取起始列索引
n <- grep('TLS_Cluster_result', colnames(tls_balanced@meta.data))
score_cols <- colnames(tls_balanced@meta.data)[(n+2):(n+5)]  # 提取目标列名

# 预处理列名：替换特殊符号为下划线（解决含连字符列名报错
clean_cols <- gsub("[^[:alnum:]_]", "_", score_cols)
colnames(tls_balanced@meta.data)[(n+2):(n+5)] <- clean_cols

# 创建绘图列表
plot_list <- list()

# 循环绘制箱线图
for (col in clean_cols) {
  p <- ggplot(tls_balanced@meta.data, 
              aes(x = .data[["TLS_Cluster_result"]], 
                  y = .data[[col]], 
                  fill = .data[["TLS_Cluster_result"]])) +
    geom_violin()+
    geom_boxplot(width = 0.3, outlier.size = 0.8,fill='white') +
    labs(title = paste((col)),  # 还原可读名称
         x = "Region", 
         y = "Geneset Score") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size=10, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1,size=10, color = "black"),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"  # 统一隐藏图例
    )
  
  plot_list[[col]] <- p
}
# 拼合图形（2行4列布局）
combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 1)

options(repr.plot.width = 14, repr.plot.height = 4)
# 输出图形
print(combined_plot)

#--------------------------------
#------Fig3.k
#------------------------------

library(GSVA)
library(limma)
library(survival)
library(survminer)
library(data.table)
library(ggplot2)
library(pROC)

# 加载数据（替换为您的实际路径）
meta_melanoma <- readRDS(file = './Genelist/melanoma/Melanoma-PRJEB23709.Response_meta.Rds')
expr_melanoma <- readRDS(file = './Genelist/melanoma/Melanoma-PRJEB23709.Response.Rds')

# **关键适配1：将表达矩阵转为样本×基因格式**
expr_matrix <- as.matrix(expr_melanoma[, -1])  # 移除首列(GENE_SYMBOL)
rownames(expr_matrix) <- expr_melanoma$GENE_SYMBOL  # 基因名作为行名
expr_matrix <- t(expr_matrix)  # 转置为样本×基因（GSVA输入要求）

# **关键适配2：创建响应分组**
# 使用response_NR列：R=响应者，N=非响应者
meta_melanoma$ResponseGroup <- ifelse(
  meta_melanoma$response_NR == "R", "Responder", "Non-Responder"
)
meta_melanoma$ResponseGroup <- factor(meta_melanoma$ResponseGroup, 
                                      levels = c("Responder", "Non-Responder"))

group <- meta_melanoma$ResponseGroup[match(rownames(expr_matrix), meta_melanoma$sample_id)]

markers <- fread(file = './Genelist/Genes_core_4groups_loop10_pct0.25_fc0.5_0.6perc_top500_GSVA_boxplot.csv')  # TLS-G3基因集
clusters <- unique(markers$cluster)  # 如c("TLS_G1", "TLS_G2")

# 初始化结果存储
plot_list <- list()
auc_results <- data.frame()

for (cluster in clusters) {
  gene_list <- markers$gene[markers$cluster == cluster]
  valid_genes <- gene_list
  # valid_genes <- intersect(gene_list, diff_genes)
  
  # **关键适配5：GSVA参数设置（输入expr_matrix）**
  gsva_param <- gsvaParam(
    exprData = t(expr_matrix),  # 需转置为基因×样本
    geneSets = list(cluster = valid_genes),
    kcdf = "Gaussian"           # 标准化表达值用高斯分布
  )
  gsva_scores <- gsva(gsva_param, verbose = FALSE)
  scores <- gsva_scores[1, ]
  
  # 绘图数据
  plot_data <- data.frame(
    Sample = names(scores),
    Score = scores,
    Group = group,
    Cluster = cluster
  )
  
  # 箱线图（含统计检验)
  p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = c("#00BA38", "#F8766D")) +  # 两组颜色
    labs(title = paste0(cluster, "\n(", length(valid_genes), " genes)"),
         y = "GSVA Score",
         x = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size=10),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1,size=10)) +
    stat_compare_means(
      method = "wilcox.test",  # 两组比较使用Wilcoxon检验
      label = "p.format",
      tip.length = 0.02,
      size = 5
    )
  
  plot_list[[cluster]] <- p
  
  # 计算AUC
  roc_obj <- roc(response = ifelse(group == "Responder", 1, 0), predictor = scores)
  auc_results <- rbind(auc_results, data.frame(
    Cluster = cluster,
    AUC = auc(roc_obj),
    Genes = length(valid_genes)
  ))
}

#------ ROC
library(GSVA)
library(limma)
library(survival)
library(survminer)
library(data.table)
library(ggplot2)

# 加载数据（替换为您的实际路径）
meta_melanoma <- readRDS(file = './Genelist/melanoma/Melanoma-PRJEB23709.Response_meta.Rds')
expr_melanoma <- readRDS(file = './Genelist/melanoma/Melanoma-PRJEB23709.Response.Rds')

# **关键适配1：将表达矩阵转为样本×基因格式**
expr_matrix <- as.matrix(expr_melanoma[, -1])  # 移除首列(GENE_SYMBOL)
rownames(expr_matrix) <- expr_melanoma$GENE_SYMBOL  # 基因名作为行名
expr_matrix <- t(expr_matrix)  # 转置为样本×基因（GSVA输入要求）

# **关键适配2：创建响应分组**
# 使用response_NR列：R=响应者，N=非响应者
meta_melanoma$ResponseGroup <- ifelse(
  meta_melanoma$response_NR == "R", "Responder", "Non-Responder"
)
meta_melanoma$ResponseGroup <- factor(meta_melanoma$ResponseGroup, 
                                      levels = c("Responder", "Non-Responder"))

group <- meta_melanoma$ResponseGroup[match(rownames(expr_matrix), meta_melanoma$sample_id)]

markers <- fread(file = './Genelist/Genes_core_4groups_loop10_pct0.25_fc0.5_0.6perc_top500_GSVA_boxplot.csv')  # 替换为实际基因集
clusters <- unique(markers$cluster)  # 如c("TLS_G1", "TLS_G2")

# 初始化结果存储
plot_list <- list()

ROC_plot_list <- list()
auc_results <- data.frame()

for (cluster in clusters) {
  gene_list <- markers$gene[markers$cluster == cluster]
  valid_genes <- gene_list
  # valid_genes <- intersect(gene_list, diff_genes)
  
  # **关键适配5：GSVA参数设置（输入expr_matrix）**
  gsva_param <- gsvaParam(
    exprData = t(expr_matrix),  # 需转置为基因×样本
    geneSets = list(cluster = valid_genes),
    kcdf = "Gaussian"           # 标准化表达值用高斯分布
  )
  gsva_scores <- gsva(gsva_param, verbose = FALSE)
  scores <- gsva_scores[1, ]
  
  # 绘图数据
  plot_data <- data.frame(
    Sample = names(scores),
    Score = scores,
    Group = group,
    Cluster = cluster
  )
  
  # 箱线图（含统计检验)
  p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = c("#00BA38", "#F8766D")) +  # 两组颜色
    labs(title = paste0(cluster, "\n(", length(valid_genes), " genes)"),
         y = "GSVA Score",
         x = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size=10),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1,size=10)) +
    stat_compare_means(
      method = "wilcox.test",  # 两组比较使用Wilcoxon检验
      label = "p.format",
      tip.length = 0.02,
      size = 5
    )
  
  plot_list[[cluster]] <- p
  
  # 计算AUC
  roc_obj <- roc(response = ifelse(group == "Responder", 1, 0), predictor = scores)
  auc_results <- rbind(auc_results, data.frame(
    Cluster = cluster,
    AUC = auc(roc_obj),
    Genes = length(valid_genes))
  )
  
  # ===================== ROC曲线绘制 =====================
  # 准备二分类标签 (Responder=1, Non-Responder=0)
  binary_labels <- ifelse(group == "Responder", 1, 0)
  
  # 计算ROC对象（静默模式）
  roc_obj <- pROC::roc(
    response = binary_labels,
    predictor = scores,
    quiet = TRUE,
    levels = c(0, 1),          # 明确指定标签顺序
    direction = "<"             # 分数越高=阳性概率越高
  )
  auc_value <- pROC::auc(roc_obj)
  
  # 创建绘图数据框
  roc_data <- data.frame(
    Sensitivity = roc_obj$sensitivities,
    Specificity = roc_obj$specificities
  )
  
  # 绘制ROC曲线（关键修改点）
  roc_p <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_line(color = "steelblue", size = 1) +  # 主曲线
    geom_segment(                             # 添加对角线
      aes(x = 0, y = 0, xend = 1, yend = 1),
      linetype = "dashed", color = "gray50"
    ) +
    labs(
      title = paste0(cluster, " (n = ", length(valid_genes), ")"),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      aspect.ratio = 1,  # 保持正方形比例
      panel.grid.minor = element_blank()
    ) +
    annotate(  # AUC标注（固定在左下角）
      "text", x = 0.7, y = 0.15, 
      label = paste("AUC =", round(auc_value, 3)),
      size = 4.5, color = "red", fontface = "bold"
    )
  
  ROC_plot_list[[cluster]] <- roc_p
}

#------生存
# 加载必要包 ----------------------------------------------------------------
library(survival)
library(survminer)
library(patchwork)
library(data.table)
library(GSVA)

outpath <- "。/Genelist/melanoma/"

# 加载黑色素瘤数据 ---------------------------------------------------------
meta_melanoma <- readRDS(file = './Genelist/melanoma/Melanoma-PRJEB23709.Response_meta.Rds')
expr_melanoma <- readRDS(file = './Genelist/melanoma/Melanoma-PRJEB23709.Response.Rds')

# 处理表达矩阵 -------------------------------------------------------------
expr_matrix <- as.matrix(expr_melanoma[, -1])
rownames(expr_matrix) <- expr_melanoma$GENE_SYMBOL
expr_matrix <- t(expr_matrix)  # 转置为样本×基因

# 加载基因集 ---------------------------------------------------------------
validate_markers <- fread("./Genelist/Genes_core_4groups_loop10_pct0.25_fc0.5_0.6perc_top500_GSVA_boxplot.csv")
clusters <- unique(validate_markers$cluster)

# 设置生存分析参数 ---------------------------------------------------------
surv_type <- "OS"                  # 使用总生存期(OS)
group_method <- "median"          # 使用最优截断值分组

# 准备生存数据 -------------------------------------------------------------
# 转换生存状态：Dead=1(事件), Alive=0(删失)
meta_melanoma$status_num <- ifelse(meta_melanoma$`vital status` == "Dead", 1, 0)
surv_data_base <- data.frame(
  ID = meta_melanoma$sample_id,
  Time = meta_melanoma$`overall survival (days)`,
  Status = meta_melanoma$status_num
)

# 循环所有基因集生成生存曲线 -----------------------------------------------
survplot_list <- list()
survtable_list <- list()

for (cluster in clusters) {
  # 提取当前基因集
  gene_list <- unique(validate_markers$gene[validate_markers$cluster == cluster])
  if (length(gene_list) == 0) next
  
  # GSVA分析
  gene_sets <- list(cluster = gene_list)
  params <- gsvaParam(
    exprData = t(expr_matrix),
    geneSets = gene_sets,
    kcdf = "Gaussian",  # 标准化数据用高斯分布
    minSize = 1
  )
  gsva_scores <- gsva(params, verbose = FALSE)
  scores <- gsva_scores[1, ]
  
  # 合并GSVA分数与生存数据
  surv_data <- surv_data_base
  surv_data$GSVA_Score <- scores[match(surv_data$ID, names(scores))]
  
  # 移除缺失值
  surv_data <- na.omit(surv_data)
  
  # 分组逻辑 ------------------------------------------------------------
  if (group_method == "median") {
    median_score <- median(surv_data$GSVA_Score)
    surv_data$GSVA_Group <- ifelse(
      surv_data$GSVA_Score >= median_score,
      "High GSVA",
      "Low GSVA"
    )
    group_title <- paste0(cluster, "\n(Median Split)")
    
  } else if (group_method == "optimal") {
    # 使用survminer寻找最优截断值
    cutpoint <- surv_cutpoint(
      data = surv_data,
      time = "Time",
      event = "Status",
      variables = "GSVA_Score",
      minprop = 0.1  # 确保每组至少有10%样本
    )
    best_cutoff <- as.numeric(summary(cutpoint)$cutpoint)
    surv_data$GSVA_Group <- ifelse(
      surv_data$GSVA_Score >= best_cutoff,
      "High GSVA",
      "Low GSVA"
    )
    group_title <- paste0(cluster, "\n(Optimal Cutoff)")
  }
  
  # 拟合生存曲线
  fit <- survfit(
    Surv(Time, Status) ~ GSVA_Group,
    data = surv_data
  )
  
  # 绘制生存曲线
  surv_plot <- ggsurvplot(
    fit,
    data = surv_data,
    pval = TRUE,             # 显示log-rank检验p值
    conf.int = FALSE,        # 不显示置信区间
    risk.table = TRUE,       # 显示风险表
    palette = c("#E74C3C", "#3498DB"),  # 红/蓝配色
    legend.labs = c("High GSVA", "Low GSVA"),
    legend.title = "GSVA Group",
    xlab = "Overall Survival (Days)",
    ylab = "Survival Probability",
    ggtheme = theme_bw(base_size = 11),
    title = group_title,
    font.title = c(12, "bold"),
    censor = TRUE,           # 显示删失点
    censor.size = 3,         # 删失点大小
    censor.shape = "|",      # 删失点形状
    risk.table.height = 0.25 # 风险表高度占比
  )
  
  # 存储图形和表格
  survplot_list[[cluster]] <- surv_plot$plot
  survtable_list[[cluster]] <- surv_plot$table
}

# 图形拼接与输出 ---------------------------------------------------------
# 拼接生存曲线 (4列布局)
if (length(survplot_list) > 0) {
  surv_combined <- wrap_plots(survplot_list, ncol = 4) + 
    plot_annotation(
      title = paste0("Melanoma TLS GeneSets Survival Analysis (", surv_type, ")"),
      subtitle = paste("Grouping method:", group_method),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )
  
  # 添加共享图例
  legend <- get_legend(
    survplot_list[[1]] + 
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1))
  )
  
  final_surv_plot <- (surv_combined / legend) + 
    plot_layout(heights = c(10, 1))
  
  # 显示图形
  options(repr.plot.width = 16, repr.plot.height = 4)
  print(final_surv_plot)
  
  # 保存图形
  ggsave(
    paste0(outpath,"Melanoma_TLS_Survival_", surv_type, "_", group_method, ".pdf"),
    final_surv_plot, 
    width = 16, 
    height = 4
  )
  
  # 保存风险表
  table_combined <- wrap_plots(survtable_list, ncol = 4)
  ggsave(
    paste0(outpath,"Melanoma_TLS_RiskTable_", surv_type, "_", group_method, ".pdf"),
    table_combined,
    width = 16,
    height = 3
  )
} else {
  warning("没有可用的生存曲线生成，请检查基因集和数据匹配情况")
}

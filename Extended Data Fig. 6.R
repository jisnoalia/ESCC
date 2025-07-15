#---------------------------------------
#--------------Extended data.6a
#----------------------------------------

library(tidyverse)
library(RColorBrewer)
library(tidyr)
library(rstatix)
library(pheatmap)
library(Seurat)
library(gridExtra)

sp <- readRDS("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/ALL_HALLMARK.rds")
meta <- read.table("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/HALLMARK_custompath_TLS_meta.txt", sep = "\t", header=T)
tls_order <- readLines("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/TLS/tls_order.txt")

meta_df <- meta[meta$tissue=="AF" | meta$tissue=="ANF" | meta$tissue=="NF" | meta$tissue=="TF", ]
TAN_subset <- subset(sp, subset=tissue %in% c("AF", "ANF", "NF", "TF"))

tls <- data.frame(key = tls_order, value = c(rep("TLS_G1", times=7), rep("TLS_G2", times=42), rep("TLS_G3", times=31)), value1 = c(rep("Peri_G1", times=7), rep("Peri_G2", times=42), rep("Peri_G3", times=31)))

load("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/TLS/rctd_full_matrix.Rdata")
df <- as.data.frame(rctd_full_matrix)
df$celltype <- rownames(df)
max_rows <- apply(df, 2, function(col) rownames(df)[which.max(col)])
max_values <- apply(df, 2, max)
epi_cancer_spots <- names(max_values)[max_values > 0.5 & (max_rows == "cancer" | max_rows == "Epithelium")]

cytotoxic = c("GNLY", "IFNG", "NKG7", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK")
icb = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "BTLA", "CD274", "CD276")
TLS_Maturation_antitumor = c("STAT1", "ISG15", "IFI27", "IFI30", "MX1", "IFITM3", "GBP1", "BST2", "TAP1", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CX3CL1", "CCL28", "C7", "CFB1", "HLA-C")
Antiviral_mucosal_immune_activity = c("MUC5B", "BPIFB1", "BPIFB2", "TFF3", "DEFB1", "SCGB3A1", "PIGR", "JCHAIN", "WFDC2", "SAA1", "AZGP1", "LCN2", "PRR4", "ZG16B", "AGR2", "LTF", "SLPI", "SAA2", "CRISP3", "CXCL17")
Chemokine = c("CCL21","CCL19","CXCL12","CXCL14","CX3CL1","CCL17","CCL3","CXCL8","CCL2","XCL1","XCL2")
Chemokine_receptors = c("CCR2","CCR6","CCR5","CXCR4","CCR7","XCR1","CCR6","CCR10","CXCR3","CXCR4","CCR1","CCR3")
TLRs_and_Recptors = c("MYD88","TICAM1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR9","TLR11")
TLS_B_maturation_slgA_transport = c("IGHA1", "IGHD", "IGHM1", "IGLC1", "IGLC3", "IGLC71", "JCHAIN", "PIGR")

marker_list = list(cytotoxic, icb, TLS_Maturation_antitumor, Antiviral_mucosal_immune_activity, Chemokine, Chemokine_receptors, TLRs_and_Recptors, TLS_B_maturation_slgA_transport)
names(marker_list) = c("Cytotoxic gene", "ICB gene", "TLS Maturation antitumor", "Antiviral mucosal immune activity", "Chemokine", "Chemokine Receptors", "TLRs and Receptors", "TLS B maturation slgA transport")

TAN_subset@meta.data$tls_group_noepicaner <- "NA"
for (tls_label in tls_order) {
  core_cells <- meta_df[meta_df$TLS_cluster_label == tls_label, ]
  peripheral_cells <- meta_df[meta_df$TLS_cluster_label == "Non-TLS" & meta_df$peripheral_label_formatted == tls_label,]
  TAN_subset@meta.data$tls_group_noepicaner[TAN_subset@meta.data$cellID %in% rownames(core_cells) & !TAN_subset@meta.data$cellID %in% epi_cancer_spots] <- tls[tls$key==tls_label,]$value
  TAN_subset@meta.data$tls_group_noepicaner[TAN_subset@meta.data$cellID %in% rownames(peripheral_cells) & !TAN_subset@meta.data$cellID %in% epi_cancer_spots] <- tls[tls$key==tls_label,]$value1
}

table(TAN_subset$tls_group_noepicaner)
meta_df1 <- TAN_subset@meta.data

available_genes <- rownames(GetAssayData(TAN_subset, slot = "data", assay="SCT"))
plots <- vector("list", 8)
plots1 <- vector("list", 8)

for (i in c(1:8)) {
  genelist <- marker_list[[i]]
  matched_genes <- intersect(genelist, available_genes)
  missing_genes <- setdiff(genelist, available_genes)
  gene_lab <- names(marker_list)[i]
  print(paste0(gene_lab, ":", missing_genes))
  expr_cy <- GetAssayData(TAN_subset, slot = "data", assay="SCT")[matched_genes, ]
  exp <- data.frame()
  exp1 <- data.frame()
  for (group in c("TLS_G1","TLS_G2", "TLS_G3") ) {
    cells <- meta_df1[meta_df1$tls_group_noepicaner == group, ]
    for ( g in matched_genes ) {
      exp[group, g] <- mean(expr_cy[g,rownames(cells)], na.rm = TRUE)
    }
  }
  texp <- t(exp)
  print(gene_lab)
  p <- pheatmap(texp, border_color = "black",angle_col = "45",main=gene_lab, legend = TRUE, cluster_rows = F, cluster_cols = F, fontsize_row = 12, fontsize_col = 12, silent = TRUE)
  gt <- p$gtable
  plots[[i]] <- gt

  for (group in c("Peri_G1","Peri_G2", "Peri_G3") ) {
    cells <- meta_df1[meta_df1$tls_group_noepicaner == group, ]
    for ( g in matched_genes ) {
      exp1[group, g] <- mean(expr_cy[g,rownames(cells)], na.rm = TRUE)
    }
  }
  texp1 <- t(exp1)
  p1 <- pheatmap(texp1, border_color = "black",angle_col = "45",main=gene_lab, legend = TRUE, cluster_rows = F, cluster_cols = F, fontsize_row = 12, fontsize_col = 12, silent = TRUE)
  gt1 <- p1$gtable
  plots1[[i]] <- gt1
}
Cairo::CairoPDF(file = "geneset_exp_TLS_Peri.pdf", width = 8, height = 8) 
layout <- grid.arrange(plots[[1]], plots[[2]], ncol = 2)
layout <- grid.arrange(plots1[[1]], plots1[[2]], ncol = 2)
layout <- grid.arrange(plots[[3]], plots[[4]], ncol = 2)
layout <- grid.arrange(plots1[[3]], plots1[[4]], ncol = 2)
layout <- grid.arrange(plots[[5]], plots[[6]], ncol = 2)
layout <- grid.arrange(plots1[[5]], plots1[[6]], ncol = 2)
layout <- grid.arrange(plots[[7]], plots[[8]], ncol = 2)
layout <- grid.arrange(plots1[[7]], plots1[[8]], ncol = 2)
dev.off()

#---------------------------------------
#----------Extended data.6b
#----------------------------------------
library(tibble)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(pheatmap)
library(viridis)
library(ggsignif)
library(rstatix)
library(dplyr)
library(tidyr)

setwd('D:/Postdoc/Projects/Single_cell/ZWM')
# ----------------加载数据。
load(file = 'E:/Data_download/Baiao_cluster/spt_TANF_subset.Rdata')
# 加载细胞组成数据
load('D:/Postdoc/Projects/Single_cell/ZWM/rctd_full_matrix.Rdata')
rctd_trans <- t(rctd_full_matrix) %>% as.data.frame()

# 将细胞组成数据添加到Seurat对象
TAN_subset <- AddMetaData(TAN_subset, metadata = rctd_trans)

metadatacell <- TAN_subset@meta.data

tls_cluster_mapping <- fread("D:/Postdoc/Projects/Single_cell/ZWM/TLS_cluster_mapping.csv") %>% as.data.frame()
colnames(tls_cluster_mapping)[1]<-'TLS_cluster_label'

metadatacell <- metadatacell[,colnames(metadatacell) %in% c('TLS_cluster_label','peripheral_label_formatted',colnames(rctd_trans))]


TLS_core_cell_data <- merge(metadatacell,tls_cluster_mapping,by='TLS_cluster_label',all.x=TRUE)
TLS_core_cell_data<-TLS_core_cell_data[TLS_core_cell_data$TLS_cluster %in% c('TLS_G1','TLS_G2','TLS_G3'),]

peri_cell_data <- merge(metadatacell,tls_cluster_mapping,by.x='peripheral_label_formatted',by.y='TLS_cluster_label',all.x=TRUE)
peri_cell_data<-peri_cell_data[peri_cell_data$TLS_cluster %in% c('TLS_G1','TLS_G2','TLS_G3'),]

library(data.table)
library(dplyr)

# 1. 数据预处理优化
metadatacell <- as.data.table(metadatacell)
selected_cols <- c('TLS_cluster_label', 'peripheral_label_formatted', colnames(rctd_trans))
metadatacell <- metadatacell[, ..selected_cols]

# 2. 映射表加载与列名
tls_cluster_mapping <- fread("D:/Postdoc/Projects/Single_cell/ZWM/TLS_cluster_mapping.csv",
                             select = c("TLS_label", "TLS_cluster")) %>%
  setnames("TLS_label", "TLS_cluster_label")

# 3. 核心区域数据处理（单次合并+筛选）
TLS_core_cell_data <- metadatacell[
  tls_cluster_mapping, 
  on = "TLS_cluster_label",
  nomatch = NULL
][TLS_cluster %in% c('TLS_G1','TLS_G2','TLS_G3')]

# 4. 周边区域数据处理（避免重复合并）
peri_cell_data <- metadatacell[
  tls_cluster_mapping, 
  on = c(peripheral_label_formatted = "TLS_cluster_label"),
  nomatch = NULL
][TLS_cluster %in% c('TLS_G1','TLS_G2','TLS_G3')]

# 1. 数据准备 - 添加细胞类型过滤
all_cell_types <- setdiff(colnames(metadatacell), c('TLS_cluster_label', 'peripheral_label_formatted'))
# 关键修改：从细胞类型列表中排除目标细胞类型
cell_types <- setdiff(all_cell_types, c("Epithelium", "Cancer"))

# 2. 计算核心区平均细胞占比（使用过滤后的细胞类型）
core_avg <- metadatacell[
  !TLS_cluster_label %in% c("", "Non-TLS", NA),
  lapply(.SD, mean, na.rm = TRUE),
  by = TLS_cluster_label, 
  .SDcols = cell_types  # 使用过滤后的细胞类型列表
]

# 3. 计算周边区平均细胞占比（使用过滤后的细胞类型）
peri_avg <- metadatacell[
  !is.na(peripheral_label_formatted) & peripheral_label_formatted != "",
  lapply(.SD, mean, na.rm = TRUE), 
  by = peripheral_label_formatted, 
  .SDcols = cell_types  # 使用过滤后的细胞类型列表
] %>% 
  rename(TLS_id = peripheral_label_formatted)

# 4. 加载TLS分组信息
tls_cluster_mapping <- fread(
  "D:/Postdoc/Projects/Single_cell/ZWM/TLS_cluster_mapping.csv",
  select = c("TLS_label", "TLS_cluster")
) %>% 
  setnames(c("TLS_label", "TLS_cluster"), c("TLS_id", "TLS_cluster"))

# 5. 合并分组信息
core_avg_grouped <- merge(
  core_avg,
  tls_cluster_mapping, 
  by.x = "TLS_cluster_label",
  by.y = "TLS_id",
  all.x = FALSE
)[TLS_cluster %in% c('TLS_G1','TLS_G2','TLS_G3')]

peri_avg_grouped <- merge(
  peri_avg,
  tls_cluster_mapping,
  by = "TLS_id",
  all.x = FALSE
)[TLS_cluster %in% c('TLS_G1','TLS_G2','TLS_G3')]

# 6. 三组间比较分析（保持不变）
compare_groups <- function(data) {
  # 这里使用过滤后的细胞类型
  cell_types_ana <- setdiff(colnames(data), c("TLS_cluster_label", "TLS_id", "TLS_cluster"))
  
  result_list <- lapply(cell_types_ana, function(ct) {
    aov_model <- aov(as.formula(paste(ct, "~ TLS_cluster")), data = data)
    p_val <- summary(aov_model)[[1]]$"Pr(>F)"[1]
    data.frame(cell_type = ct, p_value = p_val)
  })
  
  diff_results <- rbindlist(result_list) %>% 
    mutate(fdr = p.adjust(p_value, method = "BH"))
  
  return(diff_results)
}

# 核心区三组比较
core_diff <- compare_groups(core_avg_grouped) %>%
  mutate(region = "Core")

# 周边区三组比较
peri_diff <- compare_groups(peri_avg_grouped) %>%
  mutate(region = "Peri")

# 7. 结果整合与可视化（保持不变）
combined_diff <- rbind(core_diff, peri_diff) %>%
  arrange(fdr)

sig_diff <- combined_diff %>% 
  filter(fdr < 0.05) %>% 
  arrange(fdr)

# 输出TOP10显著差异
print(head(sig_diff, 10))

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

# 保存图形
ggsave("./TLS_core_celltype_abundance.pdf", core_plot,
       width = 12, height = 8)
ggsave("./TLS_peri_celltype_abundance.pdf", peri_plot,
       width = 12, height = 12)


#--------------------------------
#-----Extended_data6.e
#----------------------------------

#------1.空间韦恩图

db_spt<-fread(file = './VJC_cdr3_all_TRUST4.csv')  %>%  as.data.frame() %>%
  mutate(cell_id = gsub("\\.", "_", cell_id)) %>% filter(.,Omic == 'SPT')

db_spt$tissue <- sub("_SPT$", "", db_spt$Tissue)
db_spt <- db_spt[,c('cell_id','tissue','library','Patient','vjc_cdr3')]

# A_JNP这个样本只有A一种类型
db_spt$tissue[db_spt$library == 'A_JNP'] <- 'A'

# 加载必要的包
library(dplyr)
library(readr)
library(stringr)

# 设置工作目录
setwd("./spt_barcode_tissue_type_cor_by_ZWM")

# 1. 获取所有CSV文件列表
file_list <- list.files(pattern = "\\.csv$", full.names = TRUE)

# 2. 创建空列表存储处理后的数据
result_list <- list()

# 3. 遍历处理每个文件
for (file in file_list) {
  # 从文件名提取样本编号和组织类型
  filename <- basename(file) %>% str_remove("\\.csv$")
  parts <- unlist(strsplit(filename, "-"))
  
  sample_id <- parts[1]  # 提取"-"前的样本编号
  tissue_type <- paste(parts[-1], collapse = "-")  # 提取"-"后的组织类型
  
  # 4. 读取CSV文件并添加新列
  df <- read_csv(file, show_col_types = FALSE) %>%
    rename(Barcode = 1) %>%  # 重命名第一列为Barcode
    mutate(
      cellID = paste(sample_id, Barcode, sep = "_"),  # 创建cellID
      tissue = tissue_type,  # 添加组织类型列
      library = sample_id
    ) %>%
    select(cellID, tissue)  # 只保留需要的列
  
  # 5. 添加到结果列表
  result_list[[file]] <- df
}

# 6. 合并所有数据框
final_df <- bind_rows(result_list)
# 加载必要的包
library(dplyr)

# 步骤1：通过左连接合并数据
merged_db <- left_join(
  db_spt, 
  final_df, 
  by = c("cell_id" = "cellID")  # 指定匹配列
)

# 步骤2：条件更新tissue列
db_spt_new <- merged_db %>%
  mutate(
    tissue = ifelse(
      !is.na(tissue.y),  # 若匹配成功（final_df中存在对应值）
      tissue.y,          # 使用final_df中的tissue值
      tissue.x           # 否则保留原tissue值（去除_SPT后的值）
    )
  ) %>%
  select(-tissue.x, -tissue.y)  # 删除临时列

metadata <- fread(file = './Metadata_only_TAN_sample_TLS_mark_df_labeled.csv')

options(repr.plot.width = 12,repr.plot.height = 8)

# 提取TLS的barcode和周边区域barcode
selected_spots<-metadata$cellID[which(metadata$TLS_mark==1|!metadata$peripheral_label_formatted == '')]

for (locus in c("IGK", 'IGL', 'IGH')) {
  
  # 数据筛选
  vdata <- db_spt_new[grep(locus, db_spt_new$vjc_cdr3), ]
  vdata <- vdata[vdata$cell_id %in% selected_spots, ]
  
  # 提取三种组织的独特序列
  tissue_T <- unique(vdata$vjc_cdr3[vdata$tissue == "T"])
  tissue_N <- unique(vdata$vjc_cdr3[vdata$tissue == "N"])
  tissue_A <- unique(vdata$vjc_cdr3[vdata$tissue == "A"])
  
  # 构建命名列表
  tissue_list <- list(
    T = tissue_T,
    N = tissue_N,
    A = tissue_A
  )
  # 自动计算比例并绘图
  fit <- eulerr::euler(tissue_list)
  g<- plot(fit,
           quantities = TRUE,  # 显示交集数量
           fills = c("#E69F00", "#56B4E9", "#009E73"),
           main = paste0('Spatial ', locus," Venn Diagram")
  )
  pdf(file=paste0("D:/Postdoc/Projects/Single_cell/ZWM/SPT_", locus, "_tissue_venn.pdf"),4,4)
  print(g)
  dev.off()
}

#------2.单细胞BCR韦恩图

data_scBCR = fread(file = './Add_huaxi_scBCR_db.new.tsv')
data_scBCR$cell_id<-sub("_contig_[^_]*$", "", data_scBCR$sequence_id) %>% sub("ScBCR","ScRNA",.)

cell_type_mapping_table<-data.table::fread('./meta_data.csv')

cell_type_mapping_table_base<-cell_type_mapping_table[cell_type_mapping_table$Project %in% c('Our'),c('orig.ident','Celltype_L1','Celltype_L3_add_cnv','Tissue')] %>% as.data.frame()
cell_type_mapping_table_base$cell_id=gsub("_(?!.*_)", ".", cell_type_mapping_table_base$orig.ident,perl=TRUE)

scBCR_all<-merge(data_scBCR,cell_type_mapping_table_base,by="cell_id")

# 直接使用 $ 符号创建新列
scBCR_all$vjc_cdr3 <- paste(
  scBCR_all$v_call_10x,
  scBCR_all$j_call_10x,
  scBCR_all$c_call,
  scBCR_all$junction_10x,
  sep = "_"
)
scBCR_Bcells <- scBCR_all %>% filter(.,Celltype_L1=='B/Antibody-Secreting')
scBCR_Bcells <- scBCR_Bcells[,c('Tissue','vjc_cdr3','locus')] %>% as.data.frame()

library(VennDiagram)

for (locus in c("IGK", 'IGL', 'IGH')) {
  # 数据筛选
  vdata <- scBCR_Bcells[grep(locus, scBCR_Bcells$locus), ]
  
  # 提取三种组织的独特序列（必须去重）
  tissue_T <- unique(vdata$vjc_cdr3[vdata$Tissue == "T"])
  tissue_N <- unique(vdata$vjc_cdr3[vdata$Tissue == "N"])
  tissue_A <- unique(vdata$vjc_cdr3[vdata$Tissue == "A"])
  
  # 构建命名列表
  tissue_list <- list(
    T = tissue_T,
    N = tissue_N,
    A = tissue_A
  )
  
  # 自动计算比例并绘图
  fit <- eulerr::euler(tissue_list)
  g<- plot(fit, 
           quantities = TRUE,  # 显示交集数量
           fills = c("#E69F00", "#56B4E9", "#009E73"),
           main = paste0('ScBCR', locus," Venn Diagram")
  )
  pdf(file=paste0("./Singlecell_", locus, "_tissue_venn.pdf"),4,4)
  print(g)
  dev.off()
}
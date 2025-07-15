setwd('D:/Postdoc/Projects/Single_cell/ZWM')
#--------------------------------
#------ Sup-Fig6.a
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
# ----- sup-Fig6 b
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
# ----- sup-Fig6 c-d
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
#------Sup-Fig.e
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

# 32. GSVA分析 ---------------------------------------------------------------
maker_list <- list(
  TLS_9 = c("PTGDS", "RBP5", "EIF1AY", "CETP", "SKAP1", "LAT", "CCR6", "CD1D", "CD79B"),
  
  TLS_12 = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13"),
  
  TLS_50 = c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A",
             "CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6","ATP2A3",
             "IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3",
             "INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3",
             "HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1")
)

# clusters <- c('Non-TLS', 'peri_TLS_G3', 'TLS_G3', 'peri_TLS_G2', 'TLS_G2', 'peri_TLS_G1', 'TLS_G1')
clusters <- c('TLS_9', 'TLS_12',  'TLS_50')

# Deseq2 diff gene
dds <- readRDS(file = './Genelist/Immcohort_deseq2_dds.RDS')
res <- results(dds, alpha = 0.1)  # 放宽FDR阈值
# diff_genes <- rownames(subset(res, padj < 0.1 & abs(log2FoldChange) > 0.5))
diff_genes <- rownames(subset(res, pvalue < 0.2))

# 初始化存储结果
plot_list <- list()
auc_results <- data.frame()

for (cluster in clusters) {
  valid_genes <- unique(maker_list[[cluster]])
  
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
                    subtitle = paste(length(valid_genes), "genes"),
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

options(repr.plot.width = 5, repr.plot.height = 5)

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

ggsave(plot = final_plot,filename = './Genelist/TLS_9_12_50_immcohort_responder_boxplot.pdf',width = 14,height = 6)

#--------------------------------
#------Sup-Fig6.f 
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

maker_list <- list(
  TLS_9 = c("PTGDS", "RBP5", "EIF1AY", "CETP", "SKAP1", "LAT", "CCR6", "CD1D", "CD79B"),
  
  TLS_12 = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13"),
  
  TLS_50 = c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A",
             "CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6","ATP2A3",
             "IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3",
             "INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3",
             "HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1")
)

# clusters <- c('Non-TLS', 'peri_TLS_G3', 'TLS_G3', 'peri_TLS_G2', 'TLS_G2', 'peri_TLS_G1', 'TLS_G1')
clusters <- c('TLS_9',  'TLS_12', 'TLS_50')

# 初始化存储结果
plot_list <- list()
auc_results <- data.frame()

for (cluster in clusters) {
  valid_genes <- unique(maker_list[[cluster]])
  
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
pdf(file = "./Genelist/TLS_9_12_50_ROC_Curves.pdf", width = 16, height = 4)
print(combined_roc_plot)
dev.off()

# 5. 输出AUC结果
# fwrite(auc_results, "./Genelist/core_loop10_AUC_Results.csv")

#--------------------------------
#------Sup-Fig.g
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
#------Sup-Fig.h
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


#---------------------------
#------Sup-Fig i
#---------------------------


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

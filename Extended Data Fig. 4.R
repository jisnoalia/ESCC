#Extended Data Fig4 -----
#Extended Data Fig4a -----
fc <- read_csv("~/ESCC_codex/CN/codex_CN_subCelltype_k_16_windows_20_without_epi.csv")
neighborhood_mat <- as.matrix(fc[2:ncol(fc)])
rownames(neighborhood_mat) <- paste0("CN",rownames(fc))
palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05,3, length.out=floor(palette_length/2)))

pdf('~/ESCC_codex/3-codex_without_epi_CN_k16.pdf',width = 12,height = 6)
pheatmap::pheatmap(neighborhood_mat,scale = 'none',clustering_method='average',
                   color = my_color,
                   # display_numbers = T,
                   breaks = my_breaks,
                   cellwidth = 10,
                   cellheight = 10,
                   fontsize_row = 9,
                   fontsize_col = 9,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean")
dev.off()

#Extended Data Fig4b -----
load( "/BGFS1/home/zhangzc/ESCC_codex/cor/codex_avg.Rdata")

load("/BGFS1/home/zhangzc/ESCC_codex/cor/sp_avg.Rdata")
p2r = openxlsx::read.xlsx('/BGFS1/home/zhangzc/ESCC_codex/cor/rna_protein_df_edit.xlsx')
codex <- NormalizeData(object = codex, normalization.method = "CLR", margin = 2)
codex  = ScaleData(codex)
codex3 = subset(codex,neighborhood10_new=='other',invert=T)
Idents(codex3) = 'neighborhood10_new'
avg = AverageExpression(codex3,layer='data')
avg = avg$Akoya
codex_avg = avg 
codex_avg = as.data.frame(codex_avg)
head(p2r)
table(p2r$Protein %in% rownames(codex_avg))
p2r$Protein[p2r$Protein %in% rownames(codex_avg)]
setdiff(p2r$Protein,rownames(codex_avg))
setdiff(rownames(codex_avg),p2r$Protein)
p2r$Protein = gsub("Pan.Cytokeratin","Pan-Cytokeratin",
                   gsub("HistoneH3.pSer28","HistoneH3-pSer28",
                        gsub("Keratin8.18","Keratin8-18",
                             gsub("b-Catenin","b-Catenin1",
                                  gsub('COL.1',"COL-1",
                                       gsub("PLVAP.PV.1","PLVAP-PV-1",
                                            gsub("Keratin8.18","Keratin8-18",
                                                 gsub("E.cadherin","E-cadherin",
                                                      gsub("Bcl.2","Bcl-2",p2r$Protein)))))))))
sp_avg
sp_avg %>%
  mutate(RNA = rownames(.)) %>%
  left_join(p2r,by = "RNA") %>%
  filter(RNA !='KRT18') %>%
  column_to_rownames("Protein") %>%
  dplyr::select(-RNA) -> sp_avg2


all(rownames(sp_avg2)==rownames(codex_avg))
setdiff(rownames(sp_avg2),rownames(codex_avg))
setdiff(rownames(codex_avg),rownames(sp_avg2))

sp_avg2 = sp_avg2[match(rownames(codex_avg),rownames(sp_avg2)),]
sp_avg3 = t(scale(t(sp_avg2)))
codex_avg2 = t(scale(t(codex_avg)))
# codex_avg2 = codex_avg
num_cols_sp <- ncol(sp_avg3)
num_cols_codex <- ncol(codex_avg2)
cross_correlation_matrix <- matrix(NA, nrow = num_cols_sp, ncol = num_cols_codex)
rownames(cross_correlation_matrix) <- colnames(sp_avg3)
colnames(cross_correlation_matrix) <- colnames(codex_avg2)

for (i in 1:num_cols_sp) {
  for (j in 1:num_cols_codex) {
    # 计算 matrix_A 的第 i 列和 matrix_B 的第 j 列之间的相关性
    cross_correlation_matrix[i, j] <- cor(sp_avg3[, i], codex_avg2[, j],method = 'pearson')
  }
}
cross_correlation_matrix2 = as.data.frame(cross_correlation_matrix)
cross_correlation_matrix2 = cross_correlation_matrix2[!rownames(cross_correlation_matrix2) %in% c("CC3","CC5",'CC9',"CC10","CC11","CC12"),]

palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-0.5, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05,0.8, length.out=floor(palette_length/2)))

pdf('sp_codex_cor_pheatmap.pdf',width = 8,height = 10 )
pheatmap::pheatmap(as.matrix(cross_correlation_matrix2),scale = 'none',
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   color = my_color,
                   # annotation_col = anno_cols_codex,
                   # annotation_row = anno_cols_scRNA,
                   # annotation_colors = anno_cor_colors,
                   breaks = my_breaks,
                   fontsize_row = 6,fontsize_col = 6,
                   cellwidth  = 10,cellheight  = 10)
dev.off()

library(reshape2)
cross_correlation_matrix %>% melt() %>% group_by(Var2) %>%
  top_n(1,value) ->tmp



num_cols_scRNA <- ncol(sp_avg2)
num_cols_codex <- ncol(codex_avg)
cross_correlation_matrix <- matrix(NA, nrow = num_cols_scRNA, ncol = num_cols_codex)
rownames(cross_correlation_matrix) <- colnames(sp_avg2)
colnames(cross_correlation_matrix) <- colnames(codex_avg)


for (i in 1:num_cols_scRNA) {
  for (j in 1:num_cols_codex) {
    # 计算 matrix_A 的第 i 列和 matrix_B 的第 j 列之间的相关性
    cross_correlation_matrix[i, j] <- cor(sp_avg2[, i], codex_avg[, j],method = 'kendall')
  }
}



pdf('scRNA_codex_cor_pheatmap_data.pdf',width = 8,height = 10 )
pheatmap::pheatmap(cross_correlation_matrix,scale = 'none',cluster_cols = T,width = 4,height = 4)
dev.off()
#Extended Data Fig4c -----
library(ggplot2)
library(dplyr)

theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   # 图例键为空
    legend.text = element_text(color="black",size=8), # 定义图例文本
    legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
    legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
    legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
    legend.background=element_blank()) 
}
neighborhoods$neighborhood10 <- gsub("^CM", "CN", neighborhoods$neighborhood10)

data <- read.csv("/home/users/zhangzhichao/workspace/esca/codex_CN/cn_meta_info_subCelltype_k_16_without_epi.csv", row.names = 1)
cell_color <- setNames(cell_color$color, cell_color$ct)

ratio_subtype_sample <- neighborhoods %>%
  group_by(Sample, subCelltype) %>%
  summarise(count = n()) %>%
  mutate(Ratio = count / sum(count)) %>%
  left_join(neighborhoods %>% select(Sample, Tissue) %>% distinct(), by = "Sample")


ratio_subtype_tissue <- neighborhoods %>%
  group_by(Tissue, subCelltype) %>%
  summarise(count = n()) %>%
  mutate(Ratio = count / sum(count)) %>%
  left_join(neighborhoods %>% select(Sample, Tissue) %>% distinct(), by = "Sample")

ratio_neighbor_sample <- neighborhoods %>%
  group_by(Sample, neighborhood10) %>%
  summarise(count = n()) %>%
  mutate(Ratio = count / sum(count)) %>%
  left_join(neighborhoods %>% select(Sample, Tissue) %>% distinct(), by = "Sample")

ratio_neighbor_tissue <- neighborhoods %>%
  group_by(Tissue, neighborhood10) %>%
  summarise(count = n()) %>%
  mutate(Ratio = count / sum(count)) 


p <- ggplot(ratio_neighbor_tissue ,aes(x=Tissue,y = Ratio,fill = neighborhood10))+
  geom_bar(stat = "identity", position = "fill",  width = 0.9
  ) +
  scale_fill_manual(values = niche_cols)+
  #scale_y_continuous(labels = scales::percent_format()) +
  #scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0))+
  #coord_polar()+
  #theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size=10),
    #plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.direction = 'vertical',
    #panel.grid.major=element_blank(), 
    #panel.grid.minor=element_blank(), 
    panel.background = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(ncol = 2))+
  theme_niwot()
library(cowplot)

legend <- cowplot::get_legend(p + theme(legend.position = "right"))
plot_no_legend <- p + theme(legend.position = "none")

p <- plot_grid(plot_no_legend, legend, rel_widths = c(1, 1))  # 3:1比例

ggsave(
  filename = "Tissue_neighbor.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = p,
  device = "pdf",            # 保存为 PDF
  width = 5,               # A4 宽度 (mm)
  height = 5,              # A4 高度 (mm)
  #units = "mm",              # 单位设为毫米
  dpi = 300                  # 分辨率（默认 300 DPI）
)

ggplot(ratio_subtype_sample, aes(x = Tissue, y = Ratio, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_jitter(aes(color = "black"), width = 0.2, size = 0.1, alpha = 0.8) +
  facet_wrap(~ subCelltype, scales = "free_y") +   
  scale_fill_manual(values = tissue_cols) +
  theme_void() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 



mat <- ratio_neighbor_sample %>%  
  select(-count)%>%
  select(-Tissue)%>%
  pivot_wider(
    names_from = neighborhood10,
    values_from = Ratio,
    values_fill = 0
  ) %>% column_to_rownames("Sample")
hc <- hclust(vegan::vegdist(mat, method = 'bray'), method = 'average')
sample_order <- hc$labels[hc$order]
names(niche_cols) <- gsub("^CM", "CN", names(niche_cols))


ratio_neighbor_sample$Sample <- factor(ratio_neighbor_sample$Sample, levels = sample_order)
ratio_neighbor_sample$neighborhood10 <- factor(ratio_neighbor_sample$neighborhood10, levels = paste0("CN", 1:16))
p <- ggplot(ratio_neighbor_sample,aes(x=Sample, y = Ratio, fill = neighborhood10))+
  geom_bar(stat = "identity", position = "fill", width = 0.9
  ) +
  scale_fill_manual(values = niche_cols)+
  
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0))+
  theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8),
    legend.position = "right",
    legend.direction = 'vertical',
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background = element_blank(),
    plot.margin = margin(20, 10, 10, 10)
  ) +  theme_niwot()

tissue_bar <- ratio_neighbor_sample %>%
  distinct(Sample, Tissue) %>%
  mutate(group = "Tissue")
tissue_bar$Sample <- factor(tissue_bar$Sample, levels = sample_order)

tissue_cols <- c("Adj" = "#377EB8", "Tumor" = "#E41A1C") 
tissue_plot <- ggplot(tissue_bar, aes(x = Sample, y = group, fill = Tissue)) +
  geom_tile(width = 1, height = 0.1, position = position_dodge(width = 1)) +
  scale_fill_manual(values = tissue_cols) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(axis.title = element_blank(),
        # axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank()) +
  theme_niwot()

final_plot <- p %>% insert_bottom(tissue_plot, height = 0.03)


ggsave(
  filename = "Sample_neighborhood10_cluster_label.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = final_plot,
  device = "pdf",            # 保存为 PDF
  width = 20,               # A4 宽度 (mm)
  height = 5,              # A4 高度 (mm)
  #units = "mm",              # 单位设为毫米
  dpi = 300                  # 分辨率（默认 300 DPI）
)
write.csv(neighborhoods, "neighborhoods.csv", row.names = FALSE)

#Extended Data Fig4d -----
ratio_neighbor_sample$Tissue <- factor(ratio_neighbor_sample$Tissue, levels = c("Tumor", "Adj"))
my_comparisons <- list(c("Tumor", "Adj"))
p <- ggplot(ratio_neighbor_sample, aes(x = neighborhood10, y = Ratio, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_jitter(color = 'black',  size = 0.5, alpha = 0.8,position = position_dodge(0.8)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1))+
  scale_fill_manual(values = tissue_cols) +
  stat_compare_means(method = "wilcox.test", label = "p.format", hide.ns = FALSE, size=3, label.y.npc = 0.95) +
  theme_bw(base_size =10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(
  filename = "Sample_neighborhood10_boxplot_format.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = p,
  device = "pdf",            # 保存为 PDF
  width = 15,               # A4 宽度 (mm)
  height = 5,              # A4 高度 (mm)
  #units = "mm",              # 单位设为毫米
  dpi = 300                  # 分辨率（默认 300 DPI）
)
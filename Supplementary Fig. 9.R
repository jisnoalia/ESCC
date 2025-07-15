
######ME #######
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhcihao/spatial_pathology/sub.dist_ME.Rdata")
sub.dist.all = do.call(rbind,sub.dist.ME)

sub.dist.all$cellID = sub.dist.all$id

load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhcihao/spatial_pathology/sp_meta_info.Rdata")
sub.dist.all %>% left_join(meta,by = "cellID") -> plot_df
# table(plot_df$distinct_area)
colnames(plot_df)
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
plot_df$ME = factor(plot_df$ME,
                    levels =  c('Nor-ME','Hyp-ME','MiD-ME','MoD-ME','SD&CA-ME','ICA-ME','MCA-ME'))
hallmark_paths = c("pro_fibrotic_signature",'ECM',
                   'Anti_inflammatory','Immunosuppression',
                   'pro_metastasis',"pro_fibrotic_signature",
                   'ECM','Anti_inflammatory')
hallmark_list = lapply(hallmark_paths, function(path){
  ggplot(plot_df,aes(x=min,y=plot_df[[path]] ,color=plot_df[[path]] ) )+
    geom_point_rast(size=0.1) +
    scale_color_distiller(palette = "Spectral")+
    stat_cor(size = 2)+
    geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2",span=1, size = 0.5)+
    facet_grid(~ME)+
    labs(y=path,x="Distance to Mimer (Spot)")+theme_classic() +
    theme(
      legend.key.size = unit(2, "mm"),
      panel.spacing = unit(1, "mm") ,
      axis.line = element_line(size = 0.2),
      axis.ticks = element_line(size = 0.2),
      # 全局文本大小（包括标题、坐标轴标签等）
      text = element_text(size = 5),  # 6pt
      # 坐标轴刻度文字
      axis.text = element_text(size = 6),
      # 图例文字
      legend.text = element_text(size = 5),
      # 标题文字
      #plot.title = element_text(size = 10),
      strip.text.x = element_text(size = 6),
      strip.background = element_rect(size = 0.2)
      
    )->p1
  p1[["labels"]][["colour"]] = path
  return(p1)
})
names(hallmark_list) = hallmark_paths

library(patchwork)
wrap_plots(hallmark_list,ncol = 1)->p
ggsave(
  filename = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhcihao/spatial_pathology/distinct/plot_A4_3.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = p,
  device = "pdf",            # 保存为 PDF
  width = 210,               # A4 宽度 (mm)
  height = 297,              # A4 高度 (mm)
  units = "mm",              # 单位设为毫米
  dpi = 300                  # 分辨率（默认 300 DPI）
)


#Supplementary Fig3 -----
#Supplementary Fig3b -----
pdf("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/RCTD/ISCHIA/ISCHIA/5_spatial_plots_K15_all_new.pdf", width = 7, height = 6)
for (image_name in image_names) {
  plot <- SpatialDimPlot(sp, group.by = "cc_15", images = df[df$file==image_name,"image"],stroke = NA,image.alpha = 0,pt.size.factor = 1.5) +
    scale_fill_manual(values =  color_mapping15) +
    theme_minimal() +
    ggtitle(image_name)+
    guides(fill = guide_legend(override.aes = list(size=4)))+
    coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
    scale_x_continuous(breaks = common_xbreaks) +
    scale_y_continuous(breaks = common_xbreaks) +
    theme(aspect.ratio = 1)
  print(plot)
}
dev.off()
#Supplementary Fig3c -----
development %>%
  group_by(file,development,cc_15) %>% 
  summarise(number = n()) %>%
  ungroup() %>% 
  group_by(file,development) %>%
  mutate(sum = sum(number),
         ratio = number/sum) -> dev_cc_ratio
# dev_cc_ratio %>% left_join(total_number,by = "file") %>%
#   mutate(normalize_ratio = ratio/total) -> dev_cc_ratio

dev_cc_ratio$development = factor(dev_cc_ratio$development,
                                 levels =  c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA'))
# library(gg)

dev_cc_ratio$cc_15 = factor(dev_cc_ratio$cc_15,
                            levels = paste0("CC",1:15))
compare_list = list(
  
)

my_vector = levels(dev_cc_ratio$development)
result_list = list()
for (i in 1:(length(my_vector) - 1)) {
  # 选取当前元素和下一个元素
  pair <- my_vector[i:(i + 1)]
  result_list <- c(result_list, list(pair)) # 将这对元素作为一个子列表添加到主列表中
}
result_list


dev_cc_ratio %>%
  ggplot(.,aes(x=development,y=ratio,color= development))+
  # geom_jitter_rast()+
  geom_boxplot()+
  stat_compare_means(comparisons = result_list,show.legend = F,label = "p.signif",label.y = 0.75,step.increase = 0.1)+
  scale_color_manual(values = tissue_cols)+
  facet_wrap(~cc_15,ncol=4)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,color = "black",
                                   # vjust = 0.5,
                                   hjust = 1,
                                   size = 12),
        axis.text.y = element_text(color = "black"))+
  labs(y = 'The ratio of CC in development region',x = 'Development')

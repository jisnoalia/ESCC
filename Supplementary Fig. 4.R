
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
----------------------
#------Sup-Fig4 a
---------------------
library(Seurat)
library(ggplot2)
library(RColorBrewer)

pdf(file="TLS_featureplot.pdf", width=28, height=16)
load("/lustre/home/hycui/duxiao/spatial/spatial_pathology/sp_correct.RData")
p<-subset(sp,development=='unknown',invert=T)
df <- read.table("/lustre/home/hycui/duxiao/spatial/CellTrek/sp_table.txt", sep='\t',header = T)

TLS_12 <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
TLS_9 <- c("PTGDS", "RBP5", "EIF1AY", "CETP", "SKAP1", "LAT", "CCR6", "CD1D", "CD79B")
TLS_50 <- c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A","CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6","ATP2A3","IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3","INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3","HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1")

sp <- AddModuleScore(sp,features = list(TLS_12), name="TLS_12")
sp <- AddModuleScore(sp,features = list(TLS_9), name="TLS_9")
sp <- AddModuleScore(sp,features = list(TLS_50), name="TLS_50")

colnames(sp@meta.data)[84:86]<-c("TLS_12", "TLS_9", "TLS_50")
custom_cols <- c(rev(brewer.pal(n = 11, name = "RdYlBu")), "grey")

slices <- unique(sp@meta.data$file)
for (i in c("TLS_12", "TLS_9", "TLS_50")) {
  min <- quantile(sp@meta.data[,i], probs = c(.1, .9))[[1]]
  max <- quantile(sp@meta.data[,i], probs = c(.1, .9))[[2]]
  for (s in slices) {
    p <- SpatialFeaturePlot(sp, features = i, ncol=1, images = df[df$file==s,]$image)&scale_fill_gradientn(colors = custom_cols, limits =c(min,max),oob = scales::squish)
    print(p)
  }
}
dev.off()

#-----2.绘制圈定的TLS区域

#----2.1. 加载描绘TLS区域所需的函数。
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
# right
sel_score = c('B_score','Plasma_score','Follicular_DC_score','Tfh_score','T_score','Endothelium_score','CAF_score')
common_xlim <- c(50, 550)
common_xbreaks <- seq(50, 550, by = 100)
image_names = c("QN","QP","QT1",'QLn_PLn_N','MLN_P')

sp@meta.data %>%
  dplyr::select(c(colnames(sp@meta.data)[236:242],file)) %>%
  filter(file %in% image_names)  %>%
  pivot_longer(.,cols = -file) %>%
  # group_by(file) %>%
  summarise(min = min(value),
            q05 = quantile(value,0.05),
            q95 = quantile(value,0.95),
            max = max(value)) -> hall_range 
plot_list = list()
for (pathway in sel_score) {
    for(image_name in image_names){
    p0 <- SpatialFeaturePlot(sp, features =  pathway, images = df[df$file==image_name,"image"],stroke = NA,
                             # alpha = c(1,0.6),
                             image.alpha = 0,pt.size.factor = 1.3) +
      theme_minimal() +
      scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdYlBu"))
                           ,limits = c(as.numeric(hall_range$q05),
                                       as.numeric(hall_range$q95)),
                           oob = scales::squish
      )+
      coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
      scale_x_continuous(breaks = common_xbreaks) +
      scale_y_continuous(breaks = common_xbreaks) +
      theme(
        aspect.ratio = 1,
            legend.title.position = "top" ,
            legend.justification = "center",
            legend.position = "top")+labs(title = image_name)
    
    plot_list[[paste0(image_name ,"_",pathway)]] = p0
  }
  
}


library(patchwork)
names(plot_list)
 
wrap_plots(plot_list, ncol = 5)
ggsave(filename = "TLS_score_slide_5.pdf",width = 20,height = 30)
#------------------------------------
#--- Sup-Fig4 b
#--------------------------------
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# 读取数据
tls_cluster_mapping <- fread("./TLS_cluster_mapping.csv") %>% 
  as.data.frame() %>%
  setnames(., old = "TLS_label", new = "TLS_cluster_label")  # 确保列名正确

# 提取患者ID（如ST000118 → p0118）
tls_cluster_mapping <- tls_cluster_mapping %>%
  mutate(
    # 提取"-TLS"前的完整字符串（示例：AF01ST000118）
    prefix = str_extract(TLS_cluster_label, ".*(?=-TLS)"),
    # 截取最后6位（如ST000118），再提取连续4位数字（0001 → p0001）
    patient_code = str_sub(prefix, start = -4, end = -1),
    # 保留数字部分并添加'p'前缀（0001 → p0001）
    patient = paste0("P", str_extract(patient_code, "\\d{4}"))
  ) %>%
  select(-prefix, -patient_code)  # 删除中间列

# 加载必要包
library(dplyr)
library(ggplot2)
library(patchwork)  # 用于图形拼接

# 计算各患者不同簇的占比
cluster_percent <- tls_cluster_mapping %>%
  count(patient, TLS_cluster) %>%  # 统计患者+簇的组合频次
  group_by(patient) %>%
  mutate(
    total = sum(n),               # 患者总样本数
    percent = n / total * 100     # 簇占比(%)
  ) %>%
  ungroup()

tls_cluster_colors <- c("TLS_G1" = "#E41A1C", 
                        "TLS_G2" = "#377EB8", 
                        "TLS_G3" = "#4DAF4A")

# 计算每个患者的TLS总数
tls_count <- tls_cluster_mapping %>%
  group_by(patient) %>%
  summarise(TLS_count = n())

# 创建上方图表：患者TLS数量
top_plot <- ggplot(tls_count, aes(x = patient, y = TLS_count)) +
  geom_bar(stat = "identity", fill = "grey50", alpha = 0.8) +  # 紫色柱状图
  geom_text(aes(label = TLS_count), vjust = -0.5, size = 4) +   # 添加数值标签
  ylim(0,14)+
  labs(x = NULL, y = "TLS Count", title = "Patient TLS Distribution") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),  # 隐藏X轴标签
    axis.text.y = element_text( size = 14),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.grid.major.x = element_blank()
  )

# 创建下方图表：TLS簇百分比分布
bottom_plot <- ggplot(cluster_percent, aes(x = patient, y = percent, fill = TLS_cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = tls_cluster_colors) +
  labs(x = "Patient ID", y = "Percentage", fill = "TLS Cluster") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text( size = 12),
    legend.position = "bottom"
  )

# 上下拼接图表并设置比例 (上方图占30%，下方图占70%)
combined_plot <- top_plot / bottom_plot + 
  plot_layout(heights = c(0.3, 0.7)) +
  #plot_annotation(tag_levels = 'A') &  # 添加子图标签A/B
  theme(plot.tag = element_text(size = 14, face = "bold"))


# 显示拼接后的图表
print(combined_plot)
pdf(file = "./TLS_cluster_patient_level_count_group_perc.pdf",width = 12,height = 8)
print(combined_plot)
dev.off()
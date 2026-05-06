#1 环境初始化与路径设置====
rm(list = ls()); gc()
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(patchwork) # 新增：用于拼接小提琴图和密度图

root_dir <- getwd()
input_data_dir <- file.path(root_dir, "1_Preprocessing", "processed_matrices_stand")
lasso_dir <- file.path(root_dir, "3_Lasso", "Lasso_Results_stand")
out_dir <- file.path(root_dir, "4_Plotting", "Sensitivity_Analysis", "Rank_Evolution_Grids")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_subjects <- c("SLCH018", "BJH046")

#2 数据提取与处理====
##2.1 定义数据聚合辅助函数====
aggregate_to_region <- function(mat, channel_info, target_rois) {
  mat_df <- as.data.frame(mat)
  colnames(mat_df) <- channel_info$name
  mat_df$To_Channel <- channel_info$name
  
  long_df <- pivot_longer(mat_df, cols = -To_Channel, names_to = "From_Channel", values_to = "Interaction")
  ch_region_map <- setNames(channel_info$ROI, channel_info$name)
  
  long_df$To_Region <- ch_region_map[long_df$To_Channel]
  long_df$From_Region <- ch_region_map[long_df$From_Channel]
  
  long_df <- long_df %>% filter(To_Region %in% target_rois & From_Region %in% target_rois)
  
  region_summary <- long_df %>%
    group_by(To_Region, From_Region) %>%
    summarise(Mean_Interaction = mean(Interaction, na.rm = TRUE), .groups = 'drop')
  
  full_grid <- expand.grid(To_Region = target_rois, From_Region = target_rois, stringsAsFactors = FALSE)
  region_summary <- full_grid %>%
    left_join(region_summary, by = c("To_Region", "From_Region"))
  
  region_mat_df <- pivot_wider(region_summary, names_from = From_Region, values_from = Mean_Interaction)
  region_mat <- as.matrix(region_mat_df[, -1])
  rownames(region_mat) <- region_mat_df$To_Region
  region_mat <- region_mat[target_rois, target_rois]
  
  return(region_mat)
}

##2.2 主循环提取所有秩的结果====
all_ranks_data <- data.frame()

for (sub_id in target_subjects) {
  preproc_file <- file.path(input_data_dir, paste0(sub_id, ".rds"))
  
  if(!file.exists(preproc_file)) {
    next
  }
  
  preproc_data <- readRDS(preproc_file)
  ch_info <- preproc_data$channels
  actual_rois <- sort(unique(na.omit(ch_info$ROI)))
  
  sub_lasso_dir <- file.path(lasso_dir, sub_id)
  lasso_files <- list.files(sub_lasso_dir, pattern = "\\.rds$", full.names = TRUE)
  
  for (file in lasso_files) {
    res <- readRDS(file)
    rank_str <- str_extract(basename(file), "r_\\d+")
    r_val <- as.numeric(str_extract(rank_str, "\\d+"))
    
    pi_base_rec  <- abs(res$Baseline_Recall$Pi)
    pi_delay_rec <- abs(res$Delay_Recall$Pi)
    pi_base_man  <- abs(res$Baseline_Manip$Pi)
    pi_delay_man <- abs(res$Delay_Manip$Pi)
    
    reg_base_rec  <- aggregate_to_region(pi_base_rec, ch_info, actual_rois)
    reg_delay_rec <- aggregate_to_region(pi_delay_rec, ch_info, actual_rois)
    reg_base_man  <- aggregate_to_region(pi_base_man, ch_info, actual_rois)
    reg_delay_man <- aggregate_to_region(pi_delay_man, ch_info, actual_rois)
    
    reg_contrast <- (reg_delay_man - reg_base_man) - (reg_delay_rec - reg_base_rec)
    
    df_long <- as.data.frame(as.table(reg_contrast))
    colnames(df_long) <- c("To_Region", "From_Region", "Interaction")
    df_long$Rank <- r_val
    df_long$Subject <- sub_id
    
    df_long$Edge_Name <- paste0(df_long$From_Region, " -> ", df_long$To_Region)
    
    all_ranks_data <- bind_rows(all_ranks_data, df_long)
  }
}

#3 划分秩区与连边排序====
##3.1 划分各个受试者的拐点中心（秩区）====
all_ranks_data <- all_ranks_data %>%
  mutate(Rank_Group = case_when(
    Subject == "SLCH018" & Rank <= 25 ~ "1st",
    Subject == "SLCH018" & Rank > 25 & Rank <= 50 ~ "2nd",
    Subject == "SLCH018" & Rank > 50 ~ "3rd",
    
    Subject == "BJH046" & Rank <= 30 ~ "1st",
    Subject == "BJH046" & Rank > 30 & Rank <= 70 ~ "2nd",
    Subject == "BJH046" & Rank > 70 ~ "3rd",
    
    Subject == "COG022" & Rank <= 20 ~ "1st",
    Subject == "COG022" & Rank > 20 ~ "2nd",
    
    Subject == "COG023" & Rank <= 25 ~ "1st",
    Subject == "COG023" & Rank > 25 ~ "2nd",
    
    TRUE ~ NA_character_
  ))

##3.2 提取每个受试者的最高秩区计算SD并降序排列====
# 锁定高秩区：SLCH018 和 BJH046 为 "3rd"（第三拐点中心之后），其余为 "2nd"（第二拐点中心之后）
edge_sd_ranking <- all_ranks_data %>%
  filter(
    (Subject %in% c("SLCH018", "BJH046") & Rank_Group == "3rd") | 
      (!Subject %in% c("SLCH018", "BJH046") & Rank_Group == "2nd")
  ) %>%
  group_by(Subject, Edge_Name) %>%
  summarise(SD = sd(Interaction, na.rm = TRUE), .groups = 'drop') %>%
  filter(!is.na(SD)) %>%
  # 按标准差从大到小降序排列，以便让波动最大的排在最前
  arrange(Subject, desc(SD)) 

#4 分块绘制 8图网格(2x4)并按受试者分文件夹保存====
unique_subs <- unique(edge_sd_ranking$Subject)

for (sub in unique_subs) {
  sub_out_dir <- file.path(out_dir, sub)
  dir.create(sub_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  sub_ranked_edges <- edge_sd_ranking %>% 
    filter(Subject == sub) %>% 
    pull(Edge_Name)
  
  sub_plot_data <- all_ranks_data %>% filter(Subject == sub)
  sub_plot_data$Edge_Name <- factor(sub_plot_data$Edge_Name, levels = sub_ranked_edges)
  
  chunk_size <- 8
  num_chunks <- ceiling(length(sub_ranked_edges) / chunk_size)
  
  message(sprintf(">>> 开始绘制 %s 的折线图，共 %d 个连边，按【高秩区SD】降序排布...", sub, length(sub_ranked_edges)))
  
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(sub_ranked_edges))
    current_edges <- sub_ranked_edges[start_idx:end_idx]
    
    chunk_data <- sub_plot_data %>% filter(Edge_Name %in% current_edges)
    
    p <- ggplot(chunk_data, aes(x = Rank, y = Interaction, group = Rank_Group)) +
      geom_hline(yintercept = 0, color = "indianred", linetype = "dashed", linewidth = 0.8) +
      geom_line(color = "steelblue", linewidth = 1) +
      geom_point(color = "black", size = 1.5) +
      facet_wrap(~ Edge_Name, ncol = 2, scales = "free_y") + 
      theme_bw() +
      labs(
        title = paste0("Interaction Contrast vs. Cointegration Rank - ", sub, " (Part ", i, ")"),
        subtitle = "Sorted by standard deviation in the high-rank region (descending) | Disconnected across turning points",
        x = "Cointegration Rank",
        y = "Interaction Contrast (Manipulation - Recall)"
      ) +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, face = "italic", hjust = 0.5),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        panel.spacing = unit(1, "lines")
      )
    
    part_str <- sprintf("%02d", i)
    output_filename <- file.path(sub_out_dir, paste0(sub, "_Evolution_Grid_Part_", part_str, ".png"))
    
    ggsave(output_filename, plot = p, width = 10, height = 12, dpi = 300)
  }
}

message("\n>>> ✅ 所有受试者的分块网格图绘制完毕！")

#5 绘制各受试者不同秩区的SD分布联合图====
##5.1 计算所有秩区的SD====
sd_dist_data <- all_ranks_data %>%
  filter(!is.na(Rank_Group)) %>%
  group_by(Subject, Edge_Name, Rank_Group) %>%
  summarise(SD = sd(Interaction, na.rm = TRUE), .groups = 'drop') %>%
  filter(!is.na(SD)) %>%
  mutate(Rank_Group_Label = paste0(Rank_Group, " Turning Point"))

# 设定统一的因子顺序以保证图例和 X 轴顺序正确
label_levels <- c("1st Turning Point", "2nd Turning Point", "3rd Turning Point")
sd_dist_data$Rank_Group_Label <- factor(sd_dist_data$Rank_Group_Label, levels = label_levels)

# 设定贴近原图的颜色映射
color_map <- c(
  "1st Turning Point" = "#E07B7B", 
  "2nd Turning Point" = "#E4C05B", 
  "3rd Turning Point" = "#71A2C9"  
)

##5.2 循环绘制并保存联合图====
for (sub in unique_subs) {
  sub_sd_data <- sd_dist_data %>% filter(Subject == sub)
  
  # 上半部分：小提琴图 + 箱线图
  p_violin <- ggplot(sub_sd_data, aes(x = Rank_Group_Label, y = SD, fill = Rank_Group_Label)) +
    geom_violin(alpha = 0.6, trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
    theme_bw() +
    labs(x = NULL, y = "Standard Deviation of Interaction Contrasts") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12)
    ) +
    scale_fill_manual(values = color_map)
  
  # 下半部分：核密度估计图
  p_density <- ggplot(sub_sd_data, aes(x = SD, fill = Rank_Group_Label)) +
    geom_density(alpha = 0.6, color = "grey50") +
    theme_bw() +
    labs(
      x = "Standard Deviation of Interaction Contrasts", 
      y = "Density", 
      fill = "Turning Point"
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    scale_fill_manual(values = color_map)
  
  # 利用 patchwork 拼合上下两图
  p_combined <- p_violin / p_density +
    plot_annotation(
      title = paste0("Standard Deviation Distribution of Interaction Contrasts - ", sub),
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  
  # 保存图表
  sub_out_dir <- file.path(out_dir, sub)
  output_filename <- file.path(sub_out_dir, paste0(sub, "_SD_Distribution_Combined.png"))
  ggsave(output_filename, plot = p_combined, width = 10, height = 12, dpi = 300)
}

message("\n>>> ✅ 所有受试者的 SD 分布联合图绘制完毕！")
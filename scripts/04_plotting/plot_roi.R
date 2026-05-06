#1 环境初始化与路径设置====
rm(list = ls()); gc()
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

root_dir <- getwd()
input_data_dir <- file.path(root_dir, "1_Preprocessing", "processed_matrices_stand")
lasso_dir <- file.path(root_dir, "3_Lasso", "Lasso_Results_stand")
heatmap_dir <- file.path(root_dir, "4_Plotting", "Region_Heatmaps_stand")
bar_dir <- file.path(root_dir, "4_Plotting", "Barplots_11ROI_stand")

dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(bar_dir, showWarnings = FALSE, recursive = TRUE)

target_subjects <- list.dirs(lasso_dir, full.names = FALSE, recursive = FALSE)
target_rois <- c("A1", "STG", "STS", "MTG", "ACC", "VLPFC", "DLPFC", "M1", "HPC", "PHC", "parietal")

# 初始化存储所有受试者长格式数据的空数据框
all_subjects_data <- data.frame()

#2 定义脑区聚合辅助函数====
aggregate_to_region <- function(mat, channel_info, target_rois) {
  mat_df <- as.data.frame(mat)
  colnames(mat_df) <- channel_info$name
  mat_df$To_Channel <- channel_info$name
  
  long_df <- pivot_longer(mat_df, cols = -To_Channel, names_to = "From_Channel", values_to = "Interaction")
  ch_region_map <- setNames(channel_info$ROI, channel_info$name)
  
  long_df$To_Region <- ch_region_map[long_df$To_Channel]
  long_df$From_Region <- ch_region_map[long_df$From_Channel]
  
  long_df <- long_df %>%
    filter(To_Region %in% target_rois & From_Region %in% target_rois)
  
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

#3 定义脑区热图绘制函数====
plot_region_heatmap <- function(mat, title, filename, width = 8, height = 7) {
  regions <- rownames(mat)
  df <- as.data.frame(mat)
  df$Row <- factor(rownames(df), levels = rev(regions))
  
  df_long <- pivot_longer(df, cols = -Row, names_to = "Col", values_to = "Interaction")
  df_long$Col <- factor(df_long$Col, levels = regions)
  
  max_val <- max(abs(df_long$Interaction), na.rm = TRUE)
  if (max_val == 0) max_val <- 1 
  
  p <- ggplot(df_long, aes(x = Col, y = Row, fill = Interaction)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                         limits = c(-max_val, max_val), name = "Interaction", na.value = "grey90") +
    theme_minimal() +
    labs(title = title, x = "From Region", y = "To Region") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
}

# 辅助函数：将矩阵转为方便合并的长格式
mat_to_long <- function(mat, cond_name, sub_id) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("To_Region", "From_Region", "Interaction")
  df$Condition <- cond_name
  df$Subject <- sub_id
  return(df)
}

#4 主循环处理提取数据与画热图====
for (sub_id in target_subjects) {
  
  sub_lasso_dir <- file.path(lasso_dir, sub_id)
  
  # 获取所有 rds 文件并提取秩最大的文件
  all_lasso_files <- list.files(sub_lasso_dir, pattern = "\\.rds$", full.names = TRUE)
  ranks <- as.numeric(gsub(".*_r_(\\d+)\\.rds$", "\\1", basename(all_lasso_files)))
  lasso_file <- all_lasso_files[which.max(ranks)]
  
  preproc_file <- file.path(input_data_dir, paste0(sub_id, ".rds"))
  
  lasso_res <- readRDS(lasso_file)
  preproc_data <- readRDS(preproc_file)
  ch_info <- preproc_data$channels
  
  sub_out_dir <- file.path(heatmap_dir, sub_id)
  dir.create(sub_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  pi_base_rec  <-  abs(lasso_res$Baseline_Recall$Pi)
  pi_delay_rec <-  abs(lasso_res$Delay_Recall$Pi)
  pi_base_man  <-  abs(lasso_res$Baseline_Manip$Pi)
  pi_delay_man <-  abs(lasso_res$Delay_Manip$Pi)
  
  reg_base_rec  <- aggregate_to_region(pi_base_rec, ch_info, target_rois)
  reg_delay_rec <- aggregate_to_region(pi_delay_rec, ch_info, target_rois)
  reg_base_man  <- aggregate_to_region(pi_base_man, ch_info, target_rois)
  reg_delay_man <- aggregate_to_region(pi_delay_man, ch_info, target_rois)
  
  reg_inter_recall <- reg_delay_rec - reg_base_rec
  reg_inter_manip  <- reg_delay_man - reg_base_man
  reg_inter_contrast <- reg_inter_manip - reg_inter_recall
  
  # 记录当前受试者数据以便后续画条形图
  df_recall <- mat_to_long(reg_inter_recall, "Recall", sub_id)
  df_manip  <- mat_to_long(reg_inter_manip, "Manip", sub_id)
  df_contrast <- mat_to_long(reg_inter_contrast, "Contrast", sub_id)
  all_subjects_data <- bind_rows(all_subjects_data, df_recall, df_manip, df_contrast)
  
  # 画个体热图
  plot_region_heatmap(mat = reg_inter_recall, 
                      title = paste0("Interaction Change (Delay - Baseline) - Recall - ", sub_id), 
                      filename = file.path(sub_out_dir, paste0(sub_id, "_Region_Interaction_Recall.png")))
  
  plot_region_heatmap(mat = reg_inter_manip, 
                      title = paste0("Interaction Change (Delay - Baseline) - Manipulation - ", sub_id), 
                      filename = file.path(sub_out_dir, paste0(sub_id, "_Region_Interaction_Manip.png")))
  
  plot_region_heatmap(mat = reg_inter_contrast, 
                      title = paste0("Interaction Contrast (Manipulation - Recall) - ", sub_id), 
                      filename = file.path(sub_out_dir, paste0(sub_id, "_Region_Interaction_Contrast.png")))
}

#5 组级数据分布条形图绘制 (11 * 11 * 3 = 363张)====
conditions <- c("Recall", "Manip", "Contrast")

for (cond in conditions) {
  for (from_roi in target_rois) {
    # 按照 条件/发送区域 创建文件夹，例如 4_Plotting/Barplots_11ROI/Recall/A1
    cond_from_dir <- file.path(bar_dir, cond, from_roi)
    dir.create(cond_from_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (to_roi in target_rois) {
      # 过滤当前指定的连线数据，并剔除缺失数据的被试 (NA)
      plot_data <- all_subjects_data %>%
        filter(Condition == cond, From_Region == from_roi, To_Region == to_roi) %>%
        filter(!is.na(Interaction))
      
      # 只有当该连线有数据时才作图
      if (nrow(plot_data) > 0) {
        
        # 按照 Interaction 排序被试，使条形图呈现阶梯状分布
        plot_data <- plot_data %>%
          arrange(Interaction) %>%
          mutate(Subject = factor(Subject, levels = Subject))
        
        # 动态设置完整清晰的图表标题和 Y 轴标签
        if (cond == "Contrast") {
          plot_title <- paste0("Interaction (", from_roi, " -> ", to_roi, ") Contrast (Manipulation - Recall)")
          y_label <- "Interaction Contrast"
        } else if (cond == "Recall") {
          plot_title <- paste0("Interaction (", from_roi, " -> ", to_roi, ") Change (Delay - Baseline) - Recall")
          y_label <- "Interaction Change"
        } else if (cond == "Manip") {
          plot_title <- paste0("Interaction (", from_roi, " -> ", to_roi, ") Change (Delay - Baseline) - Manipulation")
          y_label <- "Interaction Change"
        }
        
        p <- ggplot(plot_data, aes(x = Subject, y = Interaction)) +
          geom_col(fill = "steelblue", color = "black", alpha = 0.8) +
          geom_hline(yintercept = 0, color = "indianred", linetype = "dashed", linewidth = 1) +
          theme_bw() +
          labs(title = plot_title,
               subtitle = paste0("Valid Subjects: ", nrow(plot_data)),
               x = "Subject",
               y = y_label) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9) # 倾斜受试者名字防止重叠
          )
        
        filename <- file.path(cond_from_dir, paste0(cond, "_", from_roi, "_to_", to_roi, ".png"))
        ggsave(filename, plot = p, width = 8, height = 5, dpi = 300)
      }
    }
  }
}
# 1 环境初始化与路径设置 ====
rm(list = ls()); gc()
library(ggplot2)
library(dplyr)
library(tidyr)

root_dir <- getwd()
input_data_dir <- file.path(root_dir, "1_Preprocessing", "processed_matrices_correct")
lasso_dir <- file.path(root_dir, "3_Lasso", "Lasso_Results")
output_dir <- file.path(root_dir, "4_Plotting", "Matrix_Heatmaps")

dir.create(output_dir, recursive = TRUE)

# 获取所有受试者 Lasso 结果目录
target_subjects <- list.dirs(lasso_dir, full.names = FALSE, recursive = FALSE)

# 2 定义热图绘制辅助函数 ====
plot_heatmap <- function(mat, row_labels, col_labels, title, filename, width = 12, height = 15) {
  df <- as.data.frame(mat)
  rownames(df) <- row_labels
  colnames(df) <- col_labels
  df$Row <- factor(rownames(df), levels = rev(row_labels))
  
  df_long <- pivot_longer(df, cols = -Row, names_to = "Col", values_to = "Value")
  df_long$Col <- factor(df_long$Col, levels = col_labels)
  
  max_val <- max(abs(df_long$Value), na.rm = TRUE)
  if (max_val == 0) max_val <- 1 
  
  p <- ggplot(df_long, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                         limits = c(-max_val, max_val), name = "Value") +
    theme_minimal() +
    labs(title = title, x = "", y = "") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300, limitsize = FALSE)
}


# 3 主循环 ====
for (sub_id in target_subjects) {
  
  # 3.1 直接定位和加载数据
  sub_lasso_dir <- file.path(lasso_dir, sub_id)
  lasso_file <- list.files(sub_lasso_dir, pattern = "\\.rds$", full.names = TRUE)[1] 
  preproc_file <- file.path(input_data_dir, paste0(sub_id, ".rds"))
  
  lasso_res <- readRDS(lasso_file)
  preproc_data <- readRDS(preproc_file)
  
  # 3.2 提取通道标签
  ch_info <- preproc_data$channels
  channel_labels <- paste0(ch_info$name, " (", ch_info$region, ")")
  
  sub_out_dir <- file.path(output_dir, sub_id)
  dir.create(sub_out_dir, recursive = TRUE)
  
  num_channels <- length(channel_labels)
  dynamic_height <- max(10, num_channels * 0.15)
  
  # 3.3 遍历所有实验条件绘制基础矩阵
  for (cond_name in names(lasso_res)) {
    cond_data <- lasso_res[[cond_name]]
    
    alpha_mat <- cond_data$alpha
    beta_mat <- cond_data$beta
    pi_mat <- cond_data$Pi
    
    r <- cond_data$r
    rank_labels <- paste0("Rank_", 1:r)
    
    plot_heatmap(mat = beta_mat, 
                 row_labels = channel_labels, 
                 col_labels = rank_labels, 
                 title = paste0("Beta Matrix - ", sub_id, " - ", cond_name), 
                 filename = file.path(sub_out_dir, paste0(sub_id, "_", cond_name, "_1_Beta.png")),
                 width = max(6, r * 1.5), height = dynamic_height)
    
    plot_heatmap(mat = alpha_mat, 
                 row_labels = channel_labels, 
                 col_labels = rank_labels, 
                 title = paste0("Alpha Matrix - ", sub_id, " - ", cond_name), 
                 filename = file.path(sub_out_dir, paste0(sub_id, "_", cond_name, "_2_Alpha.png")),
                 width = max(6, r * 1.5), height = dynamic_height)
    
    plot_heatmap(mat = pi_mat, 
                 row_labels = channel_labels, 
                 col_labels = channel_labels, 
                 title = paste0("PI Interaction Matrix - ", sub_id, " - ", cond_name), 
                 filename = file.path(sub_out_dir, paste0(sub_id, "_", cond_name, "_3_PI.png")),
                 width = dynamic_height, height = dynamic_height)
  }
  
  # 3.4 绘制三种差异网络矩阵 (Delay - Baseline) 以及 (Manip - Recall)
  
  # (1) Recall: Delay - Baseline
  pi_diff_recall <- lasso_res$Delay_Recall$Pi - lasso_res$Baseline_Recall$Pi
  plot_heatmap(mat = pi_diff_recall, 
               row_labels = channel_labels, 
               col_labels = channel_labels, 
               title = paste0("PI Diff (Delay - Baseline) - ", sub_id, " - Recall"), 
               filename = file.path(sub_out_dir, paste0(sub_id, "_Diff_Recall_PI.png")),
               width = dynamic_height, height = dynamic_height)
  
  # (2) Manip: Delay - Baseline
  pi_diff_manip <- lasso_res$Delay_Manip$Pi - lasso_res$Baseline_Manip$Pi
  plot_heatmap(mat = pi_diff_manip, 
               row_labels = channel_labels, 
               col_labels = channel_labels, 
               title = paste0("PI Diff (Delay - Baseline) - ", sub_id, " - Manip"), 
               filename = file.path(sub_out_dir, paste0(sub_id, "_Diff_Manip_PI.png")),
               width = dynamic_height, height = dynamic_height)
  
  # (3) Contrast: Manip 差异网络 - Recall 差异网络
  pi_diff_contrast <- pi_diff_manip - pi_diff_recall
  plot_heatmap(mat = pi_diff_contrast, 
               row_labels = channel_labels, 
               col_labels = channel_labels, 
               title = paste0("PI Contrast (Manip - Recall) - ", sub_id), 
               filename = file.path(sub_out_dir, paste0(sub_id, "_Contrast_Manip_vs_Recall_PI.png")),
               width = dynamic_height, height = dynamic_height)
}
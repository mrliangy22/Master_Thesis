library(ca4eeg)
library(urca)

# 1. 路径与参数设置 ====
processed_dir <- "1_Preprocessing/processed_matrices_downsample"
debug_dir <- "0_debug"
dir.create(debug_dir, showWarnings = FALSE, recursive = TRUE)

# 直接指定需要处理的受试者和频率序列
subjects <- c("CP40", "LL31")
fs_seq   <- c(20, 25, 40, 50, 80, 100, 200, 400)
n_fs     <- length(fs_seq)

# 2. 核心分析循环 ====
for (sub_id in subjects) {
  message(sprintf("\n>>> 开始处理受试者: %s", sub_id))
  
  # 读取一个最高频率的文件以获取通道数量
  ref_data <- readRDS(file.path(processed_dir, sprintf("%s_400Hz.rds", sub_id)))
  n_channels <- ncol(ref_data$Baseline_All_Combined$Z0)
  rm(ref_data)
  
  # 初始化每个通道的存储矩阵 [n_fs x 2] (Gamma, T-Stat)
  res_adf_by_channel <- vector("list", n_channels)
  for (ch in 1:n_channels) {
    res_adf_by_channel[[ch]] <- matrix(NA, nrow = n_fs, ncol = 2)
    colnames(res_adf_by_channel[[ch]]) <- c("Gamma", "TestStat")
  }
  
  # 遍历所有采样率文件
  for (i in seq_along(fs_seq)) {
    current_fs <- fs_seq[i]
    file_path <- file.path(processed_dir, sprintf("%s_%dHz.rds", sub_id, current_fs))
    
    data_all <- readRDS(file_path)
    block <- data_all$Baseline_All_Combined
    
    # 提取预处理阶段已经计算好的差分矩阵(Z0)和滞后矩阵(Z1)
    dy_mat <- block$Z0
    ly_mat <- block$Z1
    
    # 逐通道进行线性回归
    for (ch in 1:n_channels) {
      dy <- dy_mat[, ch]
      ly <- ly_mat[, ch]
      
      mod <- lm(dy ~ ly)
      
      res_adf_by_channel[[ch]][i, 1] <- coef(mod)[["ly"]]
      res_adf_by_channel[[ch]][i, 2] <- summary(mod)$coefficients["ly", "t value"]
    }
    
    cat(sprintf("Processed FS: %3d Hz | Ch1 Gamma: %.4f | Ch1 T-stat: %.2f\n",
                current_fs, res_adf_by_channel[[1]][i, 1], res_adf_by_channel[[1]][i, 2]))
  }
  
  # 3. 为当前受试者绘图并保存 ====
  plot_path <- file.path(debug_dir, sprintf("ADF_Sensitivity_Analysis_%s.png", sub_id))
  message(sprintf("\n>>> Plotting ADF Stats to '%s'...", plot_path))
  
  # 动态计算布局与图片高度，防止覆盖
  n_cols <- 4
  n_rows <- ceiling(n_channels / n_cols)
  img_height <- 450 * n_rows 
  
  png(plot_path, width = 1600, height = img_height, res = 120)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 4))
  
  # ADF 5% 临界值
  crit_val_5pct <- -2.86
  
  for (ch in 1:n_channels) {
    dat <- res_adf_by_channel[[ch]]
    gammas <- dat[, "Gamma"]
    stats  <- dat[, "TestStat"]
    
    plot(fs_seq, gammas, type = "o", pch = 16, col = "red",
         xlab = "Sampling Rate (Hz)", ylab = "Gamma Coeff",
         main = paste(sub_id, "- Channel", ch),
         axes = FALSE)
    axis(1)
    axis(2, col = "red", col.axis = "red")
    box()
    
    par(new = TRUE)
    ylim_stats <- range(c(stats, crit_val_5pct, -3.5), na.rm = TRUE)
    plot(fs_seq, stats, type = "o", pch = 17, col = "blue",
         xlab = "", ylab = "", axes = FALSE, ylim = ylim_stats)
    axis(4, col = "blue", col.axis = "blue")
    mtext("t-Statistic", side = 4, line = 2.5, col = "blue", cex = 0.7)
    
    abline(h = crit_val_5pct, col = "blue", lty = 2, lwd = 1.5)
    text(min(fs_seq), crit_val_5pct, "5% Crit", pos = 3, col = "blue", cex = 0.8)
    grid()
  }
  
  dev.off()
}

par(mfrow = c(1, 1))
message(">>> 全部 ADF 分析完成。")
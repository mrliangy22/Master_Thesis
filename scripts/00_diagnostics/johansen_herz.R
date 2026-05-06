#1 环境初始化与路径设置====
library(ca4eeg)
library(tsDyn)
library(ggplot2)

processed_dir <- "1_Preprocessing/processed_matrices_downsample"
debug_dir <- "0_debug"
if (!dir.exists(debug_dir)) dir.create(debug_dir, recursive = TRUE)

# 定义受试者与采样率序列
subjects <- c("CP40", "LL31")
fs_seq <- c(20, 25, 40, 50, 80, 100, 200, 400)
n_fs <- length(fs_seq)

#2 受试者批量处理主循环====
for (sub_id in subjects) {
  message(sprintf("\n>>> [Subject] 正在处理: %s", sub_id))
  
  ##2.1 探测系统维度与初始化====
  # 读取最高采样率文件获取通道数 p
  init_file <- file.path(processed_dir, paste0(sub_id, "_400Hz.rds"))
  if (!file.exists(init_file)) {
    message(sprintf("    -> [跳过] 找不到初始化文件: %s", init_file))
    next
  }
  data_tmp <- readRDS(init_file)
  n_vars <- nrow(data_tmp$channels)
  rm(data_tmp); gc()
  
  # 存储矩阵
  res_lambdas    <- matrix(NA, nrow = n_fs, ncol = n_vars) 
  res_test_stats <- matrix(NA, nrow = n_fs, ncol = n_vars) 
  
  ##2.2 提取 tsDyn 临界值 (95% 置信度)====
  # 使用 Doornik (1998) 的 Gamma 近似法获取针对 p 维系统的临界值
  # nmp 代表 p - r，这里计算从 r=0 到 r=p-1 的所有临界值
  trace_cv_mat <- tsDyn:::gamma_doornik_all(
    q = 0.95,
    nmp = n_vars:1,
    test = "H_lc", # 对应含截距项的 VECM 模型
    type = "trace"
  )
  crit_vals_95 <- as.numeric(trace_cv_mat[, "95%"])
  
  ##2.3 遍历采样率计算 Johansen====
  for(i in seq_along(fs_seq)) {
    current_fs <- fs_seq[i]
    rds_path <- file.path(processed_dir, paste0(sub_id, "_", current_fs, "Hz.rds"))
    
    if (!file.exists(rds_path)) next
    
    data_all <- readRDS(rds_path)
    block <- data_all$Baseline_All_Combined
    
    # 调用 C++ 核心进行 Johansen 迹检验计算
    fit_fs <- johansen(Y = block$Z0, Z = block$Z1, r = 1, dt = data_all$dt, 
                       intercept = TRUE, debug = FALSE)
    
    res_lambdas[i, ]    <- fit_fs$lambda[1:n_vars]
    res_test_stats[i, ] <- fit_fs$test[1:n_vars]
    
    cat(sprintf("    -> Fs: %3d Hz | p: %d | Max Stat: %.2f\n", 
                current_fs, n_vars, max(res_test_stats[i, ], na.rm=TRUE)))
    rm(data_all, block, fit_fs); gc()
  }
  
  #3 绘图与可视化保存====
  cols <- rainbow(n_vars, v = 0.85)
  
  ###3.1 特征值演变图====
  png_eig <- file.path(debug_dir, paste0(sub_id, "_Johansen_Eigenvalues.png"))
  png(png_eig, width = 800, height = 600, res = 120)
  par(mfrow = c(1, 1), mar = c(4, 5, 3, 2) + 0.1)
  
  plot(range(fs_seq), range(res_lambdas, na.rm = TRUE), type = "n",
       xlab = "Fs (Hz)", ylab = "Eigenvalue (Lambda)",
       main = paste("Eigenvalues vs Sampling Rate -", sub_id), las = 1)
  grid(col = "gray", lty = "dotted")
  for(k in 1:n_vars) {
    lines(fs_seq, res_lambdas[, k], type = "o", pch = 16, col = cols[k], lwd = 2, cex = 0.7)
  }
  legend("topright", legend = paste0("L", 1:n_vars), col = cols, lty=1, lwd=2, bty="n", cex=0.6, ncol=2)
  dev.off()
  
  ###3.2 迹统计量子图 (强化临界值标注)====
  png_trace <- file.path(debug_dir, paste0(sub_id, "_Johansen_TraceStats.png"))
  n_cols <- 3
  n_rows <- ceiling(n_vars / n_cols)
  png(png_trace, width = 1200, height = 320 * n_rows, res = 120)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  
  for(k in 1:n_vars) {
    r_val <- k - 1 
    current_cv <- crit_vals_95[k]
    
    # 动态计算纵轴，确保临界值和统计量都在视野内
    y_vals <- c(res_test_stats[, k], current_cv)
    y_limits <- c(min(y_vals, na.rm = TRUE) * 0.8, max(y_vals, na.rm = TRUE) * 1.2)
    
    plot(fs_seq, res_test_stats[, k], type = "o", pch = 17, col = cols[k], lwd = 2,
         xlab = "Fs (Hz)", ylab = "Trace Stat", ylim = y_limits,
         main = paste0(sub_id, " - r<=", r_val), las = 1)
    
    # 绘制临界值红色虚线
    abline(h = current_cv, col = "red", lty = "dashed", lwd = 1.5)
    
    # 在图上显著位置标注具体的临界值数值，防止因量级问题看不清
    text(x = mean(range(fs_seq)), y = current_cv, 
         labels = sprintf("95%% CV = %.2f", current_cv), 
         pos = 3, col = "red", cex = 0.9, font = 2)
    
    grid(col = "gray", lty = "dotted")
  }
  dev.off()
  par(mfrow = c(1, 1)) 
  
  message(paste(">>> [完成] 图表已保存至 0_debug:", sub_id))
}

message("\n>>> ✅ 所有处理已完成。")
#1 环境与参数初始化====
library(R.matlab)
library(tsDyn)
library(ggplot2)
library(dplyr)
library(tidyr)

#1 路径与目录设置====
data_dir <- file.path(getwd(), "filter_data_downsample")
debug_dir <- file.path(getwd(), "0_debug")
sub_id <- "CP40" 

# 固定采样率，设置滞后项范围
fs_targ <- 400 
lags_to_test <- 0:5 
k <- 7 

results_tsdyn_eig <- list()

#2 数据导入与前处理 (移至循环外)====
file_path <- file.path(data_dir, paste0(sub_id, "_", fs_targ, "Hz.mat"))
mat_data <- readMat(file_path)

raw_data <- mat_data[["EEGdata"]] * 10^6
time_vec <- as.vector(mat_data[["EEGtimes"]]) 
events <- mat_data[["EEGevents"]]
ch_names <- trimws(apply(mat_data[["EEGchannels"]], 1, paste, collapse = ""))

event_codes <- as.numeric(events[, 3])
code_hundreds <- floor(event_codes / 100)
idx_all_trials <- which(code_hundreds %in% 1:2)

## 提取水平序列与拼接====
t_start <- -2.0
t_end <- 0.0
t_idx <- which(time_vec >= t_start & time_vec < t_end)

y_list <- list()
for (i in seq_along(idx_all_trials)) {
  orig_trial_idx <- idx_all_trials[i]
  epoch_slice <- raw_data[, t_idx, orig_trial_idx]
  y_list[[i]] <- t(epoch_slice) 
}

y <- do.call(rbind, y_list)
colnames(y) <- ch_names

rm(mat_data, raw_data); gc()

## 记录边界位置====
lens <- vapply(y_list, nrow, integer(1))
boundary_rows_in_y <- cumsum(lens)[-length(lens)]

#3 核心计算与模型拟合 (按滞后阶数循环)====
for (p in lags_to_test) {
  
  # 计算实际方程数
  num_equations <- nrow(y) - p - 1
  
  # 直接初始化外生变量矩阵
  exo_matrix <- matrix(0, nrow = num_equations, ncol = length(boundary_rows_in_y) * (p + 1))
  col_idx <- 1
  
  for (j in seq_along(boundary_rows_in_y)) {
    bdry_idx_in_y <- boundary_rows_in_y[j]
    
    for (lag_step in 0:p) {
      # 直接计算对应的方程行号并赋值为1
      eq_idx <- bdry_idx_in_y + 1 + lag_step - p - 1
      exo_matrix[eq_idx, col_idx] <- 1
      col_idx <- col_idx + 1
    }
  }
  
  # 直接命名，不再做任何列检查
  colnames(exo_matrix) <- paste0("bdry_", 1:ncol(exo_matrix))
  ## 拟合VECM模型与特征值提取
  fit <- lineVar(y, lag = p, model = "VECM", estim = "ML", r = k,
                 include = "const", exogen = exo_matrix)
  
  lambda_vals <- sort(Re(fit$model.specific$lambda), decreasing = TRUE)
  
  results_tsdyn_eig[[paste0("Lag_", p)]] <- data.frame(
    Lag = p,
    Rank_Index = 1:length(lambda_vals),
    Eigenvalue = as.numeric(lambda_vals)
  )
  
  cat("Lag", p, "completed.\n")
}



#4 结果合并与绘图====
df_tsdyn_eig <- do.call(rbind, results_tsdyn_eig)

p1 <- ggplot(df_tsdyn_eig, aes(x = Lag, y = Eigenvalue, color = factor(Rank_Index))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = lags_to_test) +
  labs(title = paste0("tsDyn: Eigenvalues vs Lag Order (Fs = ", fs_targ, "Hz)"),
       y = "Eigenvalue", x = "Lag Order (p)", color = "Rank Index") +
  theme_minimal()

print(p1)
ggsave(file.path(debug_dir, "tsDyn_Eigenvalues_vs_Lag.png"), p1, width = 8, height = 6)
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

freqs <- c(20, 25, 40, 50, 80, 100, 200, 400)
results_tsdyn_eig <- list()
results_tsdyn_test <- list()
k <- 7 

#2 核心计算与模型拟合====
for (fs_targ in freqs) {
  ##2.1 数据导入与条件筛选====
  file_path <- file.path(data_dir, paste0(sub_id, "_", fs_targ, "Hz.mat"))
  mat_data <- readMat(file_path)
  
  raw_data <- mat_data[["EEGdata"]] * 10^6
  time_vec <- as.vector(mat_data[["EEGtimes"]]) 
  events <- mat_data[["EEGevents"]]
  ch_names <- trimws(apply(mat_data[["EEGchannels"]], 1, paste, collapse = ""))
  
  event_codes <- as.numeric(events[, 3])
  code_hundreds <- floor(event_codes / 100)
  idx_all_trials <- which(code_hundreds %in% 1:2)
  
  ##2.2 提取水平序列与拼接====
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
  
  ##2.3 构造虚拟变量====
  lens <- vapply(y_list, nrow, integer(1))
  boundary_rows_in_y <- cumsum(lens)[-length(lens)]
  
  exo <- matrix(0, nrow = nrow(y) - 1, ncol = length(boundary_rows_in_y))
  for (j in seq_along(boundary_rows_in_y)) exo[boundary_rows_in_y[j], j] <- 1
  colnames(exo) <- paste0("bdry_", seq_along(boundary_rows_in_y))
  
  ##2.4 拟合VECM模型与统计量提取====
  fit <- lineVar(y, lag = 0, model = "VECM", estim = "ML", r = k,
                 include = "const", exogen = exo)
  
  lambda_vals <- sort(Re(fit$model.specific$lambda), decreasing = TRUE)
  
  results_tsdyn_eig[[paste0(fs_targ)]] <- data.frame(
    Fs = fs_targ,
    Rank_Index = 1:length(lambda_vals),
    Eigenvalue = as.numeric(lambda_vals)
  )
  
  p_dim <- length(lambda_vals)
  r_seq <- 0:(p_dim - 1)
  
  N_cpp <- nrow(y) - 2L * length(y_list)                
  eig_asc <- sort(lambda_vals, decreasing = FALSE)
  
  trace_stats <- -N_cpp * cumsum(log(1 - eig_asc))
  trace_stats <- sort(trace_stats, decreasing = TRUE)
  
  trace_cv_mat <- tsDyn:::gamma_doornik_all(
    q = c(0.90, 0.95, 0.99),
    nmp = p_dim:1,
    test = "H_lc",
    type = "trace"
  )
  cv95_vec <- as.numeric(trace_cv_mat[, "95%"])
  
  results_tsdyn_test[[paste0(fs_targ)]] <- data.frame(
    Fs = fs_targ,
    Rank_Hypothesis = factor(paste0("r<=", r_seq), levels = paste0("r<=", r_seq)),
    Stat_Value = as.numeric(trace_stats),
    CV_95 = cv95_vec 
  )
}

#3 结果合并与绘图====
df_tsdyn_eig <- do.call(rbind, results_tsdyn_eig)
df_tsdyn_test <- do.call(rbind, results_tsdyn_test)

p1 <- ggplot(df_tsdyn_eig, aes(x = Fs, y = Eigenvalue, color = factor(Rank_Index))) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = seq(20, 400, 40)) +
  labs(title = "tsDyn: Eigenvalues vs Sampling Rate",
       y = "Eigenvalue", x = "Sampling Rate (Hz)", color = "Rank Index") +
  theme_minimal()

cv_lab <- df_tsdyn_test %>%
  group_by(Rank_Hypothesis) %>%
  summarise(
    CV_95 = first(CV_95),
    Fs = max(Fs),
    .groups = "drop"
  )

p2 <- ggplot(df_tsdyn_test, aes(x = Fs)) +
  geom_line(aes(y = Stat_Value, color = "Statistic"), linewidth = 1) +
  geom_hline(
    data = cv_lab,
    aes(yintercept = CV_95, color = "Critical Value (95%)"),
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_text(
    data = cv_lab,
    aes(x = Fs, y = CV_95, label = sprintf("CV=%.1f", CV_95),
        color = "Critical Value (95%)"),
    hjust = 1.05, vjust = -0.4, size = 3,
    show.legend = FALSE
  ) +
  facet_wrap(~ Rank_Hypothesis, scales = "free_y") +
  scale_x_continuous(breaks = seq(20, 400, 80)) +
  labs(title = "tsDyn: Trace Statistics vs 95% CV",
       y = "Value", x = "Sampling Rate (Hz)") +
  scale_color_manual(values = c("Statistic" = "blue", "Critical Value (95%)" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)
print(p2)

ggsave(file.path(debug_dir, "tsDyn_Eigenvalues.png"), p1, width = 8, height = 6)
ggsave(file.path(debug_dir, "tsDyn_TestStats.png"), p2, width = 10, height = 8)
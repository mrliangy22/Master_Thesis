# 1 环境初始化与路径设置 ====
rm(list = ls()); gc()
library(R.matlab)
library(dplyr)
root_dir <- getwd()
data_dir <- file.path(root_dir, "filter_data")
processed_dir <- file.path(root_dir, "1_Preprocessing", "processed_matrices_stand")

if (!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)

# 获取文件列表
all_mat_files <- list.files(data_dir, pattern = "\\.mat$", full.names = TRUE)
message(paste(">>> [检测到文件数]", length(all_mat_files)))


# 2 定义核心辅助函数  ====
# 函数名：build_vecm_matrix
# 输入：
#   data_3d: 原始 EEG 数据 [Channels x Timepoints x Trials]
#   trial_indices: 需要提取的试次索引列表 (如 Correct Recall 的 trial ID)
#   time_vector: 原始时间轴向量
#   t_range: 截取的时间窗 (如 c(-2.0, 0.0))
# 输出：
#   list(Z0, Z1, trial_ids, n_trials)
#   其中 Z0 为差分矩阵 (dX), Z1 为滞后矩阵 (X_t-1)，行数=总样本量，列数=通道数

build_vecm_matrix <- function(data_3d, trial_indices, time_vector, t_range) {
  # [Step A] 定位时间窗索引
  # 在 time_vector 中找到 >= t_start 且 <= t_end 的所有下标
  t_start <- t_range[1]
  t_end   <- t_range[2]
  t_idx   <- which(time_vector >= t_start & time_vector < t_end)
  
  # 初始化列表，用于暂存每个试次处理后的矩阵
  Z0_list <- list()
  Z1_list <- list()
  id_list <- list() 
  
  # [Step B] 循环遍历每一个选中的试次
  for (i in seq_along(trial_indices)) {
    orig_trial_idx <- trial_indices[i]
    
    # 1. 切片 (Slicing): 从 3D 数组中取出 [Channels, Time] 的切片
    epoch_slice <- data_3d[, t_idx, orig_trial_idx]
    
    # 2. 转置 (Transpose): 变为 [Time, Channels]
    # 原理：VECM/VAR 模型要求每一行是一个时间点，每一列是一个变量(通道)
    # 这一步至关重要，否则后续 diff 算的是通道间的差分而不是时间差分
    epoch_mat <- t(epoch_slice)
    
    # 3. 构建 VECM 回归变量
    # VECM 公式: \Delta X_t = \Pi * X_{t-1} + ...
    # Z0 (左边项): 当前时刻 - 上一时刻 (dX)
    # Z1 (右边项): 上一时刻的值 (Lag 1)
    
    d_mat <- diff(epoch_mat)                # 计算一阶差分，行数会减少 1
    l_mat <- epoch_mat[1:(nrow(epoch_mat) - 1), ] # 取 1 到 N-1 行作为滞后项
    
    # 4. 存入列表
    Z0_list[[i]] <- d_mat
    Z1_list[[i]] <- l_mat
    # 记录该段数据属于第几个试次 (方便后续做 Cross-Validation 切分数据)
    id_list[[i]] <- rep(i, nrow(d_mat))
  }
  
  # [Step C] 堆叠 (Stacking)
  # 将 list 中所有试次的矩阵在垂直方向拼接 (rbind)
  # 最终形成一个巨大的长矩阵：[Total_Samples x Channels]
  return(list(
    Z0 = do.call(rbind, Z0_list),
    Z1 = do.call(rbind, Z1_list),
    trial_ids = unlist(id_list),
    n_trials = length(trial_indices)
  ))
}



# 3 批量处理循环 ====
for (input_file_path in all_mat_files) {
  
  file_name <- basename(input_file_path)
  sub_id <- sub("\\.mat$", "", file_name)
  message(paste0("\n>>> [Processing] ", sub_id))
  
  ## 3.1 导入数据与元数据 ====
  mat_data <- readMat(input_file_path)
  raw_data <- mat_data[["EEGdata"]]#*10^6 # 转换为微伏Ch, Time, Trials
  time_vec <- as.vector(mat_data[["EEGtimes"]]) 
  events   <- mat_data[["EEGevents"]]
  
  # 物理参数
  fs <- 400
  dt <- 1/fs
  
  # 整理通道信息 (名称、脑区、坐标)
  # apply paste 用于处理 MATLAB 导入后的字符矩阵格式问题
  ch_names   <- trimws(apply(mat_data[["EEGchannels"]], 1, paste, collapse = ""))
  ch_regions <- trimws(apply(mat_data[["EEGlabels"]], 1, paste, collapse = ""))
  ch_coords  <- mat_data[["EEGcoords"]]
  
  channel_info <- data.frame(
    index  = 1:length(ch_names),
    name   = ch_names,
    region = ch_regions,
    x      = ch_coords[, 1],
    y      = ch_coords[, 2],
    z      = ch_coords[, 3],
    stringsAsFactors = FALSE
  )
  
  rm(mat_data) # 读取完毕，立即释放大对象内存
  
  # 3.2 筛选试次索引 (Condition Filtering) ====
  # 逻辑：利用 Event Code 的百位数来判断条件
  # Code 格式：1xx(Recall Correct), 2xx(Manip Correct), 3xx/4xx(Incorrect), 5xx/6xx(Null)
  event_codes <- as.numeric(events[, 3])
  code_hundreds <- floor(event_codes / 100)
  
  # 分类 A: 仅提取正确试次 (用于后续网络差异分析)
  idx_recall_correct <- which(code_hundreds == 1)
  idx_manip_correct  <- which(code_hundreds == 2)
  
  # 分类 B: 提取所有试次 (用于基线期秩估计，最大化样本量)
  idx_all_trials     <- which(code_hundreds %in% 1:2) 
  
  message(sprintf("    -> Trials: Recall(Corr)=%d | Manip(Corr)=%d | All(Base)=%d", 
                  length(idx_recall_correct), length(idx_manip_correct), length(idx_all_trials)))
  
  # 3.3 构建 VECM 矩阵 (调用外置函数) ====
  message("    -> Building Matrices...")
  
  win_baseline <- c(-2.0, 0.0)
  win_delay    <- c( 2.0, 4.0)
  
  # (A) Baseline & Delay (Correct Only) - 分开构建
  # 这些数据将用于 Step 3 Lasso 分析，计算 specific condition 的网络
  res_base_rec  <- build_vecm_matrix(raw_data, idx_recall_correct, time_vec, win_baseline)
  res_base_man  <- build_vecm_matrix(raw_data, idx_manip_correct,  time_vec, win_baseline)
  res_delay_rec <- build_vecm_matrix(raw_data, idx_recall_correct, time_vec, win_delay)
  res_delay_man <- build_vecm_matrix(raw_data, idx_manip_correct,  time_vec, win_delay)
  
  # (B) Baseline (All Trials Combined) - 合并构建
  # 这些数据将用于 Step 2 秩估计 (Rank Estimation)，提供最稳健的统计基础
  res_base_all  <- build_vecm_matrix(raw_data, idx_all_trials,     time_vec, win_baseline)
  
  # 3.4 保存结果 ====
  final_data <- list(
    subject_id = sub_id,
    fs = fs,
    dt = dt,
    channels = channel_info,
    
    # 存储四个分条件数据块
    Baseline_Recall = res_base_rec,
    Baseline_Manip  = res_base_man,
    Delay_Recall    = res_delay_rec,
    Delay_Manip     = res_delay_man,
    
    # 存储全量基线数据块
    Baseline_All_Combined = res_base_all
  )
  
  save_path <- file.path(processed_dir, paste0(sub_id, ".rds"))
  saveRDS(final_data, file = save_path)
  message(paste0("    -> Saved: ", sub_id, ".rds"))
  # 清理循环内的大变量，防止内存溢出
  rm(raw_data, events, final_data, res_base_all, res_delay_rec)
  gc()
}

message("\n>>> ✅ 批量处理完成")
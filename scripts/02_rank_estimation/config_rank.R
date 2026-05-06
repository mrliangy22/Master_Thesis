#1 环境初始化====
rm(list = ls()); gc()
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(caret)
  library(ca4eeg)
#2 全局路径设置====
root_dir <- getwd()
base_output_dir <- file.path(root_dir, "2_Rank_Estimation","Rank_Results_stand")
input_data_dir <- file.path(root_dir, "1_Preprocessing", "processed_matrices_stand")
if (!dir.exists(base_output_dir)) {
  dir.create(base_output_dir, recursive = TRUE)
}
#3 全局变量定义====
# 完整的条件列表 (Superset)
conditions_to_run <- c(
  "Baseline_Recall", 
  "Baseline_Manip", 
  "Delay_Recall", 
  "Delay_Manip", 
  "Baseline_All_Combined"
)

# 扫描受试者 (自动读取文件名)
all_subjects_list <- gsub("\\.rds$", "", list.files(input_data_dir, pattern = "\\.rds$"))

# 3.2 运行范围控制 (Scope Control)
# A. 各步骤的目标受试者 (Target Subjects) - 默认全选
target_subjects_rsc   <- all_subjects_list
target_subjects_scree <- ("COG022")
target_subjects_mse   <-all_subjects_list
target_subjects_angle <- all_subjects_list
# B. 各步骤的目标条件 (Target Conditions) - 默认跑所有条件
# [新增功能] 这里允许你为每种方法单独定制要跑的条件
target_conditions_rsc   <- conditions_to_run
target_conditions_scree <- conditions_to_run
target_conditions_mse   <- c("Baseline_All_Combined")
target_conditions_angle <- c("Baseline_All_Combined")

# [自定义示例] 
# 如果你想让 MSE 只跑 "Baseline_All_Combined"，可以在这里覆盖默认值：
#target_conditions_mse <- c("Baseline_All_Combined")
#target_conditions_angle <- c("Baseline_All_Combined")
#4 批量初始化 (目录创建 & 路径预定义)====
# 建立全局路径字典 Subject_Files
Subject_Files <- list()

message(paste(">>> [Config] 正在初始化", length(all_subjects_list), "个受试者的路径配置..."))

for (sub_id in all_subjects_list) {
  
  # A. 确保目录存在
  sub_dir <- file.path(base_output_dir, sub_id)
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir, recursive = TRUE)
  }
  
  # B. 预定义所有数据文件路径 (不含Plot)
  Subject_Files[[sub_id]] <- list(
    # 目录
    Out_Dir    = sub_dir,
    
    # 输入数据
    Input_Data = file.path(input_data_dir, paste0(sub_id, ".rds")),
    
    # 输出结果 (RDS)
    Scree_File = file.path(sub_dir, paste0(sub_id, "_Data_Scree.rds")),
    RSC_File   = file.path(sub_dir, paste0(sub_id, "_Data_RSC.rds")),
    Angle_File = file.path(sub_dir, paste0(sub_id, "_Data_Angle.rds")),
    MSE_File   = file.path(sub_dir, paste0(sub_id, "_Data_MSE.rds"))
  )
}


#5 断点续传检测====
check_resume <- function(file_path) {
  if (file.exists(file_path)) {
    message("    -> [Resume] Found existing data. Loading...")
    return(readRDS(file_path))
  } else {
    message("    -> [Calculation] Starting new calculation...")
    return(NULL)
  }
}

message(">>> [Config] 配置加载完毕。")
# 6 生成秩估计参数配置文件 (Config Generation) ====
config_csv_path <- file.path(base_output_dir, "config_rank_params.csv")
if (!file.exists(config_csv_path)) {
message(">>> [Config] 正在扫描生成参数配置...")

# 初始化列表
config_list <- list()

# 1. 遍历受试者
for (sub_id in all_subjects_list) {
  
  # 直接读取 RDS (默认文件一定存在)
  sub_data <- readRDS(Subject_Files[[sub_id]]$Input_Data)
  
  # 2. 遍历 5 个条件
  for (cond_name in conditions_to_run) {
    
    # 直接提取数据 (默认条件一定存在)
    cond_data <- sub_data[[cond_name]]
    
    # 获取核心参数
    p_channels <- ncol(cond_data$Z0)
    n_trials   <- cond_data$n_trials
    
    # 3. 计算步长
    # 逻辑：floor((p - 1) / 9)
    step_val <- floor((p_channels - 1) / 9)
    # 数学防呆：如果算出来是 0 (通道极少)，强制为 1，防止报错
    if (step_val < 1) step_val <- 1 
    
    # 4. 存入列表
    config_list[[length(config_list) + 1]] <- data.frame(
      Subject     = sub_id,
      Condition   = cond_name,
      Channels    = p_channels,
      Trials      = n_trials,
      
      # MSE 参数
      MSE_Start   = 1,
      MSE_Max     = p_channels,
      MSE_Step    = step_val,
      
      # Angle 参数
      Angle_Start = 1,
      Angle_Max   = p_channels,
      Angle_Step  = step_val,
      
      stringsAsFactors = FALSE
    )
  }
}

# 5. 合并并保存 CSV
full_config_df <- do.call(rbind, config_list)
write.csv(full_config_df, config_csv_path, row.names = FALSE)

message(paste(">>> [Success] 配置文件已生成:", config_csv_path))
print(head(full_config_df)) }# 预览前几行确保正确
#7 秩的增量续传====
get_incremental_ranks <- function(file_path, target_ranks, cond_name) {
  
  # 1. 默认状态
  df_old <- NULL
  todo   <- target_ranks
  
  # 2. 唯一判断：文件若存在
  if (file.exists(file_path)) {
    df_old <- readRDS(file_path)
    
    # [修正点] 先筛选当前 Condition 的数据，再提取已跑过的 Rank
    # 如果该 Condition 还没跑过，ranks_done 为空，todo 保持全量
    ranks_done <- df_old$Rank[df_old$Condition == cond_name]
    todo       <- setdiff(target_ranks, ranks_done)
  }
  
  # 3. 返回结果
  return(list(data = df_old, to_run = todo))
}
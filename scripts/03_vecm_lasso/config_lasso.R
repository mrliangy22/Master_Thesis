#1 环境初始化====
rm(list = ls()); gc()
library(glmnet)
library(dplyr)
library(ca4eeg)
library(utils)

#2 全局路径设置====
root_dir <- getwd()
step_dir <- file.path(root_dir, "3_Lasso")
input_data_dir <- file.path(root_dir, "1_Preprocessing", "processed_matrices_stand") # 注意：这里我也更正为了 correct 目录，与之前步骤保持一致
base_output_dir <- file.path(step_dir, "Lasso_Results_stand")
rank_info_dir   <- file.path(root_dir, "2_Rank_Estimation", "Rank_Results_stand")

# 强制创建基础输出目录
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)


#3 全局变量定义====
conditions_to_run <- c(
  "Baseline_Recall", 
  "Baseline_Manip", 
  "Delay_Recall", 
  "Delay_Manip"
)

# Lasso 参数配置
lasso_props   <- c(0.25, 0.5, 0.75, 1)
n_penalty     <- 100
n_cv          <- 7
equal_penalty <- TRUE   # TRUE=全局统一惩罚, FALSE=逐通道优化
#k_folds <- 5 # 严格指定 10 折


#4 目标受试者筛选与秩导入 (核心逻辑修改)====
# 读取 config_rank_params.csv
rank_csv_path <- file.path(rank_info_dir, "config_rank_params.csv")

# 直球读取，不搞防御性检查
rank_df <- read.csv(rank_csv_path, stringsAsFactors = FALSE)

# 筛选逻辑：
# 1. selected_rank 不为空
# 2. 仅提取 Baseline_All_Combined 行 (确保每个受试者只有一个确定的秩用于所有条件)
valid_rows <- rank_df %>% 
  filter(!is.na(selected_rank) & as.character(selected_rank) != "") %>%
  filter(Condition == "Baseline_All_Combined")

# A. 目标受试者列表（手动指定）
all_subjects <- sub("\\.rds$", "", list.files(input_data_dir, pattern="\\.rds$", full.names=FALSE))
#target_subjects_lasso <- c()
#target_subjects_lasso <- all_subjects
#target_subjects_lasso <- c("BJH049","CP40"))
target_subjects_lasso <- c("BJH046")
#target_subjects_lasso <- c()
# B. 建立 "受试者 -> 秩" 的全局映射表
# 脚本后续直接用 Target_Rank_Map[["BJH041"]] 即可拿到该受试者的秩
Target_Rank_Map <- setNames(lapply(valid_rows$selected_rank, function(x) eval(parse(text=paste0("c(", x, ")")))), valid_rows$Subject)

message(paste0(">>> [Config] 参数表读取完毕。"))
message(paste0(">>> [Config] 有效目标 (在 Combined 行填写了 Rank): ", length(target_subjects_lasso)))


#5 批量初始化 (目录 & 路径)====
Subject_Files <- list()

for (sub_id in target_subjects_lasso) {
  
  # A. 确保目录存在
  sub_dir <- file.path(base_output_dir, sub_id)
  dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
  
  # B. 预定义路径
  Subject_Files[[sub_id]] <- list(
    Out_Dir    = sub_dir,
    Input_Data = file.path(input_data_dir, paste0(sub_id, ".rds"))
  )
  
  # C. 断点续传过滤：剔除已经跑过的秩
  all_ranks <- Target_Rank_Map[[sub_id]]
  pending_ranks <- c()
  for (r in all_ranks) {
    expected_file <- file.path(sub_dir, paste0(sub_id, "_Lasso_Results_r_", r, ".rds"))
    if (!file.exists(expected_file)) {
      pending_ranks <- c(pending_ranks, r)
    }
  }
  Target_Rank_Map[[sub_id]] <- pending_ranks
}


#6 公共辅助函数====

##6.1 断点续传检测====
check_resume <- function(file_path) {
  if (file.exists(file_path)) {
    message("    -> [Resume] Found existing data. Loading...")
    return(readRDS(file_path))
  } else {
    message("    -> [Calculation] Starting new calculation...")
    return(NULL)
  }
}

message(">>> [Config] Lasso 配置加载完毕。")
# 1 环境初始化与路径设置 ====
rm(list = ls()); gc()

library(ca4eeg)

# 设置数据路径
root_dir <- getwd()
processed_dir <- file.path(root_dir, "1_Preprocessing", "processed_matrices_stand")
target_sub <- "COG022"
rds_path <- file.path(processed_dir, paste0(target_sub, ".rds"))

# 2 读取数据与准备矩阵 ====
message(paste(">>> [DEBUG] 加载受试者数据:", target_sub))
sub_data <- readRDS(rds_path)

# 选择要测试的状态 (这里以 Baseline_Recall 为例)
# 对应 preprocess.R 第 127 行存储的数据块
target_cond <- "Baseline_Recall"
cond_data <- sub_data[[target_cond]]

Z0 <- cond_data$Z0 # 差分矩阵 dX
Z1 <- cond_data$Z1 # 滞后矩阵 X_{t-1}
dt <- sub_data$dt
p  <- ncol(Z0)

message(sprintf("    -> 数据规模: %d 个样本, %d 个通道", nrow(Z0), p))

# 3 运行模型获取底层矩阵 ====
# 注意：这里我们随便指定一个秩 r=1 即可，因为 S00 和 S11 的计算不依赖于最终选择的秩
message(">>> [DEBUG] 正在运行 Johansen 估计 (debug = TRUE)...")
model <- johansen(Y = Z0, Z = Z1, r = 1, dt = dt, intercept = TRUE, debug = TRUE)

# 4 核心数值稳定性诊断 ====
message("\n>>> [DEBUG] --- 底层矩阵数值稳定性诊断 ---")

kappa_threshold <- 1e10

# 提取特征值与条件数
kappa_S00 <- kappa(model$S00, exact = TRUE)
kappa_S11 <- kappa(model$S11, exact = TRUE)
min_eigen_S11 <- min(eigen(model$S11, symmetric = TRUE, only.values = TRUE)$values)

cat(sprintf("S00 (差分项协方差) 条件数: %.2e\n", kappa_S00))
cat(sprintf("S11 (滞后项协方差) 条件数: %.2e\n", kappa_S11))
cat(sprintf("S11 最小特征值: %.2e\n", min_eigen_S11))

# 诊断输出
if (kappa_S00 > kappa_threshold) {
  message(">> [警告] S00 矩阵存在严重多重共线性，求逆极度不稳定！")
} else {
  message(">> [正常] S00 矩阵条件数在安全范围内。")
}

if (kappa_S11 > kappa_threshold || min_eigen_S11 <= .Machine$double.eps) {
  message(">> [警告] S11 矩阵非正定或极度病态，特征分解存在崩溃风险！")
} else {
  message(">> [正常] S11 矩阵状态良好，满足正定要求。")
}

message("\n>>> 诊断完毕。")
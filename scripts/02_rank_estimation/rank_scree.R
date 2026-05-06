#1 加载公共配置====
source("2_Rank_Estimation/config_rank.R")

#2 脚本特有配置====
message(paste(">>> [Scree Analysis] 目标人数:", length(target_subjects_scree)))
csv_save_path <- file.path(base_output_dir, "Summary_Scree_Long_Format.csv")

# 如果每次都是重跑，建议初始化全局变量
scree_global_summary <- data.frame() 

#3 主循环执行 (Scree)====
for (sub_id in target_subjects_scree) {
  message(paste0("\n>>> [Scree] Processing: ", sub_id))
  
  # A. 获取路径
  paths <- Subject_Files[[sub_id]]
  
  # B. 直球读取数据
  subject_data <- readRDS(paths$Input_Data)
  dt <- subject_data$dt
  
  # C. 获取通道数
  p_channels <- nrow(subject_data$channels)
  sub_title_label <- sprintf("%s (p=%d)", sub_id, p_channels)
  
  # D. 断点续传 (本次暂时屏蔽，强制重算)
  df_scree_all <- NULL # <--- 强制重置，触发重算
  
  # E. 计算模块
  if (is.null(df_scree_all)) {
    list_scree <- list()
    
    # 循环配置文件指定的条件
    for (cond in target_conditions_scree) {
      
      # 直接提取
      data_block <- subject_data[[cond]]
      n_trials <- data_block$n_trials
      cond_display_label <- sprintf("%s (N=%d)", cond, n_trials)
      
      # [核心计算] Johansen r=1, debug=TRUE 才能输出 S 矩阵
      model_scree <- johansen(Y = data_block$Z0, Z = data_block$Z1, r = 1, dt = dt, intercept = TRUE, debug = TRUE)
      
      # 1. 提取所有需要的矩阵特征值，并按降序排列
      lambdas <- model_scree$lambda  # 广义特征值 lambda
      s00_vals <- sort(eigen(model_scree$S00)$values, decreasing = TRUE)
      s11_vals <- sort(eigen(model_scree$S11)$values, decreasing = TRUE)
      
      # 计算 S11 - S10 * S00^-1 * S01 并求特征值
      # 使用 chol 分解计算 S00 的逆，假定严格正定，无安全检查
      S00_inv <- chol2inv(chol(model_scree$S00))
      diff_matrix <- model_scree$S11 - (model_scree$S10 %*% S00_inv %*% model_scree$S01)
      diff_vals <- sort(eigen(diff_matrix)$values, decreasing = TRUE)
      
      # 2. 构造长格式数据 (行数 = p)
      res_scree <- data.frame(
        Subject = sub_id,
        Condition = cond_display_label, 
        Condition_Raw = cond,           
        Rank = 1:length(lambdas),       # 1 到 p
        Lambda = lambdas,               # 广义特征值
        S00_Eigenvalue = s00_vals,      # S00特征值
        S11_Eigenvalue = s11_vals,      # S11特征值
        Diff_Eigenvalue = diff_vals     # S11 - S10*S00^-1*S01 特征值
      )
      
      # 3. 这里的 RDS 依然可以存，作为备份
      list_scree[[cond]] <- res_scree
    }
    
    # 合并保存 RDS
    df_scree_all <- do.call(rbind, list_scree)
    saveRDS(df_scree_all, paths$Scree_File)
    message("    -> [Calculation] Done and saved.")
  }
  
  ##3.2 实时汇总 CSV (每人存一次)====
  # 添加到全局表
  scree_global_summary <- rbind(scree_global_summary, df_scree_all)
  
  # [立即写入] 覆盖更新整个 CSV 文件
  write.csv(scree_global_summary, csv_save_path, row.names = FALSE)
  message(paste("    -> [CSV] Updated global summary for", sub_id))
  
  ##3.3 绘图模块 (严格还原原始格式)====
  plot_path <- file.path(paths$Out_Dir, paste0(sub_id, "_Plot_Scree.png"))
  
  # [严格还原] 只画原始 Lambda，不加任何额外线条
  p_scree <- ggplot(df_scree_all, aes(x = Rank, y = Lambda)) +
    geom_line(color = "steelblue", size = 0.8) + geom_point(size = 1.2) +
    facet_wrap(~Condition, scales = "fixed", ncol = 3) +
    labs(title = paste("Scree Plot -", sub_title_label), y = "Eigenvalue") +
    theme_minimal() + theme(panel.grid.minor = element_blank())
  
  ggsave(plot_path, p_scree, width = 12, height = 6, bg = "white")
  
  # 释放内存
  rm(subject_data, df_scree_all, model_scree); gc()
}

message("\n>>> ✅ 所有受试者处理完毕")
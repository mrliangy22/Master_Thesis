#1 加载公共配置====
source("2_Rank_Estimation/config_rank.R")
library(cowplot) 

#2 脚本特有配置====
config_params_path <- file.path(base_output_dir, "config_rank_params.csv")
config_params_df <- read.csv(config_params_path, stringsAsFactors = FALSE)
k_list <- c(4, 5)   

message(paste(">>> [Angle Analysis] 目标人数:", length(target_subjects_angle)))

#3 主循环执行 (Angle)====
for (sub_id in target_subjects_angle) {
  
  message(paste0("\n>>> [Angle] Processing: ", sub_id))
  
  # A. 获取路径与数据
  paths <- Subject_Files[[sub_id]]
  subject_data <- readRDS(paths$Input_Data)
  dt <- subject_data$dt
  p_channels <- nrow(subject_data$channels)
  
  # 初始化总表
  if (file.exists(paths$Angle_File)) {
    df_angle_all <- readRDS(paths$Angle_File)
  } else {
    df_angle_all <- data.frame()
  }
  
  ### 3.1 遍历条件 ====
  for (cond in target_conditions_angle) {
    
    data_block <- subject_data[[cond]]
    Z0 <- data_block$Z0
    Z1 <- data_block$Z1
    trial_ids <- data_block$trial_ids
    n_trials <- length(unique(trial_ids))
    
    # 获取参数
    curr_param <- config_params_df[config_params_df$Subject == sub_id & 
                                config_params_df$Condition == cond, ]
    
    # 生成完整的目标 Rank 列表
    target_ranks <- seq(from = curr_param$Angle_Start, 
                        to = curr_param$Angle_Max, 
                        by = curr_param$Angle_Step)
    target_ranks <- target_ranks[target_ranks <= p_channels]
    
    cat(sprintf("      [Angle] Condition: %-25s (Target Ranks: %d)\n", cond, length(target_ranks)))
    
    # 3.2 K折循环 ====
    for (k in k_list) {
      
      # 1. [断点续传检查]
      # 检查当前 Condition 和 K 下，哪些 Rank 已经存在于数据中
      if (nrow(df_angle_all) > 0) {
        existing_ranks <- unique(df_angle_all$Rank[
          df_angle_all$Condition == cond & 
            df_angle_all$K_Fold == k
        ])
      } else {
        existing_ranks <- numeric(0)
      }
      
      # 计算剩余需要跑的 Rank
      todo_ranks <- setdiff(target_ranks, existing_ranks)
      
      if (length(todo_ranks) == 0) {
        # 该 K 下所有 Rank 都跑完了，直接跳过
        next
      }
      
      # 预先生成 Folds (保证对所有 Rank 切分一致)
      set.seed(456)
      folds_list <- caret::createFolds(1:n_trials, k = k, list = TRUE, returnTrain = FALSE)
      
      # 2. [循环结构调整] Rank (外) -> Fold (内)
      # 这样可以确保每跑完一个 Rank，其下所有 Folds 都完成，便于整块保存
      for (r_curr in todo_ranks) {
        
        # 临时存储当前 Rank 的所有 folds 结果
        rank_folds_results <- list()
        
        # --- Fold 循环 (必须全部跑完) ---
        for (f_idx in 1:k) {
          
          # 进度条
          cat(sprintf("\r          -> Calculating: K=%d | Rank %d | Fold %d/%d          ", k, r_curr, f_idx, k))
          flush.console()
          
          # 数据切分
          test_trial_indices <- folds_list[[f_idx]]
          test_rows  <- which(trial_ids %in% test_trial_indices)
          train_rows <- setdiff(1:nrow(Z0), test_rows)
          
          Z0_tr <- Z0[train_rows, , drop = FALSE]
          Z1_tr <- Z1[train_rows, , drop = FALSE]
          Z0_te <- Z0[test_rows, , drop = FALSE]
          Z1_te <- Z1[test_rows, , drop = FALSE]
          
          # 计算基准全秩模型 (Rank = p)
          model_full <- johansen(Z0_tr, Z1_tr, r = p_channels, dt = dt, intercept = TRUE)
          Pi_full <- model_full$alpha %*% t(model_full$beta)
          
          # 计算当前 Rank 模型
          model_r <- johansen(Z0_te, Z1_te, r = r_curr, dt = dt, intercept = TRUE)
          Pi_r <- model_r$alpha %*% t(model_r$beta)
          
          # 计算角度
          theta <- ca4eeg::matrix.angle(Pi_r, Pi_full)
          
          rank_folds_results[[length(rank_folds_results) + 1]] <- data.frame(
            Condition = cond,
            K_Fold = k,
            Fold_ID = f_idx,
            Rank = r_curr,
            Angle = theta
          )
        } # End Fold Loop
        
        # 3. [保存检查点]
        # 当前 Rank 的所有 Fold 都跑完了，合并并保存到总表
        df_chunk <- do.call(rbind, rank_folds_results)
        df_angle_all <- rbind(df_angle_all, df_chunk)
        saveRDS(df_angle_all, paths$Angle_File)
        
      } # End Rank Loop
    } # End K Loop
    cat("\n")
  } # End Condition Loop
  
  
  ### 3.3 统计与绘图模块 ====
  message("    -> [Plotting] Generating plots...")
  
  # 计算汇总统计 (只留 Mean)
  df_stats <- df_angle_all %>%
    group_by(Condition, K_Fold, Rank) %>%
    summarise(
      Angle_Mean = mean(Angle),
      .groups = "drop"
    ) 
  
  unique_conditions <- unique(df_angle_all$Condition)
  
  for (cond_name in unique_conditions) {
    
    sub_raw   <- df_angle_all %>% filter(Condition == cond_name)
    sub_stats <- df_stats %>% filter(Condition == cond_name)
    
    # [新增] 获取当前 Condition 的试次数量 N
    n_trials_current <- subject_data[[cond_name]]$n_trials
    
    # 构造带 N 的标题
    plot_title <- paste0("Mean Angle - ", cond_name, " (N=", n_trials_current, ")")
    
    # 图A: 均值
    p_mean <- ggplot(sub_stats, aes(x = Rank, y = Angle_Mean, color = as.factor(K_Fold))) +
      geom_line(size = 1) + geom_point(size = 2) +
      labs(title = plot_title, x = "Rank", y = "Mean Angle", color = "K") +
      theme_minimal() + theme(legend.position = "right") 
    
    # 图B: 细节
    p_detail <- ggplot(sub_raw, aes(x = Rank, y = Angle, group = interaction(K_Fold, Fold_ID), color = as.factor(K_Fold))) +
      geom_line(alpha = 0.5, size = 0.6) +
      facet_wrap(~K_Fold, labeller = label_both) +
      labs(title = "Detail Trails", x = "Rank", y = "Angle") +
      theme_minimal() + theme(legend.position = "none")
    
    final_plot <- plot_grid(p_mean, p_detail, ncol = 1, rel_heights = c(1, 1.2))
    
    safe_cond_name <- gsub(" ", "_", cond_name)
    plot_filename <- paste0(sub_id, "_Plot_Angle_", safe_cond_name, ".png")
    plot_full_path <- file.path(paths$Out_Dir, plot_filename)
    
    ggsave(plot_full_path, final_plot, width = 10, height = 8, bg = "white")
  }
  
  # 释放内存
  rm(subject_data, df_angle_all, df_stats); gc()
}

message("\n>>> ✅ Matrix Angle 任务完成")
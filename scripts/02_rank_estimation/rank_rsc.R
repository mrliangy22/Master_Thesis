#1 加载公共配置====
source("2_Rank_Estimation/config_rank.R")

#2 脚本特有配置====
rank_step_rsc <- 4  # [配置] 步长设置为 4
message(paste(">>> [RSC Analysis] 目标人数:", length(target_subjects_rsc)))

#3 主循环执行 (RSC)====
for (sub_id in target_subjects_rsc) {
  
  message(paste0("\n>>> [RSC] Processing: ", sub_id))
  
  # A. 获取预定义路径
  paths <- Subject_Files[[sub_id]]
  
  # B. 读取数据
  subject_data <- readRDS(paths$Input_Data)
  dt <- subject_data$dt
  
  # C. 获取通道数
  p_channels <- nrow(subject_data$channels)
  sub_title_label <- sprintf("%s (p=%d)", sub_id, p_channels)
  
  # D. 最简断点续传 (直接调用 config_rank.R 中的函数)
  # 逻辑：文件存在 -> 读入并跳过计算；文件不存在 -> 返回 NULL 并开始计算
  df_rsc_all <- check_resume(paths$RSC_File)
  
  # E. 计算模块 (仅当 check_resume 返回 NULL 时执行)
  if (is.null(df_rsc_all)) {
    
    list_rsc <- list()
    
    ### 3.1 遍历条件 (config中定义的) ====
    for (cond in target_conditions_rsc) {
      
      data_block <- subject_data[[cond]]
      
      Z0 <- data_block$Z0
      Z1 <- data_block$Z1
      N <- nrow(Z0)
      p <- ncol(Z0)
      n_trials <- data_block$n_trials
      
      # 定义目标 Rank 列表 (步长为 4)
      target_ranks <- sort(unique(c(seq(1, p, by = rank_step_rsc), p)))
      
      cat(sprintf("      [RSC] Condition: %-25s (Target Ranks: %d)\n", cond, length(target_ranks)))
      
      ## 步骤 A: 计算基准全秩模型 (r=p) ----
      model_full <- johansen(Y = Z0, Z = Z1, r = p, dt = dt, intercept = TRUE)
      rss_full   <- sum(diag(model_full$Omega)) * dt
      mu <- 4 / (N - p)
      theta_val <- mu * rss_full
      
      ## 步骤 B: 循环扫描秩 ----
      # 临时列表存储当前条件的 rank 结果
      list_cond_temp <- list()
      
      for (i in seq_along(target_ranks)) {
        r_curr <- target_ranks[i]
        
        # 进度条
        cat(sprintf("\r          -> Scanning Rank %d / %d", r_curr, p))
        flush.console()
        
        model_r <- johansen(Y = Z0, Z = Z1, r = r_curr, dt = dt, intercept = TRUE)
        rss_curr <- sum(diag(model_r$Omega)) * dt
        
        score <- rss_curr + (theta_val * r_curr)
        
        # 存入临时列表
        list_cond_temp[[i]] <- data.frame(
          Rank = r_curr, 
          Score = score,
          Condition = cond,
          N_Trials = n_trials
        )
      }
      cat("\n")
      
      # 将当前条件的计算结果整合
      list_rsc[[cond]] <- do.call(rbind, list_cond_temp)
      
    } # End Condition Loop
    
    # 所有条件跑完后，一次性合并并保存 RDS
    df_rsc_all <- do.call(rbind, list_rsc)
    saveRDS(df_rsc_all, paths$RSC_File)
    message("    -> [Calculation] Done and saved.")
  }
  
  ##3.3 绘图与保存 (基于内存中的 df_rsc_all)====
  if (!is.null(df_rsc_all) && nrow(df_rsc_all) > 0) {
    
    # 找出最小值点用于标记
    best_points <- df_rsc_all %>% group_by(Condition) %>% filter(Score == min(Score))
    
    # 标签拼接 (仅用于显示)
    labels_map <- setNames(
      paste0(best_points$Condition, " (N=", best_points$N_Trials, ")\nBest: ", best_points$Rank),
      best_points$Condition
    )
    
    plot_path <- file.path(paths$Out_Dir, paste0(sub_id, "_Plot_RSC.png"))
    
    p_rsc <- ggplot(df_rsc_all, aes(x = Rank, y = Score)) +
      geom_line(color = "black", size = 0.8) + 
      geom_point(data = best_points, color = "red", size = 3) +
      facet_wrap(
        ~Condition, 
        scales = "free_y", 
        ncol = 3,
        labeller = labeller(Condition = labels_map)
      ) +
      labs(title = paste("RSC Criterion -", sub_title_label), y = "RSC Score") +
      theme_minimal()
    
    ggsave(plot_path, p_rsc, width = 12, height = 6, bg = "white")
  }
  
  message(paste0("    -> [Finished] Processed and Plotted ", sub_id))
  
  rm(subject_data, df_rsc_all); gc()
}

message("\n>>> ✅ 所有受试者 RSC 处理完毕")
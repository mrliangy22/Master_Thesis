#1 环境初始化与配置====
rm(list = ls()); gc()
# 加载基础配置 (包含路径变量等)
source("2_Rank_Estimation/config_rank.R")
library(cowplot) # 用于多图拼接
library(ggplot2) # 确保绘图支持
# 直接读取参数配置文件 (默认文件存在)
config_file <- file.path(base_output_dir, "config_rank_params.csv")
rank_config <- read.csv(config_file, stringsAsFactors = FALSE)
message(paste(">>> [CV LOTO Analysis] 目标受试者人数:", length(target_subjects_mse)))
#2 主循环: 逐个受试者处理====
for (sub_id in target_subjects_mse) {
  
  message(paste0("\n>>> [LOTO] Processing Subject: ", sub_id))
  
  # 获取路径并读取 .rds 数据
  paths <- Subject_Files[[sub_id]]
  subject_data <- readRDS(paths$Input_Data)
  dt <- subject_data$dt
  
  # 遍历要跑的条件
  for (target_cond in target_conditions_mse) {
    
    ##2.1 直接加载数据====
    data_block <- subject_data[[target_cond]]
    p <- ncol(data_block$Z0)
    n_trials <- data_block$n_trials
    
    message(sprintf("    -> [Condition: %s] p = %d | N = %d", target_cond, p, n_trials))
    
    
    ##2.2 智能增量计算策略====
    conf_row <- rank_config[rank_config$Subject == sub_id & rank_config$Condition == target_cond, ]
    target_ranks <- seq(conf_row$MSE_Start, conf_row$MSE_Max, by = conf_row$MSE_Step)
    
    # 增量检查
    check_res <- get_incremental_ranks(paths$MSE_File, target_ranks, target_cond)
    
    # [关键修改] df_accumulated 用于实时存储当前已有的所有数据
    df_accumulated <- check_res$data   
    ranks_to_run   <- check_res$to_run 
    
    
    ##2.3 核心计算模块 (LOTO Cross-Validation)====
    if (length(ranks_to_run) > 0) {
      
      cat(sprintf("    -> [Calculation] 开始计算 %d 个新秩 (实时存档模式)...\n", length(ranks_to_run)))
      
      # 提取大矩阵到局部变量
      Z0 <- data_block$Z0
      Z1 <- data_block$Z1
      trial_ids <- data_block$trial_ids
      unique_trials <- unique(trial_ids)
      
      # --- 循环 A: 遍历每一个需要计算的秩 ---
      for (i in seq_along(ranks_to_run)) {
        r_curr <- ranks_to_run[i]
        
        # 初始化存储向量
        mse_test_vec  <- numeric(n_trials) 
        mse_train_vec <- numeric(n_trials) 
        ll_vec        <- numeric(n_trials) 
        
        # --- 循环 B: 留一试次交叉验证 (N 次) ---
        for (j in seq_along(unique_trials)) {
          
          # 实时进度条
          cat(sprintf("\r          -> Rank %d (%d/%d) | Trial %d/%d          ", 
                      r_curr, i, length(ranks_to_run), j, n_trials))
          if (j == n_trials) flush.console() 
          
          target_trial_id <- unique_trials[j]
          
          # 切分训练集与测试集
          idx_test  <- which(trial_ids == target_trial_id)
          idx_train <- setdiff(1:nrow(Z0), idx_test)
          
          Z0_tr <- Z0[idx_train, , drop = FALSE]
          Z1_tr <- Z1[idx_train, , drop = FALSE]
          Z0_te <- Z0[idx_test, , drop = FALSE]
          Z1_te <- Z1[idx_test, , drop = FALSE]
          
          # 1. 模型训练 (Johansen MLE)
          model <- johansen(Z0_tr, Z1_tr, r = r_curr, dt = dt, intercept = TRUE)
          
          # 提取参数
          Pi_cont <- model$alpha %*% t(model$beta)
          mu_cont <- model$Psi 
          
          # 2. 测试集预测
          term_drift_te <- matrix(mu_cont, nrow = nrow(Z0_te), ncol = p, byrow = TRUE)
          Z0_pred_te <- (Z1_te %*% t(Pi_cont) + term_drift_te) * dt
          
          res_te <- Z0_te - Z0_pred_te
          mse_test_vec[j] <- sum(res_te^2) / nrow(res_te)
          
          # 3. 训练集回测
          term_drift_tr <- matrix(mu_cont, nrow = nrow(Z0_tr), ncol = p, byrow = TRUE)
          Z0_pred_tr <- (Z1_tr %*% t(Pi_cont) + term_drift_tr) * dt
          
          res_tr <- Z0_tr - Z0_pred_tr
          mse_train_vec[j] <- sum(res_tr^2) / nrow(res_tr)
          
          # 4. 计算对数似然 (Log-Likelihood)
          Sigma_disc <- model$Omega * dt
          
          # [Math] 计算 Log-Determinant (防溢出)
          det_obj <- determinant(Sigma_disc, logarithm = TRUE)
          log_det_val <- as.numeric(det_obj$modulus)
          
          # [Math] 计算逆矩阵 (Cholesky 硬核质检)
          R_mat <- chol(Sigma_disc)   
          inv_Sigma <- chol2inv(R_mat) 
          
          # [Math] 计算马氏距离项 (Trace Trick)
          mahal_sum <- sum((res_te %*% inv_Sigma) * res_te)
          
          # [Math] 组合公式 68
          ll_val <- -0.5 * (log_det_val + mahal_sum / nrow(res_te))
          
          ll_vec[j] <- ll_val
          
        } # End Trial Loop
        
        # --- [新增] 单个 Rank 计算完毕，立即存档与绘图 ---
        
        # 1. 构造当前 Rank 的结果行
        df_single_rank <- data.frame(
          Rank = r_curr,
          MSE_Mean = mean(mse_test_vec, na.rm=TRUE),
          MSE_Train_Mean = mean(mse_train_vec, na.rm=TRUE),
          LogLik_Mean = mean(ll_vec, na.rm=TRUE),
          Condition = target_cond
        )
        
        # 2. 合并到累积数据框
        if (is.null(df_accumulated)) {
          df_accumulated <- df_single_rank
        } else {
          df_accumulated <- rbind(df_accumulated, df_single_rank)
        }
        
        # 3. 排序并立即覆盖保存 RDS
        df_accumulated <- df_accumulated[order(df_accumulated$Rank), ]
        saveRDS(df_accumulated, paths$MSE_File)
        
        # 4. 立即更新绘图 (Real-time Monitoring)
        # 只要有数据就画，覆盖旧图
        df_plot <- df_accumulated[df_accumulated$Condition == target_cond, ]
        
        if (nrow(df_plot) > 0) {
          plot_title <- paste("LOTO -", sub_id, "-", target_cond)
          
          p_mse_test <- ggplot(df_plot, aes(x = Rank, y = MSE_Mean)) +
            geom_line(color = "firebrick", size = 1) + 
            geom_point(color = "firebrick", size = 2) +
            labs(title = plot_title, y = "Test MSE (CV)") + 
            theme_minimal() +
            theme(axis.title.x = element_blank(), plot.margin = margin(b=2, unit="pt"))
          
          p_mse_train <- ggplot(df_plot, aes(x = Rank, y = MSE_Train_Mean)) +
            geom_line(color = "darkgreen", size = 1) + 
            geom_point(color = "darkgreen", size = 1.5) +
            labs(y = "Train MSE (Fit)") + 
            theme_minimal() +
            theme(axis.title.x = element_blank(), plot.title = element_blank(), 
                  plot.margin = margin(t=2, b=2, unit="pt"))
          
          p_ll <- ggplot(df_plot, aes(x = Rank, y = LogLik_Mean)) +
            geom_line(color = "steelblue", size = 1) + 
            geom_point(color = "steelblue", size = 2) +
            labs(y = "Test Log-Likelihood", x = "Rank") + 
            theme_minimal() +
            theme(plot.title = element_blank(), plot.margin = margin(t=2, unit="pt"))
          
          p_combined <- plot_grid(p_mse_test, p_mse_train, p_ll, 
                                  ncol = 1, align = "v", axis = "lr", 
                                  rel_heights = c(1.1, 1, 1))
          
          plot_path <- file.path(paths$Out_Dir, paste0(sub_id, "_", target_cond, "_LOTO.png"))
          ggsave(plot_path, p_combined, width = 8, height = 10, bg = "white")
        }
        
        cat(paste0(" [Saved]\n"))
        
      } # End Rank Loop (循环 A 结束)
      
      message("    -> [Finished] Condition Complete.")
      
    } else {
      message("    -> [Skip] 无需计算，直接读取已有数据。")
    }
    
    ##2.4 内存深度清理====
    suppressWarnings(rm(Z0, Z1, data_block, df_accumulated, df_plot, p_combined, df_single_rank))
    gc() 
    
  } # End Condition Loop
  
  rm(subject_data)
  gc()
  
} 

message("\n>>> ✅ LOTO CV 任务全部完成。")
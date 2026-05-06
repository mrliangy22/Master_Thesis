#1 加载公共配置与函数====
source("3_Lasso/config_lasso.R")
source("3_Lasso/func_lasso_equal.R")
library(ggplot2)

#2 主循环执行====
for (sub_id in target_subjects_lasso) {
  ranks_to_run <- Target_Rank_Map[[sub_id]]
  paths <- Subject_Files[[sub_id]]
  
  for (r_selected in ranks_to_run) {
    message(paste0("\n>>> [Lasso] Processing: ", sub_id, " | Rank = ", r_selected))
    current_lasso_file <- file.path(paths$Out_Dir, paste0(sub_id, "_Lasso_Results_r_", r_selected, ".rds"))
    lasso_results_all <- check_resume(current_lasso_file)
    
    if (is.null(lasso_results_all)) {
      subject_data <- readRDS(paths$Input_Data)
      dt <- subject_data$dt
      lasso_results_all <- list()
      
      ##2.1 预处理并汇总所有条件数据====
      cond_data_list <- list()
      p <- NULL
      
      for (cond_name in conditions_to_run) {
        target_data <- subject_data[[cond_name]]
        Y <- target_data$Z0
        Z <- target_data$Z1
        N <- nrow(Y) 
        p_cond <- ncol(Y)
        p <- p_cond
        
        meanY <- apply(Y, 2, mean)
        Ystd <- Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY))
        meanZ <- apply(Z, 2, mean)
        Zstd <- Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ))
        
        fit0 <- johansen(Y = Ystd, Z = Zstd, r = r_selected, dt = dt, intercept = FALSE)
        BETA.final <- fit0$beta 
        Z.r <- matrix(Zstd %*% BETA.final, nrow = N, ncol = r_selected)
        
        # 将该条件所需的核心数据封存入列表
        cond_data_list[[cond_name]] <- list(
          Y = Y, Z = Z, N = N, p = p, 
          meanY = meanY, meanZ = meanZ, 
          Ystd = Ystd, Zstd = Zstd, 
          BETA.final = BETA.final, Z.r = Z.r
        )
      }
      
      ##2.2 执行全局参数选择与模型拟合====
      use_lasso <- (r_selected > 1) && (p > 1)
      method_used <- ifelse(use_lasso, "Lasso", "OLS")
      
      if (!use_lasso) {
        ###2.2.1 OLS模式====
        fit_results <- list()
        for (cond_name in conditions_to_run) {
          c_data <- cond_data_list[[cond_name]]
          fit_ols <- lm(c_data$Ystd ~ c_data$Z.r - 1)
          fit_results[[cond_name]] <- list(
            alpha = t(coef(fit_ols)),
            penalty = rep(0, p),
            lasso_prop = rep(0, p)
          )
        }
      } else {
        ###2.2.2 Lasso模式====
        if (equal_penalty) {
          fit_results <- estimate_alpha_equal(cond_data_list, lasso_props, n_penalty, n_cv)
        } else {
          fit_results <- estimate_alpha_unequal(cond_data_list, lasso_props, n_penalty, n_cv)
        }
      }
      
      ##2.3 结果还原与指标计算====
      for (cond_name in conditions_to_run) {
        c_data <- cond_data_list[[cond_name]]
        c_fit <- fit_results[[cond_name]]
        
        Y <- c_data$Y
        Z <- c_data$Z
        N <- c_data$N
        BETA.final <- c_data$BETA.final
        ALPHA.Sparse <- c_fit$alpha
        meanY <- c_data$meanY
        meanZ <- c_data$meanZ
        
        PI <- ALPHA.Sparse %*% t(BETA.final)
        MU <- meanY - PI %*% as.matrix(meanZ, ncol = 1)
        
        res_mat <- Y - matrix(1, nrow = N, ncol = 1) %*% t(MU) - Z %*% t(PI)
        res0 <- Y - matrix(rep(meanY, each = N), nrow = N, ncol = p)
        R2 <- 1 - sum(res_mat^2) / sum(res0^2)
        OMEGA <- (t(res_mat) %*% res_mat) / N
        
        message(sprintf("       ... Final R2 (%s): %.4f (%s)", cond_name, R2, method_used))
        
        lasso_results_all[[cond_name]] <- list(
          condition  = cond_name,
          subject_id = sub_id,
          N = N, p = p, r = r_selected,
          alpha = ALPHA.Sparse / dt, 
          beta  = BETA.final,
          Pi    = PI / dt,
          Omega = OMEGA / dt,
          Psi   = MU / dt,
          R2    = R2,
          penalty = c_fit$penalty,
          lasso_prop = c_fit$lasso_prop,
          method_used = method_used
        )
      }
      
      saveRDS(lasso_results_all, current_lasso_file)
      message("    -> [Calculation] Done and saved.")
      rm(subject_data, cond_data_list, fit_results, lasso_results_all); gc()
    }
  }
  
  ##2.4 绘制并保存 R2 变化曲线====
  # 从本地读取该受试者所有已计算的秩结果以绘制完整的 R2 曲线
  all_res_files <- list.files(path = paths$Out_Dir, pattern = paste0(sub_id, "_Lasso_Results_r_\\d+\\.rds$"), full.names = TRUE)
  if (length(all_res_files) > 0) {
    r2_plot_data <- data.frame(Rank = integer(), Condition = character(), R2 = numeric())
    for (f in all_res_files) {
      res_data <- readRDS(f)
      for (cond in names(res_data)) {
        r2_plot_data <- rbind(r2_plot_data, data.frame(
          Rank = res_data[[cond]]$r,
          Condition = cond,
          R2 = res_data[[cond]]$R2
        ))
      }
    }
    
    if (nrow(r2_plot_data) > 0) {
      p_plot <- ggplot(r2_plot_data, aes(x = Rank, y = R2, color = Condition, group = Condition)) +
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        theme_bw() +
        labs(title = paste0("R2 Variation by Rank - Subject: ", sub_id), x = "Rank (r)", y = "R-squared") +
        scale_x_continuous(breaks = unique(r2_plot_data$Rank))
      
      plot_file <- file.path(paths$Out_Dir, paste0(sub_id, "_R2_Curve.png"))
      ggsave(plot_file, plot = p_plot, width = 8, height = 6)
      message(paste0("\n    -> [Plot] R2 curve generated and saved: ", plot_file))
    }
  }
}
message("\n>>> ✅ Lasso Analysis 批量任务完成")
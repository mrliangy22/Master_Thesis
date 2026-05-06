# 文件名: func_lasso.R
# 功能: Wrappers.R 中 Equal Penalty 模式 (全局统一惩罚)
# 优化: 综合多个条件的 CV (按样本量加权) 选择全局统一参数
# 修改: 移除了按试次 (Trial-based) 的手动 fold 生成，改由 glmnet 内部自动处理 CV 折数划分

estimate_alpha_equal <- function(cond_data_list, lasso_props, n_penalty, n_cv) {
  
  cond_names <- names(cond_data_list)
  p <- ncol(cond_data_list[[1]]$Ystd)
  r <- ncol(cond_data_list[[1]]$Z.r)
  
  total_N <- sum(sapply(cond_data_list, function(x) x$N))
  
  message(">>> [Branch] 进入 Equal Penalty 模式 (跨条件全局加权联合惩罚)...")
  
  #1 初始化最优参数记录容器====
  penalty.opt <- rep(NA, length(lasso_props))
  cv.opt <- rep(NA, length(lasso_props))
  penalty.seq.list <- list() 
  
  #2 第一层循环：控制Lasso比例====
  for (i in 1:length(lasso_props)) {
    lasso.prop.i <- lasso_props[i]
    cat(sprintf("\n>>> [Equal] Testing Lasso Prop %.2f (%d/%d)...\n", 
                lasso.prop.i, i, length(lasso_props)))
    
    ##2.1 确定跨条件的全局Lambda范围====
    penalty.max <- 0
    penalty.min <- NA 
    for (cond in cond_names) {
      Ystd <- cond_data_list[[cond]]$Ystd
      Z.r <- cond_data_list[[cond]]$Z.r
      
      for (i.r in 1:p) {
        fit_path <- glmnet(y = Ystd[, i.r], x = Z.r, intercept = FALSE, family = "gaussian", alpha = lasso.prop.i)
        
        penalty.max <- max(c(penalty.max, fit_path$lambda), na.rm = TRUE)
        penalty.min <- min(c(penalty.min, fit_path$lambda), na.rm = TRUE)
      }
    }
    
    ##2.2 生成并保存当前的Lambda序列====
    penalty.seq <- exp(seq(log(penalty.max), log(penalty.min), length = n_penalty))
    penalty.seq.list[[i]] <- penalty.seq
    
    ##2.3 为了稳定性的重复试验加权计算====
    cv_acc <- rep(0, n_penalty)
    
    ###2.3.1 执行测试 (Repeated CV)====
    # 外层的 i.cv 负责多次重复 CV 过程以取平均，增加稳健性
    for (i.cv in 1:n_cv) {
      
      ####2.3.1.1 遍历条件与通道执行CV====
      for (cond in cond_names) {
        Ystd <- cond_data_list[[cond]]$Ystd
        Z.r <- cond_data_list[[cond]]$Z.r
        weight <- cond_data_list[[cond]]$N / total_N
        
        for (i.r in 1:p) {
          cat(sprintf("\r    -> [CV] Rep: %02d/%02d | Cond: %-15s | Channel: %02d/%02d          ", 
                      i.cv, n_cv, cond, i.r, p))
          flush.console()
          
          # 【关键修改】：移除 foldid 参数，完全交由 glmnet 内部随机划分
          determine_lambda <- cv.glmnet(y = Ystd[, i.r], x = Z.r, intercept = FALSE, 
                                        family = "gaussian", lambda = penalty.seq, 
                                        alpha = lasso.prop.i)              
          
          current_cvm <- determine_lambda$cvm
          if(length(current_cvm) < n_penalty) {
            current_cvm <- c(current_cvm, rep(tail(current_cvm, 1), n_penalty - length(current_cvm)))
          }
          cv_acc <- cv_acc + (weight / p / n_cv) * current_cvm
        }
      }
    }
    cat("\n") 
    
    ##2.4 记录当前Prop下的最佳结果====
    penalty.opt[i] <- penalty.seq[which.min(cv_acc)]
    cv.opt[i] <- min(cv_acc)
  }
  
  #3 模型全局最优参数选择====
  best_idx <- which.min(cv.opt)
  lasso.prop.opt <- lasso_props[best_idx]
  penalty.final.val <- penalty.opt[best_idx]
  penalty.seq.opt <- penalty.seq.list[[best_idx]]
  
  message(sprintf(">>> 全局联合最佳参数: Prop = %.2f, Lambda = %.4f", lasso.prop.opt, penalty.final.val))
  message(">>> [Final Fit] Fitting models for all conditions...")
  
  #4 为每个条件建立最终模型====
  results_list <- list()
  
  for (cond in cond_names) {
    Ystd <- cond_data_list[[cond]]$Ystd
    Z.r <- cond_data_list[[cond]]$Z.r
    
    ALPHA.Sparse <- matrix(0, nrow = p, ncol = r)
    penalty.final <- rep(penalty.final.val, p)
    lasso.prop.final <- rep(lasso.prop.opt, p)
    
    for (i.r in 1:p) {
      cat(sprintf("\r    -> [Final Fit] Cond: %-15s | Channel: %02d/%02d          ", cond, i.r, p))
      flush.console()
      
      LASSOfinal <- glmnet(y = Ystd[, i.r], x = Z.r, intercept = FALSE,
                           lambda = penalty.seq.opt, 
                           family = "gaussian", alpha = lasso.prop.opt)
      
      # 提取系数时注意 glmnet 返回结果包含截距项(即使intercept=FALSE有时也会占位)，这里[-1]去除第一项
      coefs <- matrix(coef(LASSOfinal, s = penalty.final.val), nrow = 1)[-1]
      ALPHA.Sparse[i.r, ] <- coefs
    }
    cat("\n")
    
    results_list[[cond]] <- list(
      alpha = ALPHA.Sparse,
      penalty = penalty.final,
      lasso_prop = lasso.prop.final
    )
  }
  
  return(results_list)
}

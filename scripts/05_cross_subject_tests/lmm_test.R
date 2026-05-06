#1 线性混合模型(LMM)检验====
##1.1 环境与核心数据准备====
###1.1.1 加载必需的宏包====
# 若未安装需执行 install.packages("lmerTest")
library(lmerTest) 
library(dplyr)

###1.1.2 提取与清洗观测数据====
# 剔除 Contrast 数据，保留模型实际所需的连续观测因变量
lmm_data <- all_subjects_data %>%
  filter(Condition %in% c("Recall", "Manip")) %>%
  filter(!is.na(Interaction)) 

# 将 Condition 转换为因子，明确 Recall 为基准组
# 此时模型估计出的 ConditionManip 系数即代表 Manipulation 相比 Recall 的增量
lmm_data$Condition <- factor(lmm_data$Condition, levels = c("Recall", "Manip"))
lmm_data$Subject <- as.factor(lmm_data$Subject)

##1.2 批量执行模型拟合====
###1.2.1 初始化结果容器====
target_rois <- c("A1", "STG", "STS", "MTG", "ACC", "VLPFC", "DLPFC", "M1", "HPC", "PHC", "parietal")
lmm_results <- data.frame(
  From_Region = character(),
  To_Region = character(),
  N_Subjects = integer(),
  Estimate = numeric(),
  Std_Error = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

###1.2.2 遍历所有 ROI 连接矩阵====
for (from_roi in target_rois) {
  for (to_roi in target_rois) {
    
    ###1.2.3 过滤单条连线的有效被试====
    df_roi <- lmm_data %>%
      filter(From_Region == from_roi, To_Region == to_roi)
    
    # 确保被试在当前连线上同时拥有 Recall 和 Manip 的数据，形成有效配对
    valid_subjects <- df_roi %>%
      group_by(Subject) %>%
      summarise(n = n(), .groups = "drop") %>%
      filter(n == 2) %>%
      pull(Subject)
    
    df_roi <- df_roi %>% filter(Subject %in% valid_subjects)
    n_sub <- length(valid_subjects)
    
    # 样本量极小的连接无法支撑混合效应模型收敛，予以跳过
    if (n_sub < 3) {
      next
    }
    
    ###1.2.4 拟合单条连线的混合效应模型====
    # Interaction ~ Condition + (1 | Subject)
    # (1 | Subject) 捕获个体的基线交互强度差异（随机截距）
    model <- tryCatch({
      lmer(Interaction ~ Condition + (1 | Subject), data = df_roi)
    }, error = function(e) {
      return(NULL) 
    })
    
    ###1.2.5 提取固定效应参数====
    if (!is.null(model)) {
      coef_summary <- summary(model)$coefficients
      
      if ("ConditionManip" %in% rownames(coef_summary)) {
        lmm_results <- bind_rows(lmm_results, data.frame(
          From_Region = from_roi,
          To_Region = to_roi,
          N_Subjects = n_sub,
          Estimate = coef_summary["ConditionManip", "Estimate"],
          Std_Error = coef_summary["ConditionManip", "Std. Error"],
          t_value = coef_summary["ConditionManip", "t value"],
          p_value = coef_summary["ConditionManip", "Pr(>|t|)"] 
        ))
      }
    }
  }
}

##1.3 统计校正与结果输出====
###1.3.1 多重比较校正====
# 按原始 p 值排序，并计算 FDR
lmm_results <- lmm_results %>% arrange(p_value)
lmm_results$p_adj_fdr <- p.adjust(lmm_results$p_value, method = "fdr")

###1.3.2 打印与保存结果====
print(head(lmm_results, 15))
write.csv(lmm_results, file.path(root_dir, "LMM_Interaction_Results.csv"), row.names = FALSE)
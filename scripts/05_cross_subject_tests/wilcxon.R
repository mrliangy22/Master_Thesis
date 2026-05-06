#1 跨被试组级统计检验（支持多重比较校正选项）====
library(dplyr)
library(stats)

##1.1 参数与数据初始化====
# 设置检验方法: 
# "Wilcoxon" (非参数检验，即威尔科克森符号秩检验)
# "t-test" (参数检验，即单样本 t 检验)
test_method <- "t-test" 

# 设置校正模式: 
# "None" (不进行校正)
#"FDR" (不严格，采用 Benjamini-Hochberg 控制错误发现率)
# "Bonferroni" (最严格，控制族系整体错误率)
correction_mode <- "FDR" 

contrast_data <- all_subjects_data %>%
  filter(Condition == "Contrast") %>%
  filter(!is.na(Interaction))

##1.2 批量执行假设检验====
roi_pairs <- unique(paste(contrast_data$From_Region, contrast_data$To_Region, sep = "->"))
results_list <- list()

for (pair in roi_pairs) {
  regions <- unlist(strsplit(pair, "->"))
  from_r <- regions[1]
  to_r <- regions[2]
  
  pair_data <- contrast_data %>%
    filter(From_Region == from_r, To_Region == to_r)
  
  ###1.2.1 样本量检查与统计检验====
  # 剔除样本量过少的连线（t检验至少需要2个样本计算方差）
  if (nrow(pair_data) >= 2) { 
    
    # 根据设定的模式执行对应的检验
    if (test_method == "t-test") {
      test_res <- t.test(pair_data$Interaction, mu = 0)
    } else {
      test_res <- wilcox.test(pair_data$Interaction, mu = 0, exact = TRUE)
    }
    
    results_list[[pair]] <- data.frame(
      From_Region = from_r,
      To_Region = to_r,
      N_Subjects = nrow(pair_data),
      Mean_Contrast = mean(pair_data$Interaction),
      Median_Contrast = median(pair_data$Interaction),
      p_value = test_res$p.value
    )
  }
}

##1.3 合并结果并执行多重校正====
stat_results <- bind_rows(results_list)

if (correction_mode == "Bonferroni") {
  # 最严格的校正方式
  stat_results$p_adj <- p.adjust(stat_results$p_value, method = "bonferroni")
} else if (correction_mode == "FDR") {
  # 相对不严格的校正方式
  stat_results$p_adj <- p.adjust(stat_results$p_value, method = "BH")
} else {
  # 不进行多重比较校正
  stat_results$p_adj <- stat_results$p_value
}

##1.4 筛选显著边并输出====
# 按照原始未校正 p 值从小到大排列，并过滤 p_value < 0.8
significant_edges <- stat_results %>%
  filter(p_value < 0.8) %>%
  arrange(p_value)

print(significant_edges)
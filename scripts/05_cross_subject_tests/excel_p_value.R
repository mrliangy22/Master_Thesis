#1 跨被试组级统计检验（支持多重比较校正选项并导出Excel）====
library(dplyr)
library(stats)

# 如果未安装 writexl 包，请取消下一行的注释进行安装
# install.packages("writexl")
library(writexl)

##1.1 数据初始化====
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
    
    # 同时执行 T 检验 和 Wilcoxon 检验
    test_res_t <- t.test(pair_data$Interaction, mu = 0)
    
    # 注意：如果数据存在大量相同值(ties)，exact=TRUE 可能会触发警告。
    # 如果运行报错，可以将 exact = TRUE 改为 exact = FALSE
    test_res_w <- wilcox.test(pair_data$Interaction, mu = 0, exact = TRUE)
    
    results_list[[pair]] <- data.frame(
      From_Region = from_r,
      To_Region = to_r,
      N_Subjects = nrow(pair_data),
      Mean_Contrast = mean(pair_data$Interaction),
      Median_Contrast = median(pair_data$Interaction),
      p_value_ttest = test_res_t$p.value,
      p_value_wilcox = test_res_w$p.value
    )
  }
}

##1.3 合并结果并执行多重校正====
stat_results <- bind_rows(results_list)

# 对 t检验和 Wilcoxon 检验的结果分别进行 FDR (BH) 和 Bonferroni 校正
stat_results <- stat_results %>%
  mutate(
    p_adj_ttest_FDR = p.adjust(p_value_ttest, method = "BH"),
    p_adj_ttest_Bonferroni = p.adjust(p_value_ttest, method = "bonferroni"),
    p_adj_wilcox_FDR = p.adjust(p_value_wilcox, method = "BH"),
    p_adj_wilcox_Bonferroni = p.adjust(p_value_wilcox, method = "bonferroni")
  ) %>%
  # 默认按照 T 检验的未校正 P 值从小到大排序，方便查看
  arrange(p_value_ttest)

##1.4 将所有脑区对应的结果保存至 Excel 并打印显著边====
# 输出完整数据到 Excel 文件在当前工作目录
excel_filename <- "All_Regions_Statistics_Results.xlsx"
write_xlsx(stat_results, path = excel_filename)
cat("所有检验结果已成功保存至：", getwd(), "/", excel_filename, "\n\n", sep="")

# 在控制台打印初步筛选的结果（这里以 t 检验或 wilcox 检验 p<0.8 为例展示）
significant_edges <- stat_results %>%
  filter(p_value_ttest < 0.8 | p_value_wilcox < 0.8)

print(significant_edges)
# 环境初始化
rm(list = ls()); gc()
library(ggplot2)
library(dplyr)
library(tidyr)

# 路径设置
lasso_dir <- "3_Lasso/Lasso_Results_stand"
output_dir <- "4_Plotting/R2_Evaluation_stand"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 获取所有受试者文件夹
target_subjects <- list.dirs(lasso_dir, full.names = FALSE, recursive = FALSE)

r2_data <- data.frame()

# 循环读取并直接提取 R2
for (sub_id in target_subjects) {
  # 直接读取该受试者目录下的 rds 文件
  lasso_file <- list.files(file.path(lasso_dir, sub_id), pattern = "\\.rds$", full.names = TRUE)[1]
  lasso_res <- readRDS(lasso_file)
  
  # 直接提取你保存在列表里的 R2 标量值
  temp_df <- data.frame(
    Subject = sub_id,
    Baseline_Recall = lasso_res$Baseline_Recall$R2,
    Delay_Recall = lasso_res$Delay_Recall$R2,
    Baseline_Manip = lasso_res$Baseline_Manip$R2,
    Delay_Manip = lasso_res$Delay_Manip$R2
  )
  
  r2_data <- bind_rows(r2_data, temp_df)
}

# 转换为长格式以便分组画条形图
r2_long <- pivot_longer(r2_data, cols = -Subject, names_to = "Condition", values_to = "R2")

# 按总体 R2 均值给受试者排序，让直方图呈现递增阶梯状，便于观察群体分布
sub_order <- r2_data %>%
  mutate(Mean_R2 = (Baseline_Recall + Delay_Recall + Baseline_Manip + Delay_Manip) / 4) %>%
  arrange(Mean_R2) %>%
  pull(Subject)

r2_long$Subject <- factor(r2_long$Subject, levels = sub_order)
# 固定四个条件的显示顺序
r2_long$Condition <- factor(r2_long$Condition, levels = c("Baseline_Recall", "Delay_Recall", "Baseline_Manip", "Delay_Manip"))

# 画图
p <- ggplot(r2_long, aes(x = Subject, y = R2, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c("Baseline_Recall" = "#82B0D2", "Delay_Recall" = "#8E2043", 
                               "Baseline_Manip" = "#8CB369", "Delay_Manip" = "#F4A261")) +
  theme_bw() +
  labs(title = "R-squared across Subjects and Conditions", 
       x = "Subject", 
       y = expression(R^2)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank()
  )

# 导出图片和 CSV
ggsave(file.path(output_dir, "All_Subjects_R2_Barplot_high_rank.png"), plot = p, width = 14, height = 6, dpi = 300)
write.csv(r2_data, file.path(output_dir, "All_Subjects_R2_Values_high_rank.csv"), row.names = FALSE)
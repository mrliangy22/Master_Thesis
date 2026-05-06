library(dplyr)
library(readxl)

# 读取你刚生成的包含 ROI 的 Excel 文件
excel_data <- read_excel("1_Preprocessing/channel_label_summary_final.xlsx")

# 获取所有 rds 文件的路径
rds_files <- list.files("1_Preprocessing/processed_matrices_stand", pattern = "\\.rds$", full.names = TRUE)

# 循环更新每个 rds 文件
for (file in rds_files) {
  data <- readRDS(file)
  sub_id <- data$subject_id
  
  # 提取该受试者的 ChannelName 和 ROI
  sub_roi <- excel_data %>% 
    filter(SubjectID == sub_id) %>% 
    select(ChannelName, ROI)
  
  # 把 ROI 匹配合并到 rds 的 channels 信息中
  data$channels <- data$channels %>%
    left_join(sub_roi, by = c("name" = "ChannelName"))
  
  # 直接覆盖保存更新后的 rds
  saveRDS(data, file)
}

print("所有 rds 文件的 ROI 列已添加并覆盖保存完毕！")

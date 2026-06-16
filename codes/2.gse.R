# GSE281917清洗合并 ------------------------------------------------------------
library(dplyr)
# GSE 文件夹路径
gse_dir <- "./GSE"  # 替换成你自己的路径
# 找到所有样本子文件夹
sample_dirs <- list.dirs(gse_dir, recursive = FALSE)
# 初始化列表存储每个样本的数据
expr_list <- list()
for (sd in sample_dirs) {
  # 找到 txt 文件
  txt_files <- list.files(sd, pattern = "_log2CPM\\.txt$", full.names = TRUE)
  
  if (length(txt_files) == 0) next  # 没有文件就跳过
  
  f <- txt_files[1]  # 每个样本文件夹只有一个 txt
  
  # 样本名，从文件名提取 mucCRCXXX
  sample_name <- sub(".*_(mucCRC\\d+)_log2CPM\\.txt$", "\\1", basename(f))
  
  # 读取文件
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 确保列名为 probe_ID 和 sample
  colnames(df) <- c("probe_ID", sample_name)
  
  # 添加到列表
  expr_list[[sample_name]] <- df
}
# 合并所有样本（按 probe_ID）
expr_df <- Reduce(function(x, y) merge(x, y, by = "probe_ID", all = TRUE), expr_list)
# 查看结果
dim(expr_df)
head(expr_df)

expr_gse <- expr_df %>% 
  dplyr::rename(gene=probe_ID)

saveRDS(expr_gse,file = "expr_gse.RDS")

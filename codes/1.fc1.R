
# 1.fc1阈值清洗，统一命名 ---------------------------------------------------------
load("result_degs_fc1.Rda")
result_fc1 <- result
result_degs_fc1 <- result_degs
save(result_fc1,result_degs_fc1,
        file = "fc1/result_degs_fc1.Rda")

# 2.提取交集基因 ----------------------------------------------------------------
load("fc1/result_degs_fc1.Rda")

intersect_genes <- result_degs_fc1[[1]]
intersect_num <- c()
for(i in 2:length(result_degs_fc1)){
  
  intersect_genes <- intersect(
    intersect_genes,
    result_degs_fc1[[i]]
  )
  
  intersect_num[i] <- length(intersect_genes)
}
intersect_num

degs_freq_fc1 <- data.frame(
  intersection = c(2:24),
  genes = c(602,435, 383, 358 ,339, 314 ,302, 294, 286, 275, 273, 270, 269, 268,
            262 ,262 ,255, 254 ,254, 253, 253, 246 ,246
  )
) %>% 
  arrange(desc(intersection))

save(degs_freq_fc1,
     file = "fc1/degs_freq_fc1.Rda")


# 3.交集折线图 -----------------------------------------------------------------

library(ggview)
library(ggrepel)
library(tidyverse)
library(randomcoloR)
load("fc1/degs_freq_fc1.Rda")
# intersection,genes

Color <- randomColor(23)

p <- ggplot(degs_freq_fc1, aes(x = genes, y = reorder(as.character(intersection), 
                                                  -genes))) +
  geom_point(aes(color = as.character(intersection)), size = 5) +
  geom_point(aes(color = as.character(intersection)), size = 7, shape = 21, fill = NA) +
  geom_segment(aes(x = 0, xend = genes, y = reorder(as.character(intersection), -genes), 
                   yend = reorder(as.character(intersection), -genes), 
                   color = as.character(intersection)), 
               linewidth = 1) +
  scale_color_manual(values = Color) +
  geom_text(aes(label = genes), hjust = -1) +
  expand_limits(x = max(degs_freq_fc1$genes) * 1.1)+
  scale_color_manual(values = Color) +
  labs(y = "Number of intersecting sets", 
       x = "Number of genes")+
  theme_bw() + 
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0.5),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 14), # 调整 x 轴标题大小
        axis.title.y = element_text(size = 14), # 调整 y 轴标题大小
        panel.grid.major = element_blank(), # 去除主要网格线
        panel.grid.minor = element_blank()) + # 去除次要网格线 +
  canvas(width = 8, height = 6,dpi = 1000)
p

save_ggplot(p, "fc1/bar_chart_fc1.png")


# 4.提取基因 ------------------------------------------------------------------

library(tidyverse)
library(ComplexHeatmap)
library(randomcoloR)

load("fc1/result_degs_fc1.Rda")
load("fc1/degs_freq_fc1.Rda")
cols <- distinctColorPalette(24) #差异明显的60种

#生成组合矩阵
m <- make_comb_mat(result_degs_fc1)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
#输入集合的大小
set_size(m)
#组合集合的大小
comb_size(m)
#组合集合的度
comb_degree(m)
#转置组合集合
t(m)

# #提取核心基因
# genes_core <- extract_comb(m,"111111111111111111111111")

#提取非|核心集因
names <- names(result_degs_fc1[1:24])
genes_all <- data.frame(gene=as.character())
for (i in names){
  dat <- data.frame(gene=result_degs_fc1[[i]])
  genes_all <- rbind(genes_all,dat)
}
genes_all_fc1 <- genes_all %>%
  count(gene)

genes_core_fc1 <- genes_all %>% 
  filter(n==24)

genes_nocore_fc1 <- genes_all %>% 
  filter(n!=24)

save(genes_all_fc1,genes_core_fc1,genes_nocore_fc1,
     file = "fc1/genes_core_nocore_fc1.Rda")


# 5.非核心基因富集分析 -------------------------------------------------------------

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

load("genes_core_nocore_fc1.Rda")
# write.csv(genes_all,file = "genes_all.csv")
# write.csv(genes_core,file = "genes_core.csv")
# write.csv(genes_nocore,file = "genes_nocore.csv")
# GO
go_result <- enrichGO(gene = genes_nocore_fc1$gene,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)

# source("getGoTerm.R")# 公众号收藏
# GO_DATA<-get_GO_data("org.Hs.eg.db","ALL","SYMBOL")
# save(GO_DATA,file="GO_DATA.RData")
findGO<-function(pattern,method="key"){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  if(method=="key"){
    pathways=cbind(GO_DATA$PATHID2NAME[grep(pattern,GO_DATA$PATHID2NAME)])
  }else if(method=="gene"){
    pathways=cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
  }
  
  colnames(pathways)="pathway"
  
  if(length(pathways)==0){
    cat("No results!\n")
  } else{
    return(pathways)
  }
}
getGO<-function(ID){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME=names(GO_DATA$PATHID2EXTID)
  if(ID%in%allNAME){
    geneSet=GO_DATA$PATHID2EXTID[ID]
    names(geneSet)=GO_DATA$PATHID2NAME[ID]
    return(geneSet)
  }else{
    cat("No results!\n")
  }
}

load("GO_DATA.RData")#载入数据GO_DATA
# findGO("insulin")#寻找含有指定关键字的pathway name的pathway
# findGO("INS",method="gene")#寻找含有指定基因名的pathway
# getGO("GO:0045229")#获取指定GO ID的gene set

go_genes <- data.frame(gene=as.character())
go_ID <- go_data$ID
for(i in go_ID){
  p <- getGO(i)
  genes <- data.frame(gene=p[[1]])
  go_genes <- rbind(go_genes,genes) %>% 
    distinct(gene,.keep_all = T)
}

# top5 <- go_data %>%
#   group_by(ONTOLOGY) %>%
#   arrange(pvalue) %>%
#   slice_head(n = 5)
# df_top5 <- top5
# df_top5$ONTOLOGY <- factor(df_top5$ONTOLOGY, levels=c('CC', 'MF', "BP"))
# df_top5$Description <- factor(df_top5$Description, levels = rev(top5$Description))
# 
# mycol3 <- c('#6BA5CE', '#F5AA5F',"#8BA7BA")
# cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")
# ggplot(data = df_top5, aes(x = Count, y = Description, fill=ONTOLOGY)) +
#   geom_bar(width = 0.5,stat = 'identity') +
#   theme_classic() + 
#   scale_x_continuous(expand = c(0,0.5)) +
#   scale_fill_manual(values = alpha(mycol3, 0.66))+
#   geom_text(data = df_top5,
#             aes(x = 0.1, y = Description, label = Description),
#             size = 4.8,
#             hjust = 0)


# KEGG
diff_entrez<-bitr(
  genes_nocore_fc1$gene,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)

KEGG_result<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",#物种Homosapiens
                                         pvalueCutoff=0.05,#pvalue阈值
                                         qvalueCutoff=0.05,#qvalue阈值
                                         pAdjustMethod="BH",#p值矫正方法
                                         #"hochberg","hommel",
                                         #"bonferroni","BH",
                                         #"BY","fdr","none"
                                         minGSSize=10,#富集分析中考虑的最小基因集合大小
                                         maxGSSize=500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_result<-setReadable(KEGG_result,
                         OrgDb=org.Hs.eg.db,
                         keyType='ENTREZID')
#提取KEGG富集结果表格
KEGG_data <- as.data.frame(KEGG_result)

getKEGG<-function(ID){
  
  library("KEGGREST")
  
  gsList=list()
  for(xID in ID){
    
    gsInfo=keggGet(xID)[[1]]
    if(!is.null(gsInfo$GENE)){
      geneSetRaw=sapply(strsplit(gsInfo$GENE,";"),function(x)x[1])
      xgeneSet=list(geneSetRaw[seq(2,length(geneSetRaw),2)])
      NAME=sapply(strsplit(gsInfo$NAME,"-"),function(x)x[1])
      names(xgeneSet)=NAME
      gsList[NAME]=xgeneSet
    }else{
      cat("",xID,"No corresponding gene set in specific database.\n")
    }
  }
  return(gsList)
}
# getKEGG('hsa04930')
# 从KEGG富集结果中直接提取基因
KEGG_genes <- KEGG_data %>%
  dplyr::select(geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  distinct(geneID, .keep_all = TRUE) %>%
  rename(gene = geneID)

head(KEGG_genes)

genes_exclude <- rbind(KEGG_genes,go_genes) %>% 
  distinct(gene,.keep_all = T)

genes_exclude_fc1 <- genes_exclude
go_genes_fc1 <- go_genes
KEGG_genes_fc1 <- KEGG_genes
go_data_fc1 <- go_data
KEGG_data_fc1 <- KEGG_data

save(genes_exclude_fc1,go_genes_fc1,KEGG_genes_fc1,
     go_data_fc1,KEGG_data_fc1,
     file = "fc1/genes_exclude_fc1.Rda")


# 6.test1 -----------------------------------------------------------------
library(ComplexHeatmap)
library(ggvenn)
library(DESeq2)
library(tidyverse)
load("fc1/genes_core_nocore_fc1.Rda")
load("expr.Rda")
load("fc1/genes_exclude_fc1.Rda")
load("fc1/result_degs_fc1.Rda")

train_name <- names(result_degs_fc1)
name <-expr_srcc %>% 
  column_to_rownames(var = "gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  dplyr::select(name)
name <- name$name
test_name <- name[!name %in% train_name]


# test 1
expr_srcc <- expr_srcc %>% 
  arrange(gene) %>% 
  column_to_rownames(var = "gene")
expr_srcc <- 2^expr_srcc-1
expr_srcc <- round(expr_srcc)  # 将数据四舍五入为整数

expr_tumor <- expr_tumor %>% 
  arrange(gene) %>% 
  column_to_rownames(var = "gene")
expr_tumor <- 2^expr_tumor-1
expr_tumor<- round(expr_tumor)  # 将数据四舍五入为整数

srcc_sample <- test_name[[1]]

srcc_test1 <- expr_srcc %>% 
  dplyr::select(srcc_sample)
expr_diff1 <- cbind(expr_tumor,srcc_test1)
expr_diff1 <- round(expr_diff1 )
expr_diff1[expr_diff1<0] <- 0

expr_diff1_group <- expr_diff1 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::select(sample) %>% 
  mutate(group=ifelse(str_sub(sample,1,20)==srcc_sample,"SRCC","Tumor")) %>% 
  column_to_rownames(var = "sample")
# DIFF1
expr_diff1_group$group <- as.factor(expr_diff1_group$group)
colData <- data.frame(row.names = colnames(expr_diff1),                      
                      condition = expr_diff1_group$group)
dds <- DESeqDataSetFromMatrix(countData = expr_diff1,                              
                              colData = colData,                              
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition",rev(levels(expr_diff1_group$group))))
DEG <- res[order(res$pvalue),] %>% 
  as.data.frame()
k1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -1)
k2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > 1)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
# head(DEG)
DEG <- DEG %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene,change,pvalue,log2FoldChange)
DEG_genes <- DEG %>% 
  filter(change != "NOT")
test1_DEG_genes <- DEG_genes
genes_test1 <- DEG_genes$gene
genes_exclude_fc1 <- genes_exclude_fc1$gene

list <- list(test=genes_test1,exclude=genes_exclude_fc1)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_remain1 <- extract_comb(m,"10")

list <- list(test=test_remain1,core=genes_core_fc1$gene)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_core1 <- extract_comb(m,"11")

test_core1_fc1 <- test_core1
test_remain1_fc1 <- test_remain1
test1_DEG_genes_fc1 <- test1_DEG_genes
genes_test1_fc1 <- genes_test1

save(test_core1_fc1,test_remain1_fc1,test1_DEG_genes_fc1
     ,genes_test1_fc1,
     file = "fc1/test1_fc1.Rda")


# 7.test2 -------------------------------------------------------------
library(ComplexHeatmap)
library(ggvenn)
library(DESeq2)
library(tidyverse)
load("fc1/genes_core_nocore_fc1.Rda")
load("expr.Rda")
load("fc1/genes_exclude_fc1.Rda")
load("fc1/result_degs_fc1.Rda")

train_name <- names(result_degs_fc1)
name <-expr_srcc %>% 
  column_to_rownames(var = "gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  dplyr::select(name)
name <- name$name
test_name <- name[!name %in% train_name]


# test 1
expr_srcc <- expr_srcc %>% 
  arrange(gene) %>% 
  column_to_rownames(var = "gene")
expr_srcc <- 2^expr_srcc-1
expr_srcc <- round(expr_srcc)  # 将数据四舍五入为整数

expr_tumor <- expr_tumor %>% 
  arrange(gene) %>% 
  column_to_rownames(var = "gene")
expr_tumor <- 2^expr_tumor-1
expr_tumor<- round(expr_tumor)  # 将数据四舍五入为整数

srcc_sample <- test_name[[2]]

srcc_test2 <- expr_srcc %>% 
  dplyr::select(srcc_sample)
expr_diff2 <- cbind(expr_tumor,srcc_test2)
expr_diff2 <- round(expr_diff2 )
expr_diff2[expr_diff2<0] <- 0

expr_diff2_group <- expr_diff2 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::select(sample) %>% 
  mutate(group=ifelse(str_sub(sample,1,20)==srcc_sample,"SRCC","Tumor")) %>% 
  column_to_rownames(var = "sample")
# DIFF1
expr_diff2_group$group <- as.factor(expr_diff2_group$group)
colData <- data.frame(row.names = colnames(expr_diff2),                      
                      condition = expr_diff2_group$group)
dds <- DESeqDataSetFromMatrix(countData = expr_diff2,                              
                              colData = colData,                              
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition",rev(levels(expr_diff2_group$group))))
DEG <- res[order(res$pvalue),] %>% 
  as.data.frame()
k1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -1)
k2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > 1)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
# head(DEG)
DEG <- DEG %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene,change,pvalue,log2FoldChange)
DEG_genes <- DEG %>% 
  filter(change != "NOT")
test2_DEG_genes <- DEG_genes
genes_test2 <- DEG_genes$gene
genes_exclude_fc1 <- genes_exclude_fc1$gene

list <- list(test=genes_test2,exclude=genes_exclude_fc1)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_remain2 <- extract_comb(m,"10")

list <- list(test=test_remain2,core=genes_core_fc1$gene)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_core2 <- extract_comb(m,"11")

test_core2_fc1 <- test_core2
test_remain2_fc1 <- test_remain2
test2_DEG_genes_fc1 <- test2_DEG_genes
genes_test2_fc1 <- genes_test2

save(test_core2_fc1,test_remain2_fc1,test2_DEG_genes_fc1
     ,genes_test2_fc1,
     file = "fc1/test2_fc1.Rda")

# save(test_core1,file="genes.Rda")


# 8.交集 --------------------------------------------------------------------

load("fc1/test2_fc1.Rda")
load("fc1/test1_fc1.Rda")
load("fc1/genes_core_nocore_fc1.Rda")

library(tidyverse)
library(randomcoloR)
library(ggvenn)
library(ggview)

cols <- distinctColorPalette(3)
a <- list(sample_1=test_core1_fc1,
          sample_2=test_core2_fc1)
p <- ggvenn(a,c("sample_1","sample_2"),
            fill_color = cols,
            show_percentage=F,
            set_name_size=15,
            text_size = 15)+
  ggview::canvas(width = 16,height = 16,dpi = 1000)
p
save_ggplot(p, "fc1/venn_fc1.png")

# 9.富集分析 ------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
load("fc1/test2_fc1.Rda")
load("fc1/test1_fc1.Rda")

go_result <- enrichGO(gene = test_core1_fc1,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)

diff_entrez<-bitr(
  test_core1_fc1,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)

KEGG_result<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",#物种Homosapiens
                                         pvalueCutoff=0.05,#pvalue阈值
                                         qvalueCutoff=0.05,#qvalue阈值
                                         pAdjustMethod="BH",#p值矫正方法
                                         #"hochberg","hommel",
                                         #"bonferroni","BH",
                                         #"BY","fdr","none"
                                         minGSSize=5,#富集分析中考虑的最小基因集合大小
                                         maxGSSize=500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_result<-setReadable(KEGG_result,
                         OrgDb=org.Hs.eg.db,
                         keyType='ENTREZID')
#提取KEGG富集结果表格
KEGG_data <- as.data.frame(KEGG_result)

go_data_fc1 <- go_data
KEGG_data_fc1 <- KEGG_data

write_xlsx(go_data_fc1, path = "fc1/go_data_fc1.xlsx")
write_xlsx(KEGG_data_fc1, path = "fc1/KEGG_data_fc1.xlsx")

# devtools::install_github("dxsbiocc/gground")
library(gground)
library(ggprism)
library(randomcoloR)
library(ggview)

use_pathway <- group_by(go_data, ONTOLOGY) %>%
  rbind(KEGG_data%>%
          mutate(ONTOLOGY = 'KEGG')
  ) %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, 
                           levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

width <- 0.5 # 左侧分类标签和基因数量点图的宽度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1 # x 轴长度
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  ) # 左侧分类标签数据

pal <- distinctColorPalette()

p <- ggplot(use_pathway,aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) + # 绘制圆角柱状图
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 5
  ) + # 添加描述文本
  geom_text(
    aes(x = 0.1, label = geneID, colour = ONTOLOGY), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) + # 添加基因ID文本
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) + # 绘制基因数量点图
  geom_text(
    aes(x = -width, label = Count)
  ) + # 添加基因数量文本
  scale_size_continuous(name = 'Count', range = c(5, 16)) +# 设置点大小的比例尺
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +# 绘制分类标签矩形
  geom_text(
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes = FALSE
  ) + # 添加分类标签文本
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +# 绘制 x 轴线段
  labs(y = NULL) +
  scale_fill_manual(name = 'Category', values = pal) +# 设置填充颜色比例尺
  scale_colour_manual(values = pal) +# 设置线条颜色比例尺
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 0.5), 
    expand = expansion(c(0, 0))
  ) +# 设置 x 轴刻度
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    panel.background = element_blank()
  )+
  labs(x="-log10(p.adjust)")+
  canvas(width = 9, height = 7,dpi = 1000)
p

save_ggplot(p, "fc1/enrichment.png")


# 10.Rank Product-RankProd ------------------------------------------------
# BiocManager::install("RankProd")
library(RankProd)
library(tidyverse)

expr_srcc2 <- expr_srcc %>%
  column_to_rownames("gene")
expr_tumor2 <- expr_tumor %>%
  column_to_rownames("gene")

expr_all <- cbind(expr_srcc2, expr_tumor2)
expr_all <- as.matrix(expr_all)

cl <- c(
  rep(0, ncol(expr_srcc2)),
  rep(1, ncol(expr_tumor2))
)

RP.out <- RP(
  expr_all,
  cl = cl,
  logged = TRUE,
  gene.names = rownames(expr_all),
  rand = 123
)

up <- topGene(
  RP.out,
  cutoff = 0.05,
  method = "pfp",
  logged = TRUE
)$Table1 %>% 
  as.data.frame()

down <- topGene(
  RP.out,
  cutoff = 0.05,
  method = "pfp",
  logged = TRUE
)$Table2 %>% 
  as.data.frame()

names(up)
up <- up[up$pfp < 0.05, ]
down <- down[down$pfp < 0.05, ]
up <- up[up$`FC:(class1/class2)` >1, ]
down <- down[down$`FC:(class1/class2)` > 1, ]

up_gene <- data.frame(
    gene = rownames(up),
    change = "UP"
  )
down_gene <- data.frame(
  gene = rownames(down),
  change = "DOWN"
)
degs_rp <- rbind(up_gene, down_gene)
head(degs_rp)

save(degs_rp,RP.out,
     file = "fc1/degs_rp.rda")

# 11.RRA ------------------------------------------------------------------
# install.packages("RobustRankAggreg")
load("fc1/result_degs_fc1.Rda")
library(RobustRankAggreg)

deg_rank_list <- list()
for(i in 1:length(result_fc1)){
  
  res <- result_fc1[[i]]
  res <- res[!is.na(res$padj), ]
  res <- res[order(res$padj,
                   -abs(res$log2FoldChange)), ]
  deg_rank_list[[i]] <- res$gene
}
deg_rank_list

rra_result <- aggregateRanks(
  glist = deg_rank_list,
  method = "RRA"
)
rra_result

degs_rra <- rra_result[rra_result$Score < 0.05, ]

save(degs_rra,
     file = "fc1/degs_rra.rda")

# 12.三种方法取交集 --------------------------------------------------------------
load("fc1/degs_rra.rda")
load("fc1/degs_rp.rda")
load("fc1/test1_fc1.Rda")

library(tidyverse)
library(randomcoloR)
library(ggvenn)
library(ggview)

cols <- distinctColorPalette(3)
a <- list(RP=degs_rp$gene,
          RRA=degs_rra$Name,
          SM=test_core1_fc1)
p <- ggvenn(a,c("RP","RRA","SM"),
            fill_color = cols,
            show_percentage=F,
            set_name_size=15,
            text_size = 15)+
  ggview::canvas(width = 16,height = 16,dpi = 1000)
p

save_ggplot(p, "fc1/venn_rp_rra_sm.png")

# 13.rp进行富集分析 ---------------------------------------------------------
library(readxl)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
load("fc1/degs_rp.Rda")

go_result <- enrichGO(gene = degs_rp$gene,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)

diff_entrez<-bitr(
  degs_rp$gene,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)

KEGG_result<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",#物种Homosapiens
                                         pvalueCutoff=0.05,#pvalue阈值
                                         qvalueCutoff=0.05,#qvalue阈值
                                         pAdjustMethod="BH",#p值矫正方法
                                         #"hochberg","hommel",
                                         #"bonferroni","BH",
                                         #"BY","fdr","none"
                                         minGSSize=5,#富集分析中考虑的最小基因集合大小
                                         maxGSSize=500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_result<-setReadable(KEGG_result,
                         OrgDb=org.Hs.eg.db,
                         keyType='ENTREZID')
#提取KEGG富集结果表格
KEGG_data <- as.data.frame(KEGG_result)

go_data_fc1 <- go_data
KEGG_data_fc1 <- KEGG_data

write_xlsx(go_data_fc1, path = "fc1/go_data_rp.xlsx")
write_xlsx(KEGG_data_fc1, path = "fc1/KEGG_data_rp.xlsx")

# devtools::install_github("dxsbiocc/gground")
library(gground)
library(ggprism)
library(randomcoloR)
library(ggview)

go_top <- go_data %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()
kegg_top <- KEGG_data %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  mutate(ONTOLOGY = "KEGG")
kegg_top <- kegg_top[, colnames(go_top)]
use_pathway <- rbind(
  go_top,
  kegg_top
)
use_pathway <- use_pathway %>%
  mutate(
    ONTOLOGY = factor(
      ONTOLOGY,
      levels = rev(c("BP", "CC", "MF", "KEGG"))
    )
  )
use_pathway <- use_pathway %>%
  arrange(ONTOLOGY, p.adjust)
use_pathway <- use_pathway %>%
  mutate(
    Description = factor(
      Description,
      levels = Description
    )
  )
use_pathway <- use_pathway %>%
  rowid_to_column("index")
head(use_pathway)
table(use_pathway$ONTOLOGY)

width <- 0.5 # 左侧分类标签和基因数量点图的宽度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1 # x 轴长度
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  ) # 左侧分类标签数据

pal <- distinctColorPalette(4)

p <- ggplot(use_pathway,aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) + # 绘制圆角柱状图
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 5
  ) + # 添加描述文本
  geom_text(
    aes(x = 0.1, label = geneID, colour = ONTOLOGY), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) + # 添加基因ID文本
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) + # 绘制基因数量点图
  geom_text(
    aes(x = -width, label = Count)
  ) + # 添加基因数量文本
  scale_size_continuous(name = 'Count', range = c(5, 16)) +# 设置点大小的比例尺
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +# 绘制分类标签矩形
  geom_text(
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes = FALSE
  ) + # 添加分类标签文本
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +# 绘制 x 轴线段
  labs(y = NULL) +
  scale_fill_manual(name = 'Category', values = pal) +# 设置填充颜色比例尺
  scale_colour_manual(values = pal) +# 设置线条颜色比例尺
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 0.5), 
    expand = expansion(c(0, 0))
  ) +# 设置 x 轴刻度
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    panel.background = element_blank()
  )+
  labs(x="-log10(p.adjust)")+
  canvas(width = 9, height = 7,dpi = 1000)
p

save_ggplot(p, "fc1/enrichment_rp.png")

# 14.rra进行富集分析 ---------------------------------------------------------
library(readxl)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
load("fc1/degs_rra.Rda")

go_result <- enrichGO(gene = degs_rra$Name,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)

diff_entrez<-bitr(
  degs_rra$Name,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)

KEGG_result<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",#物种Homosapiens
                                         pvalueCutoff=0.05,#pvalue阈值
                                         qvalueCutoff=0.05,#qvalue阈值
                                         pAdjustMethod="BH",#p值矫正方法
                                         #"hochberg","hommel",
                                         #"bonferroni","BH",
                                         #"BY","fdr","none"
                                         minGSSize=5,#富集分析中考虑的最小基因集合大小
                                         maxGSSize=500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_result<-setReadable(KEGG_result,
                         OrgDb=org.Hs.eg.db,
                         keyType='ENTREZID')
#提取KEGG富集结果表格
KEGG_data <- as.data.frame(KEGG_result)

go_data_fc1 <- go_data
KEGG_data_fc1 <- KEGG_data

write_xlsx(go_data_fc1, path = "fc1/go_data_rra.xlsx")
write_xlsx(KEGG_data_fc1, path = "fc1/KEGG_data_rra.xlsx")

# devtools::install_github("dxsbiocc/gground")
library(gground)
library(ggprism)
library(randomcoloR)
library(ggview)

go_top <- go_data %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()
kegg_top <- KEGG_data %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  mutate(ONTOLOGY = "KEGG")
kegg_top <- kegg_top[, colnames(go_top)]
use_pathway <- rbind(
  go_top,
  kegg_top
)
use_pathway <- use_pathway %>%
  mutate(
    ONTOLOGY = factor(
      ONTOLOGY,
      levels = rev(c("BP", "CC", "MF", "KEGG"))
    )
  )
use_pathway <- use_pathway %>%
  arrange(ONTOLOGY, p.adjust)
use_pathway <- use_pathway %>%
  mutate(
    Description = factor(
      Description,
      levels = Description
    )
  )
use_pathway <- use_pathway %>%
  rowid_to_column("index")
head(use_pathway)
table(use_pathway$ONTOLOGY)

width <- 0.5 # 左侧分类标签和基因数量点图的宽度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1 # x 轴长度
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  ) # 左侧分类标签数据

pal <- distinctColorPalette(4)

p <- ggplot(use_pathway,aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) + # 绘制圆角柱状图
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 5
  ) + # 添加描述文本
  geom_text(
    aes(x = 0.1, label = geneID, colour = ONTOLOGY), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) + # 添加基因ID文本
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) + # 绘制基因数量点图
  geom_text(
    aes(x = -width, label = Count)
  ) + # 添加基因数量文本
  scale_size_continuous(name = 'Count', range = c(5, 16)) +# 设置点大小的比例尺
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +# 绘制分类标签矩形
  geom_text(
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes = FALSE
  ) + # 添加分类标签文本
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +# 绘制 x 轴线段
  labs(y = NULL) +
  scale_fill_manual(name = 'Category', values = pal) +# 设置填充颜色比例尺
  scale_colour_manual(values = pal) +# 设置线条颜色比例尺
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 0.5), 
    expand = expansion(c(0, 0))
  ) +# 设置 x 轴刻度
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    panel.background = element_blank()
  )+
  labs(x="-log10(p.adjust)")+
  canvas(width = 9, height = 7,dpi = 1000)
p

save_ggplot(p, "fc1/enrichment_rra.png")

# 15.gse折线图 ----------------------------------------------------------------
load("C:/Users/21118/Desktop/1.A novel biological function-based method for mining core genes in rare disease with limited cases/gse/degs_freq_gse.Rda")

library(ggview)
library(ggrepel)
library(tidyverse)
library(randomcoloR)
# intersection,genes

Color <- randomColor(7)

p <- ggplot(degs_freq, aes(x = genes, y = reorder(as.character(intersection), 
                                                      -genes))) +
  geom_point(aes(color = as.character(intersection)), size = 5) +
  geom_point(aes(color = as.character(intersection)), size = 7, shape = 21, fill = NA) +
  geom_segment(aes(x = 0, xend = genes, y = reorder(as.character(intersection), -genes), 
                   yend = reorder(as.character(intersection), -genes), 
                   color = as.character(intersection)), 
               linewidth = 1) +
  scale_color_manual(values = Color) +
  geom_text(aes(label = genes), hjust = -1) +
  expand_limits(x = max(degs_freq$genes) * 1.1)+
  scale_color_manual(values = Color) +
  labs(y = "Number of intersecting sets", 
       x = "Number of genes")+
  theme_bw() + 
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0.5),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 14), # 调整 x 轴标题大小
        axis.title.y = element_text(size = 14), # 调整 y 轴标题大小
        panel.grid.major = element_blank(), # 去除主要网格线
        panel.grid.minor = element_blank()) + # 去除次要网格线 +
  canvas(width = 8, height = 6,dpi = 1000)
p

save_ggplot(p, "gse/bar_chart_gse.png")


# 16. 提取基因 -----------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(randomcoloR)

load("gse/result_degs_gse.Rda")
load("gse/degs_freq_gse.Rda")
cols <- distinctColorPalette(24) #差异明显的60种

#生成组合矩阵
m <- make_comb_mat(result_degs)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
#输入集合的大小
set_size(m)
#组合集合的大小
comb_size(m)
#组合集合的度
comb_degree(m)
#转置组合集合
t(m)

# #提取核心基因
# genes_core <- extract_comb(m,"111111111111111111111111")

#提取非|核心集因
names <- names(result_degs[1:8])
genes_all <- data.frame(gene=as.character())
for (i in names){
  dat <- data.frame(gene=result_degs[[i]])
  genes_all <- rbind(genes_all,dat)
}
genes_all_gse <- genes_all %>%
  count(gene)

genes_core_gse <- genes_all_gse %>% 
  filter(n==8)

genes_nocore_gse <- genes_all_gse %>% 
  filter(n!=8)

save(genes_all_gse,genes_core_gse,genes_nocore_gse,
     file = "gse/genes_core_nocore_gse.Rda")

# 17.非核心基因富集分析 -------------------------------------------------------------

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

load("gse/genes_core_nocore_gse.Rda")
# write.csv(genes_all,file = "genes_all.csv")
# write.csv(genes_core,file = "genes_core.csv")
# write.csv(genes_nocore,file = "genes_nocore.csv")
# GO
go_result <- enrichGO(gene = genes_nocore_gse$gene,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)

# source("getGoTerm.R")# 公众号收藏
# GO_DATA<-get_GO_data("org.Hs.eg.db","ALL","SYMBOL")
# save(GO_DATA,file="GO_DATA.RData")
findGO<-function(pattern,method="key"){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  if(method=="key"){
    pathways=cbind(GO_DATA$PATHID2NAME[grep(pattern,GO_DATA$PATHID2NAME)])
  }else if(method=="gene"){
    pathways=cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
  }
  
  colnames(pathways)="pathway"
  
  if(length(pathways)==0){
    cat("No results!\n")
  } else{
    return(pathways)
  }
}
getGO<-function(ID){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME=names(GO_DATA$PATHID2EXTID)
  if(ID%in%allNAME){
    geneSet=GO_DATA$PATHID2EXTID[ID]
    names(geneSet)=GO_DATA$PATHID2NAME[ID]
    return(geneSet)
  }else{
    cat("No results!\n")
  }
}

load("GO_DATA.RData")#载入数据GO_DATA
# findGO("insulin")#寻找含有指定关键字的pathway name的pathway
# findGO("INS",method="gene")#寻找含有指定基因名的pathway
# getGO("GO:0045229")#获取指定GO ID的gene set

go_genes <- data.frame(gene=as.character())
go_ID <- go_data$ID
for(i in go_ID){
  p <- getGO(i)
  genes <- data.frame(gene=p[[1]])
  go_genes <- rbind(go_genes,genes) %>% 
    distinct(gene,.keep_all = T)
}

# top5 <- go_data %>%
#   group_by(ONTOLOGY) %>%
#   arrange(pvalue) %>%
#   slice_head(n = 5)
# df_top5 <- top5
# df_top5$ONTOLOGY <- factor(df_top5$ONTOLOGY, levels=c('CC', 'MF', "BP"))
# df_top5$Description <- factor(df_top5$Description, levels = rev(top5$Description))
# 
# mycol3 <- c('#6BA5CE', '#F5AA5F',"#8BA7BA")
# cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")
# ggplot(data = df_top5, aes(x = Count, y = Description, fill=ONTOLOGY)) +
#   geom_bar(width = 0.5,stat = 'identity') +
#   theme_classic() + 
#   scale_x_continuous(expand = c(0,0.5)) +
#   scale_fill_manual(values = alpha(mycol3, 0.66))+
#   geom_text(data = df_top5,
#             aes(x = 0.1, y = Description, label = Description),
#             size = 4.8,
#             hjust = 0)


# KEGG
diff_entrez<-bitr(
  genes_nocore_gse$gene,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)

KEGG_result<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",#物种Homosapiens
                                         pvalueCutoff=0.05,#pvalue阈值
                                         qvalueCutoff=0.05,#qvalue阈值
                                         pAdjustMethod="BH",#p值矫正方法
                                         #"hochberg","hommel",
                                         #"bonferroni","BH",
                                         #"BY","fdr","none"
                                         minGSSize=10,#富集分析中考虑的最小基因集合大小
                                         maxGSSize=500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_result<-setReadable(KEGG_result,
                         OrgDb=org.Hs.eg.db,
                         keyType='ENTREZID')
#提取KEGG富集结果表格
KEGG_data <- as.data.frame(KEGG_result)

getKEGG<-function(ID){
  
  library("KEGGREST")
  
  gsList=list()
  for(xID in ID){
    
    gsInfo=keggGet(xID)[[1]]
    if(!is.null(gsInfo$GENE)){
      geneSetRaw=sapply(strsplit(gsInfo$GENE,";"),function(x)x[1])
      xgeneSet=list(geneSetRaw[seq(2,length(geneSetRaw),2)])
      NAME=sapply(strsplit(gsInfo$NAME,"-"),function(x)x[1])
      names(xgeneSet)=NAME
      gsList[NAME]=xgeneSet
    }else{
      cat("",xID,"No corresponding gene set in specific database.\n")
    }
  }
  return(gsList)
}
# getKEGG('hsa04930')
# 从KEGG富集结果中直接提取基因
KEGG_genes <- KEGG_data %>%
  dplyr::select(geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  distinct(geneID, .keep_all = TRUE) %>%
  dplyr::rename(gene = geneID)

head(KEGG_genes)

genes_exclude <- rbind(KEGG_genes,go_genes) %>% 
  distinct(gene,.keep_all = T)

genes_exclude_gse <- genes_exclude
go_genes_gse <- go_genes
KEGG_genes_gse <- KEGG_genes
go_data_gse <- go_data
KEGG_data_gse <- KEGG_data

save(genes_exclude_gse,go_genes_gse,KEGG_genes_gse,
     go_data_gse,KEGG_data_gse,
     file = "gse/genes_exclude_gse.Rda")

# 18.test1 -----------------------------------------------------------------
library(ComplexHeatmap)
library(ggvenn)
library(DESeq2)
library(tidyverse)
load("expr.rda")
load("gse/genes_core_nocore_gse.Rda")
expr_gse <- readRDS("C:/Users/21118/Desktop/1.A novel biological function-based method for mining core genes in rare disease with limited cases/gse/expr_gse.RDS")
load("gse/genes_exclude_gse.Rda")
load("gse/result_degs_gse.Rda")

train_name <- names(result_degs)
name <-expr_gse %>% 
  column_to_rownames(var = "gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  dplyr::select(name)
name <- name$name
test_name <- name[!name %in% train_name]

# test 1
expr_gse <- expr_gse %>% 
  arrange(gene) %>% 
  column_to_rownames(var = "gene")
expr_gse <- 2^expr_gse-1
expr_gse <- round(expr_gse)  # 将数据四舍五入为整数

expr_tumor <- expr_tumor %>% 
  arrange(gene) %>% 
  column_to_rownames(var = "gene")
expr_tumor <- 2^expr_tumor-1
expr_tumor<- round(expr_tumor)  # 将数据四舍五入为整数

gse_sample <- test_name[[1]]

gse_test1 <- expr_gse %>% 
  dplyr::select(gse_sample)

common_gene <- intersect(
  rownames(expr_tumor),
  rownames(gse_test1)
)
expr_tumor <- expr_tumor[common_gene, ]
gse_test1 <- gse_test1[common_gene, , drop = FALSE]
expr_diff1 <- cbind(
  expr_tumor,
  gse_test1
)
expr_diff1 <- round(expr_diff1 )
expr_diff1[expr_diff1<0] <- 0

expr_diff1_group <- expr_diff1 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::select(sample) %>% 
  mutate(group=ifelse(str_sub(sample,1,9)==gse_sample,"GSE","Tumor")) %>% 
  column_to_rownames(var = "sample")
# DIFF1
expr_diff1_group$group <- as.factor(expr_diff1_group$group)
colData <- data.frame(row.names = colnames(expr_diff1),                      
                      condition = expr_diff1_group$group)
dds <- DESeqDataSetFromMatrix(countData = expr_diff1,                              
                              colData = colData,                              
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition",rev(levels(expr_diff1_group$group))))
DEG <- res[order(res$pvalue),] %>% 
  as.data.frame()
k1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -1)
k2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > 1)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
# head(DEG)
DEG <- DEG %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene,change,padj,log2FoldChange)
DEG_genes <- DEG %>% 
  filter(change != "NOT")
test1_DEG_genes <- DEG_genes
genes_test1 <- DEG_genes$gene
genes_exclude_gse <- genes_exclude_gse$gene

list <- list(test=genes_test1,exclude=genes_exclude_gse)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_remain1 <- extract_comb(m,"10")

list <- list(test=test_remain1,core=genes_core_gse$gene)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_core1 <- extract_comb(m,"11")

test_core1_gse <- test_core1
test_remain1_gse <- test_remain1
test1_DEG_genes_gse <- test1_DEG_genes
genes_test1_gse <- genes_test1

save(test_core1_gse,test_remain1_gse,test1_DEG_genes_gse
     ,genes_test1_gse,
     file = "gse/test1_gse.Rda")


# 19.test2 -------------------------------------------------------------
library(ComplexHeatmap)
library(ggvenn)
library(DESeq2)
library(tidyverse)
load("expr.rda")
load("gse/genes_core_nocore_gse.Rda")
expr_gse <- readRDS("C:/Users/21118/Desktop/1.A novel biological function-based method for mining core genes in rare disease with limited cases/gse/expr_gse.RDS")
load("gse/genes_exclude_gse.Rda")
load("gse/result_degs_gse.Rda")

train_name <- names(result_degs)
name <-expr_gse %>% 
  column_to_rownames(var = "gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  dplyr::select(name)
name <- name$name
test_name <- name[!name %in% train_name]

# test 2
expr_gse <- expr_gse %>% 
  arrange(gene) %>% 
  column_to_rownames(var = "gene")
expr_gse <- 2^expr_gse-1
expr_gse <- round(expr_gse)  # 将数据四舍五入为整数

expr_tumor <- expr_tumor %>% 
  arrange(gene) %>% 
  column_to_rownames(var = "gene")
expr_tumor <- 2^expr_tumor-1
expr_tumor<- round(expr_tumor)  # 将数据四舍五入为整数

gse_sample <- test_name[[2]]

gse_test2 <- expr_gse %>% 
  dplyr::select(gse_sample)

common_gene <- intersect(
  rownames(expr_tumor),
  rownames(gse_test2)
)
expr_tumor <- expr_tumor[common_gene, ]
gse_test2 <- gse_test2[common_gene, , drop = FALSE]
expr_diff1 <- cbind(
  expr_tumor,
  gse_test2
)
expr_diff1 <- round(expr_diff1 )
expr_diff1[expr_diff1<0] <- 0

expr_diff1_group <- expr_diff1 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::select(sample) %>% 
  mutate(group=ifelse(str_sub(sample,1,9)==gse_sample,"GSE","Tumor")) %>% 
  column_to_rownames(var = "sample")
# DIFF1
expr_diff1_group$group <- as.factor(expr_diff1_group$group)
colData <- data.frame(row.names = colnames(expr_diff1),                      
                      condition = expr_diff1_group$group)
dds <- DESeqDataSetFromMatrix(countData = expr_diff1,                              
                              colData = colData,                              
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition",rev(levels(expr_diff1_group$group))))
DEG <- res[order(res$pvalue),] %>% 
  as.data.frame()
k1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -1)
k2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > 1)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
# head(DEG)
DEG <- DEG %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene,change,padj,log2FoldChange)
DEG_genes <- DEG %>% 
  filter(change != "NOT")
test2_DEG_genes <- DEG_genes
genes_test2 <- DEG_genes$gene
genes_exclude_gse <- genes_exclude_gse$gene

list <- list(test=genes_test2,exclude=genes_exclude_gse)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_remain2 <- extract_comb(m,"10")

list <- list(test=test_remain2,core=genes_core_gse$gene)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_core2 <- extract_comb(m,"11")

test_core2_gse <- test_core2
test_remain2_gse <- test_remain2
test2_DEG_genes_gse <- test2_DEG_genes
genes_test2_gse <- genes_test2

save(test_core2_gse,test_remain2_gse,test2_DEG_genes_gse
     ,genes_test2_gse,
     file = "gse/test2_gse.Rda")


# 20.交集 --------------------------------------------------------------------

load("gse/test2_gse.Rda")
load("gse/test1_gse.Rda")
load("gse/genes_core_nocore_gse.Rda")

library(tidyverse)
library(randomcoloR)
library(ggvenn)
library(ggview)

cols <- distinctColorPalette(2)
a <- list(sample_1=test_core1_gse,
          sample_2=test_core2_gse)
p <- ggvenn(a,c("sample_1","sample_2"),
            fill_color = cols,
            show_percentage=F,
            set_name_size=15,
            text_size = 15)+
  ggview::canvas(width = 16,height = 16,dpi = 1000)
p
save_ggplot(p, "gse/venn_gse.png")

# 21.富集分析 ------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
load("gse/test2_gse.Rda")
load("gse/test1_gse.Rda")

go_result <- enrichGO(gene = test_core2_gse,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)

diff_entrez<-bitr(
  test_core2_gse,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)

KEGG_result<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",#物种Homosapiens
                                         pvalueCutoff=0.05,#pvalue阈值
                                         qvalueCutoff=0.05,#qvalue阈值
                                         pAdjustMethod="BH",#p值矫正方法
                                         #"hochberg","hommel",
                                         #"bonferroni","BH",
                                         #"BY","fdr","none"
                                         minGSSize=5,#富集分析中考虑的最小基因集合大小
                                         maxGSSize=500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_result<-setReadable(KEGG_result,
                         OrgDb=org.Hs.eg.db,
                         keyType='ENTREZID')
#提取KEGG富集结果表格
KEGG_data <- as.data.frame(KEGG_result)

go_data_gse <- go_data
KEGG_data_gse <- KEGG_data

write_xlsx(go_data_gse, path = "gse/go_data_gse.xlsx")
write_xlsx(KEGG_data_gse, path = "gse/KEGG_data_gse.xlsx")

# devtools::install_github("dxsbiocc/gground")
library(gground)
library(ggprism)
library(randomcoloR)
library(ggview)

use_pathway <- group_by(go_data, ONTOLOGY) %>%
  rbind(KEGG_data%>%
          mutate(ONTOLOGY = 'KEGG')
  ) %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, 
                           levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

width <- 0.5 # 左侧分类标签和基因数量点图的宽度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1 # x 轴长度
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  ) # 左侧分类标签数据

pal <- distinctColorPalette(1)

p <- ggplot(use_pathway,aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) + # 绘制圆角柱状图
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 5
  ) + # 添加描述文本
  geom_text(
    aes(x = 0.1, label = geneID, colour = ONTOLOGY), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) + # 添加基因ID文本
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) + # 绘制基因数量点图
  geom_text(
    aes(x = -width, label = Count)
  ) + # 添加基因数量文本
  scale_size_continuous(name = 'Count', range = c(5, 16)) +# 设置点大小的比例尺
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +# 绘制分类标签矩形
  geom_text(
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes = FALSE
  ) + # 添加分类标签文本
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +# 绘制 x 轴线段
  labs(y = NULL) +
  scale_fill_manual(name = 'Category', values = pal) +# 设置填充颜色比例尺
  scale_colour_manual(values = pal) +# 设置线条颜色比例尺
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 0.5), 
    expand = expansion(c(0, 0))
  ) +# 设置 x 轴刻度
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    panel.background = element_blank()
  )+
  labs(x="-log10(p.adjust)")+
  canvas(width = 9, height = 7,dpi = 1000)
p

save_ggplot(p, "gse/enrichment_gse.png")


# 22.初始差异基因富集分析图tcga的数据 -----------------------------------------------------------

library(readxl)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
load("fc1/test1_fc1.Rda")

go_result <- enrichGO(gene = genes_test1_fc1,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)

diff_entrez<-bitr(
  genes_test1_fc1,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)

KEGG_result<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",#物种Homosapiens
                                         pvalueCutoff=0.05,#pvalue阈值
                                         qvalueCutoff=0.05,#qvalue阈值
                                         pAdjustMethod="BH",#p值矫正方法
                                         #"hochberg","hommel",
                                         #"bonferroni","BH",
                                         #"BY","fdr","none"
                                         minGSSize=5,#富集分析中考虑的最小基因集合大小
                                         maxGSSize=500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_result<-setReadable(KEGG_result,
                         OrgDb=org.Hs.eg.db,
                         keyType='ENTREZID')
#提取KEGG富集结果表格
KEGG_data <- as.data.frame(KEGG_result)

go_data_gse <- go_data
KEGG_data_gse <- KEGG_data

# write_xlsx(go_data_gse, path = "gse/go_data_gse.xlsx")
# write_xlsx(KEGG_data_gse, path = "gse/KEGG_data_gse.xlsx")

# devtools::install_github("dxsbiocc/gground")
library(gground)
library(ggprism)
library(randomcoloR)
library(ggview)

use_pathway <- group_by(go_data, ONTOLOGY) %>%
  rbind(KEGG_data%>%
          mutate(ONTOLOGY = 'KEGG')
  ) %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, 
                           levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

width <- 0.5 # 左侧分类标签和基因数量点图的宽度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1 # x 轴长度
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  ) # 左侧分类标签数据

pal <- distinctColorPalette(13)

p <- ggplot(use_pathway,aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) + # 绘制圆角柱状图
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 5
  ) + # 添加描述文本
  geom_text(
    aes(x = 0.1, label = geneID, colour = ONTOLOGY), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) + # 添加基因ID文本
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) + # 绘制基因数量点图
  geom_text(
    aes(x = -width, label = Count)
  ) + # 添加基因数量文本
  scale_size_continuous(name = 'Count', range = c(5, 16)) +# 设置点大小的比例尺
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +# 绘制分类标签矩形
  geom_text(
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes = FALSE
  ) + # 添加分类标签文本
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +# 绘制 x 轴线段
  labs(y = NULL) +
  scale_fill_manual(name = 'Category', values = pal) +# 设置填充颜色比例尺
  scale_colour_manual(values = pal) +# 设置线条颜色比例尺
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 0.5), 
    expand = expansion(c(0, 0))
  ) +# 设置 x 轴刻度
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    panel.background = element_blank()
  )+
  labs(x="-log10(p.adjust)")+
  canvas(width = 9, height = 7,dpi = 1000)
p

save_ggplot(p, "fc1/enrichment_差异基因.png")

# 23.去掉非核心基因但是没和初始核心基因取交集富集分析图tcga的数据 -----------------------------------------------------------

library(readxl)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
load("fc1/test1_fc1.Rda")

go_result <- enrichGO(gene = test_remain1_fc1,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)

diff_entrez<-bitr(
  test_remain1_fc1,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)

KEGG_result<-clusterProfiler::enrichKEGG(gene=diff_entrez$ENTREZID,
                                         organism="hsa",#物种Homosapiens
                                         pvalueCutoff=0.05,#pvalue阈值
                                         qvalueCutoff=0.05,#qvalue阈值
                                         pAdjustMethod="BH",#p值矫正方法
                                         #"hochberg","hommel",
                                         #"bonferroni","BH",
                                         #"BY","fdr","none"
                                         minGSSize=5,#富集分析中考虑的最小基因集合大小
                                         maxGSSize=500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_result<-setReadable(KEGG_result,
                         OrgDb=org.Hs.eg.db,
                         keyType='ENTREZID')
#提取KEGG富集结果表格
KEGG_data <- as.data.frame(KEGG_result)

go_data_gse <- go_data
KEGG_data_gse <- KEGG_data

# write_xlsx(go_data_gse, path = "gse/go_data_gse.xlsx")
# write_xlsx(KEGG_data_gse, path = "gse/KEGG_data_gse.xlsx")

# devtools::install_github("dxsbiocc/gground")
library(gground)
library(ggprism)
library(randomcoloR)
library(ggview)

use_pathway <- group_by(go_data, ONTOLOGY) %>%
  rbind(KEGG_data%>%
          mutate(ONTOLOGY = 'KEGG')
  ) %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, 
                           levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

width <- 0.5 # 左侧分类标签和基因数量点图的宽度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1 # x 轴长度
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  ) # 左侧分类标签数据

pal <- distinctColorPalette(13)

p <- ggplot(use_pathway,aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) + # 绘制圆角柱状图
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 5
  ) + # 添加描述文本
  geom_text(
    aes(x = 0.1, label = geneID, colour = ONTOLOGY), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) + # 添加基因ID文本
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) + # 绘制基因数量点图
  geom_text(
    aes(x = -width, label = Count)
  ) + # 添加基因数量文本
  scale_size_continuous(name = 'Count', range = c(5, 18)) +# 设置点大小的比例尺
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +# 绘制分类标签矩形
  geom_text(
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes = FALSE
  ) + # 添加分类标签文本
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +# 绘制 x 轴线段
  labs(y = NULL) +
  scale_fill_manual(name = 'Category', values = pal) +# 设置填充颜色比例尺
  scale_colour_manual(values = pal) +# 设置线条颜色比例尺
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 0.5), 
    expand = expansion(c(0, 0))
  ) +# 设置 x 轴刻度
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    panel.background = element_blank()
  )+
  labs(x="-log10(p.adjust)")+
  canvas(width = 9, height = 7,dpi = 1000)
p

save_ggplot(p, "fc1/enrichment_初始核心基因取交集.png")


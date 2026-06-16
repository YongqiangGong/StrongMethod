# 1.packages --------------------------------------------------------------

library(tidyverse)
library(DESeq2)#差异分析
library(sva)#批次效应
library(data.table)
library(ggvenn)

# 2.renmin datasets -------------------------------------------------------

renmin_data <- epxr_renmin %>% 
  dplyr::rename(gene=gene_symbol) %>% 
  aggregate(.~gene,mean) %>% 
  column_to_rownames(var="gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="sample") %>% 
  mutate(sample=str_sub(sample,12,30)) %>% 
  filter(str_sub(sample,6,6)!="P") %>% 
  filter(str_sub(sample,6,6)!="N")

sample_name <- str_c("TUMC_T",83:108)

expr_srcc <- renmin_data %>% 
  filter(sample %in% sample_name)
expr_srcc <- expr_srcc %>% 
  column_to_rownames(var = "sample")
expr_srcc <- log2(expr_srcc+1)

expr_tumor <- renmin_data %>% 
  filter(!sample %in% sample_name)

# 3.TCGA and article datasets --------------------------------------------------------------

# TCGA

expr_tcga <- expr_tcga %>% 
  rownames_to_column(var = "gene") %>% 
  aggregate(.~gene,mean) %>%
  column_to_rownames(var="gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  filter(str_sub(sample,14,15)=="01") %>% 
  column_to_rownames(var = "sample") %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "gene")
  
# Article
expr_article <- fread("expr.txt")%>% 
  column_to_rownames(var="SYMBOL") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  filter(str_sub(sample,14,14)=="T") %>% 
  column_to_rownames(var = "sample") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  aggregate(.~gene,mean)

expr_article <- expr_article %>% 
  column_to_rownames(var="gene")
expr_article <- log2(expr_article+1)
expr_article <- expr_article %>% 
  rownames_to_column(var="gene")

# 批次效应 tcga and article

batch<-c(rep("TCGA",614),rep("Article",1063))#批次信息

all_expr <- expr_tcga %>% 
  inner_join(.,expr_article,by="gene") %>% 
  column_to_rownames(var = "gene")

all_expr <-all_expr %>% 
  mutate(
    zero_count = apply(., 1, function(x) sum(x == 0))
  ) %>%
  dplyr::filter(zero_count <= (ncol(.)-1) / 2) %>%  # 只保留0的个数不超过一半的基因
  dplyr::select(-zero_count)

combat_data<-ComBat(all_expr,batch=batch)

expr_tcga_article <- as.data.frame(combat_data)

# 批次效应 expr_tcga_article and srcc

expr_tumor <- expr_tumor %>% 
  column_to_rownames(var = "sample") %>% 
  t() %>% 
  as.data.frame()
expr_tumor <- log2(expr_tumor+1)  
expr_tumor <- expr_tumor %>% 
  rownames_to_column(var = "gene")
expr_tcga_article <- expr_tcga_article %>% 
  rownames_to_column(var = "gene")

batch<-c(rep("TCGA_article",1677),rep("tumor",28))

all_expr <- expr_tcga_article %>% 
  inner_join(.,expr_tumor,by="gene") %>% 
  column_to_rownames(var = "gene")

all_expr <-all_expr %>% 
  mutate(
    zero_count = apply(., 1, function(x) sum(x == 0))
  ) %>%
  dplyr::filter(zero_count <= (ncol(.)-1) / 2) %>%  # 只保留0的个数不超过一半的基因
  dplyr::select(-zero_count)

combat_data<-ComBat(all_expr,batch=batch)

expr_tumor <- as.data.frame(combat_data)


# 4. last data ------------------------------------------------------------

expr_tumor_genes <- expr_tumor %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene)

expr_srcc_genes <- expr_srcc %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene)

genes <- expr_srcc_genes %>% 
  inner_join(.,expr_tumor_genes,by="gene")

expr_tumor <- expr_tumor %>% 
  rownames_to_column(var = "gene") %>% 
  filter(gene %in% genes$gene)

expr_srcc <- expr_srcc %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  filter(gene %in% genes$gene)

save(expr_srcc,expr_tumor,
     file = "expr.Rda")

# 5.test ------------------------------------------------------------------

load("expr.Rda")

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

srcc_sample <- expr_srcc %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample")
srcc_sample <- srcc_sample$sample

srcc_test1 <- expr_srcc %>% 
  dplyr::select(TUMC_T83)
srcc_test2 <- expr_srcc%>% 
  dplyr::select(TUMC_T84)
srcc_test3 <- expr_srcc%>% 
  dplyr::select(TUMC_T85)

expr_diff1 <- cbind(expr_tumor[,1:100],srcc_test1)
expr_diff2 <- cbind(expr_tumor[,1:100],srcc_test2)
expr_diff3 <- cbind(expr_tumor[,1:100],srcc_test3)
expr_diff1 <- round(expr_diff1 )
expr_diff2 <- round(expr_diff2)
expr_diff3 <- round(expr_diff3)
expr_diff1[expr_diff1<0] <- 0
expr_diff2[expr_diff2<0] <- 0
expr_diff3[expr_diff3<0] <- 0

expr_diff1_group <- expr_diff1 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::select(sample) %>% 
  mutate(group=ifelse(str_sub(sample,1,8)=="TUMC_T83","SRCC","Tumor")) %>% 
  column_to_rownames(var = "sample")


expr_diff2_group <- expr_diff2 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::select(sample) %>% 
  mutate(group=ifelse(str_sub(sample,1,8)=="TUMC_T84","SRCC","Tumor"))%>% 
  column_to_rownames(var = "sample")

expr_diff3_group <- expr_diff3 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::select(sample) %>% 
  mutate(group=ifelse(str_sub(sample,1,8)=="TUMC_T85","SRCC","Tumor"))%>% 
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
k1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -2)
k2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > 2)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
DEG1 <- DEG %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene,change)
DEG1_genes <- DEG1 %>% 
  filter(change != "NOT")
DEG1_nogenes <- DEG1 %>% 
  filter(change == "NOT")

# DIFF2
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
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -2)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > 2)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
DEG2 <- DEG %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene,change)
DEG2_genes <- DEG2 %>% 
  filter(change != "NOT")
DEG2_nogenes <- DEG2 %>% 
  filter(change == "NOT")

# DIFF3
expr_diff3_group$group <- as.factor(expr_diff3_group$group)
colData <- data.frame(row.names = colnames(expr_diff3),                      
                      condition = expr_diff3_group$group)
dds <- DESeqDataSetFromMatrix(countData = expr_diff3,                              
                              colData = colData,                              
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition",rev(levels(expr_diff3_group$group))))
DEG <- res[order(res$pvalue),] %>% 
  as.data.frame()
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -2)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > 2)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
DEG3 <- DEG %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene,change)
DEG3_genes <- DEG3 %>% 
  filter(change != "NOT")
DEG3_nogenes <- DEG3 %>% 
  filter(change == "NOT")

a=list(DEG1=DEG1_genes$gene,
       DEG2=DEG2_genes$gene,
       DEG3=DEG3_genes$gene)
ggvenn(a,c("DEG1","DEG2","DEG3"),fill_color = c("#FF4500", "#00B2EE","blue"))



# 6.cycle -----------------------------------------------------------------


library(DESeq2)
library(tidyverse)
load("expr.Rda")

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

srcc_sample <- expr_srcc %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  arrange(sample)

srcc_sample <- srcc_sample$sample
srcc_sample <- srcc_sample[1:24]

result <- list()
for (i in srcc_sample){
  srcc_test1 <- expr_srcc %>% 
    dplyr::select(i)
  expr_diff1 <- cbind(expr_tumor,srcc_test1)
  expr_diff1 <- round(expr_diff1 )
  expr_diff1[expr_diff1<0] <- 0
  
  expr_diff1_group <- expr_diff1 %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "sample") %>% 
    dplyr::select(sample) %>% 
    mutate(group=ifelse(str_sub(sample,1,20)==i,"SRCC","Tumor")) %>% 
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
  head(DEG)
  DEG <- DEG %>% 
    rownames_to_column(var = "gene") %>% 
    dplyr::select(gene,change,padj,log2FoldChange)
  DEG_genes <- DEG %>% 
    filter(change != "NOT")
  result[[i]] <- DEG_genes
}

result_degs <- list()
for (i in srcc_sample[1:24]){
  data <- result[[i]]$gene
  result_degs[[i]] <- data
}

save(result,result_degs,
     file = "result_degs_fc1.Rda")


# 7.upset -----------------------------------------------------------------
library(UpSetR)
library(randomcoloR)

cols <- distinctColorPalette(4) #差异明显的60种
load("result_degs_fc1.Rda")

names <- names(result_degs[1:24])
upset(fromList(result_degs[1:24]),
      order.by="freq",#排序方式
      nsets=24,#展示几个集合，按照数量从大到小排列，或者使用sets参数指定集合名字
      mb.ratio=c(0.55,0.45),#条形图和矩阵的相对比例
      number.angles=0,#条形图上面数字角度
      point.size=1.5,#点的大小
      line.size=1.2,#线条粗细
      mainbar.y.label="size of intersection",#上面条形图的标题
      sets.x.label = "the number of each degs",#坐标条形图的标题
      text.scale = c(1, 1, 1, 1, 1, 1),#元素大小
      matrix.color="firebrick",
      main.bar.color="steelblue",
      sets.bar.color="grey70",
      queries=list(
        list(query=intersects,params=list(names),color="orange",active=T)
      )
)

intersect_genes <- result_degs[[1]]
intersect_num <- c()
for(i in 2:length(result_degs)){
  
  intersect_genes <- intersect(
    intersect_genes,
    result_degs[[i]]
  )
  
  intersect_num[i] <- length(intersect_genes)
}
intersect_num

degs_freq <- data.frame(
  intersection = c(2:24),
  genes = c(602,435, 383, 358 ,339, 314 ,302, 294, 286, 275, 273, 270, 269, 268,
            262 ,262 ,255, 254 ,254, 253, 253, 246 ,246
             )
) %>% 
  arrange(desc(intersection))

save(degs_freq,
     file = "degs_freq_fc1.Rda")

# 8.交集基因折线图 and 条形图 -------------------------------------------------------------------
# remotes::install_github("idmn/ggview",force=T)

library(ggview)
library(ggrepel)
library(tidyverse)
library(randomcoloR)
load("degs_freq_fc1.Rda")
# intersection,genes

Color <- randomColor(23)

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

save_ggplot(p, "bar_chart.png")

# 9.提取基因 ------------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(randomcoloR)

load("result_degs_fc1.Rda")
load("degs_freq_fc1.Rda")
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
names <- names(result_degs[1:24])
genes_all <- data.frame(gene=as.character())
for (i in names){
  dat <- data.frame(gene=result_degs[[i]])
  genes_all <- rbind(genes_all,dat)
}
genes_all <- genes_all %>%
  count(gene)

genes_core <- genes_all %>% 
  filter(n==24)

genes_nocore <- genes_all %>% 
  filter(n!=24)

save(genes_all,genes_core,genes_nocore,
     file = "genes_core_nocore_fc1.Rda")

# 10.enrichment -----------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

load("genes_core_nocore_fc1.Rda")
# write.csv(genes_all,file = "genes_all.csv")
# write.csv(genes_core,file = "genes_core.csv")
# write.csv(genes_nocore,file = "genes_nocore.csv")
# GO
go_result <- enrichGO(gene = genes_nocore$gene,
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
  genes_nocore$gene,
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

KEGG_genes <- data.frame(gene=as.character())
KEGG_ID <- KEGG_data$ID
for(i in KEGG_ID){
  p <- getKEGG(i)
  genes <- data.frame(gene=p[[1]])
  KEGG_genes <- rbind(KEGG_genes,genes) %>% 
    distinct(gene,.keep_all = T)
}

genes_exclude <- rbind(KEGG_genes,go_genes) %>% 
  distinct(gene,.keep_all = T)

save(genes_exclude,go_genes,KEGG_genes,go_data,KEGG_data,
     file = "genes_exclude.Rda")
# 11.机器学习筛选(没用到) ---------------------------------------------------------------
library(tidyverse)
library(Mime1)

load("genes_core_nocore.Rda")
load("expr.Rda")

expr_srcc <- expr_srcc %>% 
  column_to_rownames(var = "gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(Var="Y") %>% 
  dplyr::select(ID,Var,everything())
expr_tumor <- expr_tumor %>% 
  column_to_rownames(var = "gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(Var="N") %>% 
  dplyr::select(ID,Var,everything()) #%>% 
#dplyr::filter(str_sub(ID,1,4)=="TUMC")

sample_tumor <- sample(1:nrow(expr_tumor), 187)  # 随机抽取600个样本
training_tumor <- expr_tumor[sample_tumor,]  # 创建新的数据框包含随机抽取的样本
sample_srcc <- sample(1:nrow(expr_srcc), 13)
training_srcc <- expr_srcc[sample_srcc,]
training <- rbind(training_srcc,training_tumor)
validation <- rbind(expr_tumor[1001:1187,],expr_srcc[-sample_srcc,])

list_train_vali_Data <- list(training=training,
                             validation=validation)

genelist <- genes_core$gene

res.ici <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$training,
                                    list_train_vali_Data = list_train_vali_Data,
                                    candidate_genes = c("AARSD1","ABCB6"),
                                    methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                    seed = 5201314,
                                    cores_for_parallel = 60
)

auc_vis_category_all(res.ici,dataset = c("training","validation"),
                     order= c("training","validation"))

plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.ici,model_name = i,dataset = c("training","validation"),
                                   order= c("training","validation"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)

# 12.ROC_genes_core (没用到)-------------------------------------------------------
library(plotROC)
library(ggsci)
library(tidyverse)

load("genes_core_nocore.Rda")
load("expr.Rda")

genes_core <- genes_core$gene
expr_srcc <- expr_srcc %>% 
  column_to_rownames(var = "gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(Var="SRCC") %>% 
  dplyr::select(ID,Var,everything())
expr_tumor <- expr_tumor %>% 
  column_to_rownames(var = "gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(Var="Tumor") %>% 
  dplyr::select(ID,Var,everything()) %>% 
  dplyr::filter(str_sub(ID,1,4)=="TUMC")

roc_data <- rbind(expr_srcc,expr_tumor) %>% 
  column_to_rownames(var = "ID") %>% 
  dplyr::select(Var,genes_core) %>% 
  dplyr::rename(group=Var)

long_roc_data <- melt_roc(roc_data, "group", genes_core)
p <- ggplot(long_roc_data, aes(d = factor(D), m = M)) +
  geom_roc(n.cuts = 0)+
  facet_wrap(~ name,nrow = 2)
roc_auc <- as.data.frame(calc_auc(p)) %>% 
  dplyr::select(name,AUC) %>% 
  dplyr::filter(AUC>0.75|AUC<0.25)



# 13.验证集test1 ------------------------------------------------------------------
library(ComplexHeatmap)
library(ggvenn)
library(DESeq2)
library(tidyverse)
load("genes_core_nocore.Rda")
load("expr.Rda")
load("genes_exclude.Rda")
load("result_degs.Rda")

train_name <- names(result_degs)
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
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -2)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > 2)
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
genes_exclude <- genes_exclude$gene

list <- list(test=genes_test1,exclude=genes_exclude)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_remain1 <- extract_comb(m,"10")

list <- list(test=test_remain1,core=genes_core$gene)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_core1 <- extract_comb(m,"11")

save(test_core1,test_remain1,test1_DEG_genes,genes_test1,
     file = "test1.Rda")


# 14.验证集test2 -------------------------------------------------------------
library(ComplexHeatmap)
library(ggvenn)
library(DESeq2)
library(tidyverse)
load("genes_core_nocore.Rda")
load("expr.Rda")
load("genes_exclude.Rda")
load("result_degs.Rda")

train_name <- names(result_degs)
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
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -2)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > 2)
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
genes_exclude <- genes_exclude$gene

list <- list(test=genes_test2,exclude=genes_exclude)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_remain2 <- extract_comb(m,"10")

list <- list(test=test_remain2,core=genes_core$gene)
m <- make_comb_mat(list)
#获取输入集合的名字
set_name(m)
#获取组合集合的名字，111代表在a，b，c都存在的元素
comb_name(m)
test_core2 <- extract_comb(m,"11")

save(test_core2,test_remain2,test2_DEG_genes,genes_test2,
     file = "test2.Rda")

# save(test_core1,file="genes.Rda")

# 15.The ranks of test1  --------------------------------------------------
load("result_degs.Rda")
load("test1.Rda")
load("genes_core_nocore.Rda")

library(tidyverse)

genes_core <- genes_core$gene

dat_p <- data.frame()
for (i in 1:24){
  dat <- result[[i]] %>% 
    filter(gene %in% genes_core) %>% 
    column_to_rownames(var = "gene") %>% 
    dplyr::select(pvalue) %>% 
    t() %>% 
    as.data.frame()
  dat_p <- rbind(dat_p,dat)
}
genes_p <- apply(dat_p, 2, median, na.rm = TRUE) %>%  # 计算每列的中位数
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  rename(p = ".") %>% 
  arrange(p) %>% 
  mutate(n = 1:nrow(.))  # 自动根据行数生成序号


dat_log2FC <- data.frame()
for (i in 1:24){
  dat <- result[[i]] %>% 
    filter(gene %in% genes_core) %>% 
    column_to_rownames(var = "gene") %>% 
    dplyr::select(log2FoldChange) %>% 
    t() %>% 
    as.data.frame()
  dat_log2FC <- rbind(dat_log2FC,dat)
}
genes_log2FC <- apply(dat_log2FC, 2, median, na.rm = TRUE) %>%  # 计算每列的中位数
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  rename(log2FC = ".") %>% 
  mutate(log2FC = abs(log2FC)) %>%  # 取中位数的绝对值
  arrange(desc(log2FC)) %>%  # 按中位数的绝对值降序排列
  mutate(n = 1:nrow(.))  # 自动生成序号列

test1_p <- data.frame(gene=test_core1) %>% 
  inner_join(genes_p,.,by="gene") %>% 
  dplyr::select(gene,n)

test1_log2FC <- data.frame(gene=test_core1) %>% 
  inner_join(genes_log2FC,.,by="gene") %>% 
  dplyr::select(gene,n)

save(genes_p,genes_log2FC,
     test1_log2FC,test1_p,
     file = "test1_ranks.Rda")

# 16.The eventual result --------------------------------------------------
load("test2.Rda")
load("test1.Rda")
load("genes_core_nocore.Rda")

library(tidyverse)
library(randomcoloR)
library(ggvenn)
library(ggview)

cols <- distinctColorPalette(3)
a <- list(test_remain1=test_remain1,
          test_remain2=test_remain2,
          genes_core=genes_core$gene)
ggvenn(a,c("test_remain1","test_remain2","genes_core"),
       fill_color = cols)

cols <- distinctColorPalette(2)
a <- list(sample_1=test_core1,
          sample_2=test_core2)
p <- ggvenn(a,c("sample_1","sample_2"),
       fill_color = cols,
       show_percentage=F,
       set_name_size=15,
       text_size = 15)+
  ggview::canvas(width = 16,height = 16,dpi = 1000)
p
save_ggplot(p, "venn.png")

# 17.富集分析 -----------------------------------------------------------------
library(readxl)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
load("test2.Rda")
load("test1.Rda")

go_result <- enrichGO(gene = test_core1,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05)
go_data <- as.data.frame(go_result)



diff_entrez<-bitr(
  test_core1,
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

# write_xlsx(go_data, path = "go_data.xlsx")
# write_xlsx(KEGG_data, path = "KEGG_data.xlsx")

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

save_ggplot(p, "enrichment.png")
  

# 18.Rank Product-RankProd ------------------------------------------------

load("C:/Users/21118/Desktop/1.A novel biological function-based method for mining core genes in rare disease with limited cases/result_degs.Rda")
view(result$TUMC_T103)


# 19.GSE281917 ------------------------------------------------------------
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

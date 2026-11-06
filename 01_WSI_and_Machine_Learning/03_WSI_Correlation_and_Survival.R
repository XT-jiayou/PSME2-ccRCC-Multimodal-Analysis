rm(list=ls())
## 读取TCGA-KIRC病理图像特征、UPS相关基因表达和基因集评分数据
setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/')
load('./features_504.Rdata')
load('./gene.Rdata')
load('./gsva.Rdata')
setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/55/model/')
load('./final_feature.Rdata')
## 以下为单因素Cox筛选特征的备用代码（当前未执行）
#final_feature2=read.csv('./unicox.csv')
#final_feature=final_feature2$feature


features=features[,final_feature]


identical(rownames(features),rownames(gene))
identical(rownames(features),rownames(gsva))


# 加载相关性分析及热图绘制所需R包
library(corrplot)  # 绘制Spearman相关性热图
library(Hmisc)     # 计算相关系数及对应P值

# 根据共同患者编号整合病理特征和UPS相关基因表达数据
#df<-cbind(features,gene)
# 1. 获取病理特征、基因表达和基因集评分数据的共同患者编号
common_samples <- Reduce(
  intersect,
  list(rownames(features), rownames(gene), rownames(gsva))
)

# 2. 按共同患者编号提取数据并统一样本顺序
features <- features[common_samples, , drop = FALSE]
gene <- gene[common_samples, , drop = FALSE]
gsva <- gsva[common_samples, , drop = FALSE]

# 3. 检查各数据矩阵的维度和样本顺序
dim(features)
dim(gene)
dim(gsva)

identical(rownames(features), rownames(gene))
identical(rownames(features), rownames(gsva))

# 4. 合并病理图像特征与UPS相关基因表达矩阵
df <- cbind(features, gene)

# 5. 计算合并矩阵的Spearman相关系数和P值
cor_results <- Hmisc::rcorr(as.matrix(df), type = "spearman")

# 备用相关性计算语句（当前未执行）
#cor_results <- rcorr(as.matrix(df), type = "spearman")


# 提取Spearman相关系数矩阵和P值矩阵
cor_df <- cor_results$r
pval_df <- cor_results$P


cor_df=cor_df[colnames(features),colnames(gene)]
# 为改善热图可读性，将特征名和基因名截取为前10个字符
colnames(cor_df) <- stringr::str_sub(colnames(cor_df), 1, 10)
rownames(cor_df) <- stringr::str_sub(rownames(cor_df), 1, 10)

pval_df=pval_df[colnames(features),colnames(gene)]
# 对P值矩阵采用相同的名称截取方式，保持热图标签一致
colnames(pval_df) <- stringr::str_sub(colnames(pval_df), 1, 10)
rownames(pval_df) <- stringr::str_sub(rownames(pval_df), 1, 10)


mycol<- colorRampPalette(c("#06a7cd", "white", "#e74a32"), alpha = TRUE)

# 绘制病理图像特征与UPS相关基因表达的Spearman相关性热图
corrplot(cor_df, p.mat = pval_df,
         method = "pie", col = mycol(100),
         insig = "blank",  # 不显示未达到显著性阈值的相关性
         cl.lim = c(-1, 1),  # 设置相关系数颜色映射范围
         title = "Correlation Heatmap (Spearman)",  # 热图标题
         addCoef.col = "black",  # 显示相关系数数值
         number.cex = 0.4, cl.cex = 0.8, # 设置相关系数与图例字体大小
         tl.cex = 1,   # 设置行列标签字体大小
)


dim(features)
colnames(features)[1:10]
ls(pattern = "features")
library(dplyr)
library(reshape2)
## 构建病理图像特征与UPS相关基因之间的相关性网络
# 重新计算用于网络分析的Spearman相关系数和P值
cor_results <- rcorr(as.matrix(df), type = "spearman")


# 提取相关系数矩阵和P值矩阵
cor_df <- cor_results$r
pval_df <- cor_results$P


cor_df=cor_df[colnames(features_net),colnames(gene)]
pval_df=pval_df[colnames(features_net),colnames(gene)]


cor_long <- melt(cor_df)
pval_long <- melt(pval_df)

# 整合每组病理特征与基因对应的相关系数和P值
long_df <- merge(cor_long, pval_long, by = c("Var1", "Var2"))
colnames(long_df) <- c("Feature", "Gene", "Correlation", "P_Value")

# 保留名义P值不大于0.05的特征-基因关联
significant_df <- long_df %>% filter(P_Value <= 0.05)

# 检查筛选结果
head(significant_df)
library(stringr)
library(tidygraph)

### 构建相关性网络的节点与边数据
net=significant_df[,c(1,2,3)]
colnames(net)=c('to','from','weight')
# 按ResNet和CellProfiler来源对病理特征节点分类
group1 <- net
group1$Group <- str_sub(group1$to, 1, 6)
group1$Group <- ifelse(group1$Group == 'resnet', 'resnet', 'cellprofiler')

group1 <- group1[, c(1, 4)]
group1= group1[!duplicated(group1$to),]
colnames(group1)<- c('ID', 'Group')

group2 <- net
group2$Group <- 'Gene'
group2 <- group2[, c(2, 4)]
group2= group2[!duplicated(group2$from),]
colnames(group2)<- c('ID', 'Group')
group=rbind(group1,group2)

net$weight <- abs(net$weight)

# 2. 根据特征-基因关联构建无向网络对象
graph <- tbl_graph(edges = net, nodes = group, directed = FALSE)

dev.off()
library(ggraph)
# 3. 绘制ResNet、CellProfiler特征与UPS相关基因的关联网络
ggraph(graph, layout = "fr") + 
  geom_edge_link(aes(width = abs(weight)), alpha = 0.6) +
  geom_node_point(aes(color = Group, size = centrality_degree())) +
  geom_node_text(aes(label =ID), repel = TRUE, size = 3) +  # 标注网络节点名称
  scale_edge_width(range = c(0.2, 3)) +
  scale_size(range = c(3, 10)) +
  theme_void() +
  labs(title = "Network Graph") +
  scale_color_manual(values = c("resnet" = "firebrick", 
                                "cellprofiler" = "steelblue", 
                                "Gene" = "grey")


library(dplyr)
library(ggplot2)

load("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/55/model/riskscore.Rdata")

ls()
str(features_train2, max.level = 1)
str(features_vali, max.level = 1)

head(rownames(features_train2))
head(rownames(features_vali))

colnames(features_train2)
colnames(features_vali)
# 1. 设置训练集的预后风险分组顺序
# train_data应包含PSME2表达量和预后风险分组信息
train_data$RiskGroup <- factor(
  train_data$RiskGroup,
  levels = c("Low", "High")
)

# 2. 设置验证集的预后风险分组顺序
# vali_data应包含PSME2表达量和预后风险分组信息
vali_data$RiskGroup <- factor(
  vali_data$RiskGroup,
  levels = c("Low", "High")
)

# 3. 使用Wilcoxon检验比较训练集不同风险组的PSME2表达
wilcox.test(PSME2 ~ RiskGroup, data = train_data)

# 4. 使用Wilcoxon检验比较验证集不同风险组的PSME2表达
wilcox.test(PSME2 ~ RiskGroup, data = vali_data)

# 5. 计算训练集及验证集各风险组的PSME2表达中位数
aggregate(
  PSME2 ~ RiskGroup,
  data = train_data,
  median
)

aggregate(
  PSME2 ~ RiskGroup,
  data = vali_data,
  median
)

setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/55/')
## TCGA-KIRC队列中PSME2表达与总体生存的关联分析

## 读取UPS相关基因表达矩阵和TCGA-KIRC生存数据
#setwd('~/pathology2/')
load('./gene.Rdata')
library(data.table)
suv=fread('./TCGA-KIRC.survival.tsv.gz',check.names = F,data.table = F)

colnames(suv)[4]='Sample'

suv=suv[!duplicated(suv$Sample),]
rownames(suv)=suv$Sample

## 匹配基因表达矩阵与临床生存数据的共同患者
sss=intersect(rownames(gene),rownames(suv))
suv=suv[sss,]
suv=suv[,c(2,3)]
gene=gene[sss,]

## 合并生存信息与基因表达矩阵
rt=cbind(suv,gene)
library(survival)
library(survminer)
i='PSME2'  ## 指定本研究关注的PSME2基因
res.cut=surv_cutpoint(rt, time="OS.time", event="OS", variables=i)
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(rt[,i]<= cutoff, "Low", "High")
Type=ifelse(rt[,i]<= median(rt[,i]), "Low", "High")

data=rt
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(OS.time, OS) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
## 根据PSME2高低表达分组拟合Kaplan–Meier生存曲线
fit <- survfit(Surv(OS.time, OS) ~ group, data = data)
set.seed(1234)

color=RColorBrewer::brewer.pal(8,'RdPu')
color
ggsurvplot(fit, 
           data=data,
           pval.size=4,
           palette = color[c(3,8)],
           pval = T,
           legend = c(0.8, 0.2),
           font.legend=8,
           xlab="Time(Months)",
           surv.median.line = "hv",
           risk.table=T,
           cumevents=F,
           risk.table.height=.25,ggtheme = theme_minimal()+
             theme(panel.border = element_rect(color = "black", fill = NA, size = 1)  # 为生存曲线添加黑色边框
             ))# 使用简洁主题显示Kaplan–Meier曲线

rm(list=ls())
setwd("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/")

total=read.csv('Total_MyExpt_Image.csv',header = T,check.names = F)
object=read.csv('Object_MyExpt_Image.csv',header = T,check.names = F)


gran=grep('^Granularity',colnames(total),value = T)
inten=grep('^Intensity_',colnames(total),value = T)
texture=grep('^Texture_',colnames(total),value = T)
gran
inten
texture

total2=total[,c('Metadata_FileLocation',gran,inten,texture)]
max(total2[,2:ncol(total2)])

colnames(total2)[2:ncol(total2)]=paste0('Total_',colnames(total2)[2:ncol(total2)])

colnames(object)
texture2=grep('^Texture_',colnames(object),value = T)
Mean=grep('^Mean_',colnames(object),value = T)
Median=grep('^Median_',colnames(object),value = T)
StDev=grep('^StDev_',colnames(object),value = T)


object2=object[,c('Metadata_FileLocation',Mean,Median,StDev,texture2)]
colnames(object2)[2:ncol(object2)]=paste0('Object_',colnames(object2)[2:ncol(object2)])


library(dplyr)
cp=merge(total2,object2,by='Metadata_FileLocation')


library(stringr)
# 从图像文件路径中提取TCGA患者编号及图像块编号
cp <- cp %>%
  mutate(
    Sample = str_extract(Metadata_FileLocation, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"),
    tile = str_extract(Metadata_FileLocation, "tile_\\d+")
  )
head(cp$Sample)

# 检查提取的患者编号、图像块编号及特征列
head(cp)

cp$Metadata_FileLocation=NULL

cp=select(cp,'Sample','tile',everything())

cp$tile=NULL
library(dplyr)

# 对同一患者的CellProfiler图像块特征计算均值，生成患者级特征
cp2 <- cp %>%
  group_by(Sample) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%  # 对数值型特征计算均值并忽略缺失值
  ungroup()  # 取消分组，保留患者级特征表

cp2=as.data.frame(cp2)
cp2=cp2[!duplicated(cp2$Sample),]
rownames(cp2)=cp2$Sample
cp2$Sample=NULL



resnet=read.csv('./resnet50_features.csv',header = T,check.names = F)

colnames(resnet)[1]="Metadata_FileLocation"

library(stringr)
# 从图像文件路径中提取TCGA患者编号及图像块编号
resnet<-resnet %>%
  mutate(
    Sample = str_extract(Metadata_FileLocation, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"),
    tile = str_extract(Metadata_FileLocation, "tile_\\d+")
  )


# 检查提取的患者编号、图像块编号及特征列
head(resnet)

resnet$Metadata_FileLocation=NULL

resnet=select(resnet,'Sample','tile',everything())

# 对同一患者的ResNet-50图像块特征计算均值
resnet2 <- resnet %>%
  group_by(Sample) %>%  # 按TCGA患者编号分组
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%  # 对数值型特征计算均值并忽略缺失值
  ungroup()  # 取消分组，保留患者级特征表

resnet2=as.data.frame(resnet2)
resnet2=resnet2[!duplicated(resnet2$Sample),]

rownames(resnet2)=resnet2$Sample

resnet2$Sample=NULL

colnames(resnet2)=paste0('resnet',colnames(resnet2))

## 整合CellProfiler与ResNet-50提取的患者级病理图像特征

# 识别两类病理图像特征共有的TCGA患者
common_samples <- intersect(rownames(cp2), rownames(resnet2))

# 检查可匹配的患者数量
length(common_samples)

# 根据相同患者编号对齐两类特征矩阵
cp2_new <- cp2[common_samples, , drop = FALSE]
resnet2_new <- resnet2[common_samples, , drop = FALSE]

# 合并CellProfiler形态学特征与ResNet-50深度学习特征
features <- cbind(cp2_new, resnet2_new)

# 检查合并后特征矩阵的维度
dim(features)

# 保存整合后的患者级病理特征矩阵
save(features, file = "features_504.Rdata")
write.csv(features, file = "features_504.csv", quote = FALSE)
## 整合CellProfiler与ResNet-50提取的患者级病理图像特征
#features=cbind(cp2,resnet2)
#save(features,file ='features.Rdata')
#write.csv(features,file ='features.csv',quote = F)

## TCGA-KIRC转录组数据整理与基因注释
# 加载转录组数据处理所需的R包
#setwd('~/pathology2/')
#例如setwd('e:/work/aggre')

#install.packages('data.table')
library(data.table)
## 读取TCGA-KIRC转录组TPM表达矩阵
expr= fread(file='TCGA-KIRC.star_tpm.tsv.gz',header = T,sep = '\t')
#https://xenabrowser.net/datapages/

library(tibble)                 
expr=column_to_rownames(expr,var = 'Ensembl_ID')
max(expr)


library(tibble)
exprSet = rownames_to_column(expr, var ='gene_id')
library(tidyr)
exprSet <- exprSet %>% 
  tidyr::separate(gene_id,into = c("gene_id"),sep="\\.")

## 读取用于Ensembl基因ID与基因符号映射的注释文件
load('gtf_df.Rdata')


### 提取蛋白编码基因的表达数据
mRNA_exprSet <- gtf_df %>% 
  dplyr::filter(type=="gene",gene_biotype=="protein_coding") %>% 
  dplyr::select(c(gene_name,gene_id,gene_biotype)) %>% 
  dplyr::inner_join(exprSet,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_biotype,sep = " | ")

exprSet_mRNA <-mRNA_exprSet %>%
  tidyr::separate(gene_id,c('gene_name','gene_id','gene_biotype'),
                  sep = " \\| ")
exprSet_mRNA<- exprSet_mRNA[,-c(2,3)]
exprSet_mRNA<-exprSet_mRNA[exprSet_mRNA$gene_name != 'NA',]
rownames(exprSet_mRNA)=NULL
# 对映射至同一基因符号的表达值取平均
exprSet_tcga_mRNA<-aggregate(x=exprSet_mRNA[,2:(ncol(exprSet_mRNA))],by=list(exprSet_mRNA$gene_name),FUN=mean)
exprSet_tcga_mRNA <- column_to_rownames(exprSet_tcga_mRNA, var='Group.1')
max(exprSet_tcga_mRNA)
# 完成TCGA-KIRC蛋白编码基因表达矩阵整理
exprSet=exprSet_tcga_mRNA
save(exprSet,file = 'KIRC_tpm.Rdata')
#load('KIRC_tpm.Rdata')

### 计算TCGA-KIRC肿瘤样本的UPS相关基因集评分
load('KIRC_tpm.Rdata')

library(stringr)
k1 = str_starts(colnames(exprSet),"TCGA")
k2 = as.numeric(str_sub(colnames(exprSet),14,15))<10
## 运行时可能出现非数值转换产生NA的警告
## 检查TCGA样本编号对应的样本类型
table(k1&k2)


group_list = ifelse(k1&k2,"tumor","normal")
group_list = factor(group_list,levels = c("normal","tumor"))

## 筛选TCGA-KIRC肿瘤样本用于后续分析
exprSet=exprSet[,group_list=='tumor']

#geneset=read.table('./geneset.txt')
geneset <- read.csv("UPS_KEGG_genes.csv", header = TRUE)
# 提取UPS相关基因的基因符号
geneset <- geneset$Gene

# 将UPS基因集合整理为GSVA所需的列表格式
geneset <- list(geneset)

# 将基因集命名为geneset
names(geneset) <- "geneset"
## 使用与该脚本调用方式兼容的GSVA版本进行ssGSEA评分
library(GSVA)
gsva_matrix<- gsva(as.matrix(exprSet), geneset,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

gsva_matrix=as.data.frame(t(gsva_matrix))


gsva=rownames_to_column(gsva_matrix,var='Sample')
gsva$Sample=str_sub(gsva$Sample,1,12)
## 将同一患者的多个肿瘤样本评分汇总为患者级均值
gsva<-aggregate(x=gsva[,2:(ncol(gsva))],by=list(gsva$Sample),FUN=mean)


rownames(gsva)=gsva$Group.1
gsva$Group.1=NULL
colnames(gsva)='geneset'
save(gsva,file ='gsva.Rdata')

### 提取UPS相关基因在TCGA-KIRC肿瘤样本中的表达矩阵
load('KIRC_tpm.Rdata')

library(stringr)
k1 = str_starts(colnames(exprSet),"TCGA")
k2 = as.numeric(str_sub(colnames(exprSet),14,15))<10
## 运行时可能出现非数值转换产生NA的警告
## 检查TCGA样本编号对应的样本类型
table(k1&k2)


group_list = ifelse(k1&k2,"tumor","normal")
group_list = factor(group_list,levels = c("normal","tumor"))

## 筛选TCGA-KIRC肿瘤样本用于后续分析
exprSet=exprSet[,group_list=='tumor']

#geneset=read.table('./geneset.txt')
#geneset=geneset$V1
geneset <- read.csv("UPS_KEGG_genes.csv", header = TRUE)
# 提取UPS相关基因的基因符号
geneset <- geneset$Gene
gene=exprSet[rownames(exprSet) %in% geneset,]

gene=as.data.frame(t(gene))

gene=rownames_to_column(gene,var='Sample')
gene$Sample=str_sub(gene$Sample,1,12)
gene<-aggregate(x=gene[,2:(ncol(gene))],by=list(gene$Sample),FUN=mean)

rownames(gene)=gene$Group.1
gene$Group.1=NULL

save(gene,file ='gene.Rdata')


### 病理特征与UPS基因集评分的关联分析
## 读取并整合病理图像特征与基因集评分
#setwd('~/pathology2/')
#load("./features.Rdata")
#load('./gsva.Rdata')
setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/55/')
### 匹配病理图像特征与UPS基因集评分的共有患者
ss=intersect(rownames(features),rownames(gsva))


# 固定随机种子以便复现样本划分或交叉验证结果
set.seed(123456)

# 按当前代码中的比例随机抽取训练样本（此处实际为50%，需与最终Methods核对）
train_set <- sample(ss, size = floor(length(ss) * 0.5))

# 将其余样本划分为内部验证集
validation_set <- setdiff(ss, train_set)


features_train=features[train_set,]
gsva_train=gsva[train_set,,drop=F]

features_vali=features[validation_set,]

save(features_vali,file='features_vali.Rdata')


### 在训练集中筛选与UPS基因集评分相关的病理图像特征
# 初始化Spearman相关系数及P值结果表
result <- data.frame(
  Feature = colnames(features_train),  # 病理图像特征名称
  Correlation = numeric(ncol(features_train)),
  P_Value = numeric(ncol(features_train))
)

# 逐一计算训练集病理特征与UPS基因集评分的Spearman相关性
for (i in 1:ncol(features_train)) {
  test <- cor.test(features_train[, i], as.vector(gsva_train$geneset), method = "spearman")
  result$Correlation[i] <- test$estimate  # Spearman相关系数
  result$P_Value[i] <- test$p.value  # 双侧相关性检验P值
}


result= na.omit(result)

result$FDR <- p.adjust(result$P_Value, method = "fdr")
result$Bonferroni <- p.adjust(result$P_Value, method = "bonferroni")


# 检查相关性分析结果
head(result)

# 导出病理特征与UPS基因集评分的相关性结果
write.csv(result, "correlation_results.csv", row.names = FALSE)

## 按当前脚本的名义P值阈值筛选相关特征（另计算了FDR供核对）
result2=result[result$P_Value<0.05,]
#result2=result[result$FDR<0.05,]


### 构建经过相关性筛选的训练集特征矩阵
features_train2=features_train[,result2$Feature]
save(features_train2,result2,file ='features_train_cor.Rdata')


library(tidygraph)
library(ggraph)
library(dplyr)
library(stringr)


## 绘制UPS基因集评分与病理图像特征的相关性网络
# 整理网络边表与节点信息
net <- result2[,1:2]
net$Score <- 'Geneset'
colnames(net) <- c('to', 'weight', 'from')

# 将相关系数转换为数值类型
net$weight <- as.numeric(net$weight)

# 根据特征来源区分ResNet-50、CellProfiler和UPS基因集节点
group <- net
group$Group <- str_sub(group$to, 1, 6)
group$Group <- ifelse(group$Group == 'resnet', 'resnet', 'cellprofiler')
group <- group[, c(1, 4)]
colnames(group) <- c('ID', 'Group')
group <- rbind(group, data.frame(ID = 'Geneset', Group = 'Geneset'))
net$weight <- abs(net$weight)

# 根据节点与边表构建无向网络
graph <- tbl_graph(edges = net, nodes = group, directed = FALSE)

dev.off()

# 绘制特征相关性网络并设置节点和边的显示方式
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
                                "Geneset" = "grey"))





### 使用Elastic Net进一步筛选UPS相关病理图像特征
library(glmnet)

y <- as.matrix(gsva_train)

# 将训练集病理图像特征转换为模型输入矩阵
X <- as.matrix(features_train2)

# 设置Elastic Net混合参数alpha
alpha_value <- 0.2  # 当前分析使用alpha = 0.2

# 通过交叉验证选择正则化参数lambda
set.seed(123)  # 固定随机种子以便复现样本划分或交叉验证结果
cv_model <- cv.glmnet(X, y, alpha = alpha_value)

# 使用最优lambda拟合Elastic Net模型
best_lambda <- cv_model$lambda.min
final_model <- glmnet(X, y, alpha = alpha_value, lambda = best_lambda)

# 提取模型回归系数
coef_values <- coef(final_model)
coef_values <- as.data.frame(as.matrix(coef_values))

# 提取非零系数对应的候选病理特征
selected_features <- coef_values[coef_values$s0 != 0, , drop = FALSE]

# 检查Elastic Net筛选得到的特征
print(selected_features)

# 提取筛选后特征的名称
selected_feature_names <- rownames(selected_features)
print(selected_feature_names)
selected_feature_names=selected_feature_names[2:length(selected_feature_names)]

result3=result2[result2$Feature %in% selected_feature_names,]


## 绘制UPS基因集评分与病理图像特征的相关性网络
# 整理网络边表与节点信息
net <- result3[,1:2]
net$Score <- 'Geneset'
colnames(net) <- c('to', 'weight', 'from')

# 将相关系数转换为数值类型
net$weight <- as.numeric(net$weight)

# 根据特征来源区分ResNet-50、CellProfiler和UPS基因集节点
group <- net
group$Group <- str_sub(group$to, 1, 6)
group$Group <- ifelse(group$Group == 'resnet', 'resnet', 'cellprofiler')
group <- group[, c(1, 4)]
colnames(group) <- c('ID', 'Group')
group <- rbind(group, data.frame(ID = 'Geneset', Group = 'Geneset'))
net$weight <- abs(net$weight)

# 根据节点与边表构建无向网络
graph <- tbl_graph(edges = net, nodes = group, directed = FALSE)

dev.off()

# 绘制特征相关性网络并设置节点和边的显示方式
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
                                "Geneset" = "green"))



### 保存Elastic Net筛选后的训练集病理图像特征矩阵
features_train3=features_train2[,result3$Feature]
#features_train3=features_train2
#result3=result2
save(features_train3,result3,file ='features_train_enet.Rdata')

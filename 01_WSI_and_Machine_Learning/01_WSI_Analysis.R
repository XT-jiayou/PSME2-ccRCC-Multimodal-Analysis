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
# 提取字符串
cp <- cp %>%
  mutate(
    Sample = str_extract(Metadata_FileLocation, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"),
    tile = str_extract(Metadata_FileLocation, "tile_\\d+")
  )
head(cp$Sample)

# 查看新的列
head(cp)

cp$Metadata_FileLocation=NULL

cp=select(cp,'Sample','tile',everything())

cp$tile=NULL
library(dplyr)

# 计算每个样本下各个tile的平均值
cp2 <- cp %>%
  group_by(Sample) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%  # 使用 lambda 函数明确传递 mean 函数
  ungroup()  # 解除分组

cp2=as.data.frame(cp2)
cp2=cp2[!duplicated(cp2$Sample),]
rownames(cp2)=cp2$Sample
cp2$Sample=NULL



resnet=read.csv('./resnet50_features.csv',header = T,check.names = F)

colnames(resnet)[1]="Metadata_FileLocation"

library(stringr)
# 提取字符串
resnet<-resnet %>%
  mutate(
    Sample = str_extract(Metadata_FileLocation, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"),
    tile = str_extract(Metadata_FileLocation, "tile_\\d+")
  )


# 查看新的列
head(resnet)

resnet$Metadata_FileLocation=NULL

resnet=select(resnet,'Sample','tile',everything())

# 计算每个样本下各个瓷砖的平均值
resnet2 <- resnet %>%
  group_by(Sample) %>%  # 按 Sample 分组
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%  # 使用 lambda 函数明确传递 mean 函数
  ungroup()  # 解除分组

resnet2=as.data.frame(resnet2)
resnet2=resnet2[!duplicated(resnet2$Sample),]

rownames(resnet2)=resnet2$Sample

resnet2$Sample=NULL

colnames(resnet2)=paste0('resnet',colnames(resnet2))

## 合并特征

# 1. 找到两个数据框共有的患者编号
common_samples <- intersect(rownames(cp2), rownames(resnet2))

# 2. 查看共有多少个患者
length(common_samples)

# 3. 按照相同患者编号提取数据
cp2_new <- cp2[common_samples, , drop = FALSE]
resnet2_new <- resnet2[common_samples, , drop = FALSE]

# 4. 合并CellProfiler和ResNet特征
features <- cbind(cp2_new, resnet2_new)

# 5. 检查结果
dim(features)

# 6. 保存数据
save(features, file = "features_504.Rdata")
write.csv(features, file = "features_504.csv", quote = FALSE)
## 合并特征
#features=cbind(cp2,resnet2)
#save(features,file ='features.Rdata')
#write.csv(features,file ='features.csv',quote = F)

## 整理测序数据-----------------
# 先装包
#setwd('~/pathology2/')
#例如setwd('e:/work/aggre')

#install.packages('data.table')
library(data.table)
##载入LIHC表达谱文件（无法获取的请联系管理员获取）
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

##载入GRh38.104文件（无法获取的请联系管理员获取），104和103相去甚远
load('gtf_df.Rdata')


###提取mRNA
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
#相同探针取平均值
exprSet_tcga_mRNA<-aggregate(x=exprSet_mRNA[,2:(ncol(exprSet_mRNA))],by=list(exprSet_mRNA$gene_name),FUN=mean)
exprSet_tcga_mRNA <- column_to_rownames(exprSet_tcga_mRNA, var='Group.1')
max(exprSet_tcga_mRNA)
#TCGA最新版本mRNA注释完毕
exprSet=exprSet_tcga_mRNA
save(exprSet,file = 'KIRC_tpm.Rdata')
#load('KIRC_tpm.Rdata')

### 基因集评分-----
load('KIRC_tpm.Rdata')

library(stringr)
k1 = str_starts(colnames(exprSet),"TCGA")
k2 = as.numeric(str_sub(colnames(exprSet),14,15))<10
## Warning: NAs introduced by
## coercion
table(k1&k2)


group_list = ifelse(k1&k2,"tumor","normal")
group_list = factor(group_list,levels = c("normal","tumor"))

## 只保留肿瘤样本
exprSet=exprSet[,group_list=='tumor']

#geneset=read.table('./geneset.txt')
geneset <- read.csv("UPS_KEGG_genes.csv", header = TRUE)
# 提取Gene列
geneset <- geneset$Gene

# 转换成list
geneset <- list(geneset)

# 命名基因集
names(geneset) <- "geneset"
## GSVA：得装旧版!!!!不然疯狂报错
library(GSVA)
gsva_matrix<- gsva(as.matrix(exprSet), geneset,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

gsva_matrix=as.data.frame(t(gsva_matrix))


gsva=rownames_to_column(gsva_matrix,var='Sample')
gsva$Sample=str_sub(gsva$Sample,1,12)
## 取平均
gsva<-aggregate(x=gsva[,2:(ncol(gsva))],by=list(gsva$Sample),FUN=mean)


rownames(gsva)=gsva$Group.1
gsva$Group.1=NULL
colnames(gsva)='geneset'
save(gsva,file ='gsva.Rdata')

### gene
load('KIRC_tpm.Rdata')

library(stringr)
k1 = str_starts(colnames(exprSet),"TCGA")
k2 = as.numeric(str_sub(colnames(exprSet),14,15))<10
## Warning: NAs introduced by
## coercion
table(k1&k2)


group_list = ifelse(k1&k2,"tumor","normal")
group_list = factor(group_list,levels = c("normal","tumor"))

## 只保留肿瘤样本
exprSet=exprSet[,group_list=='tumor']

#geneset=read.table('./geneset.txt')
#geneset=geneset$V1
geneset <- read.csv("UPS_KEGG_genes.csv", header = TRUE)
# 提取Gene列
geneset <- geneset$Gene
gene=exprSet[rownames(exprSet) %in% geneset,]

gene=as.data.frame(t(gene))

gene=rownames_to_column(gene,var='Sample')
gene$Sample=str_sub(gene$Sample,1,12)
gene<-aggregate(x=gene[,2:(ncol(gene))],by=list(gene$Sample),FUN=mean)

rownames(gene)=gene$Group.1
gene$Group.1=NULL

save(gene,file ='gene.Rdata')


###开始！
## feature
#setwd('~/pathology2/')
#load("./features.Rdata")
#load('./gsva.Rdata')
setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/55/')
### 病理特征和基因集评分
ss=intersect(rownames(features),rownames(gsva))


# 设置随机种子以确保可重复性
set.seed(123456)

# 随机取 70% 样本作为训练集
train_set <- sample(ss, size = floor(length(ss) * 0.5))

# 剩余的样本作为验证集
validation_set <- setdiff(ss, train_set)


features_train=features[train_set,]
gsva_train=gsva[train_set,,drop=F]

features_vali=features[validation_set,]

save(features_vali,file='features_vali.Rdata')


### 看看哪些病理特征与基因集评分有关
# 初始化保存结果的空数据框
result <- data.frame(
  Feature = colnames(features_train),  # 特征的列名
  Correlation = numeric(ncol(features_train)),
  P_Value = numeric(ncol(features_train))
)

# 使用 apply 计算 Spearman 相关性和 P 值（病理特征和基因集）
for (i in 1:ncol(features_train)) {
  test <- cor.test(features_train[, i], as.vector(gsva_train$geneset), method = "spearman")
  result$Correlation[i] <- test$estimate  # 相关系数
  result$P_Value[i] <- test$p.value  # P 值
}


result= na.omit(result)

result$FDR <- p.adjust(result$P_Value, method = "fdr")
result$Bonferroni <- p.adjust(result$P_Value, method = "bonferroni")


# 查看结果
head(result)

# 保存结果为 CSV 文件
write.csv(result, "correlation_results.csv", row.names = FALSE)

## 显著的
result2=result[result$P_Value<0.05,]
#result2=result[result$FDR<0.05,]


### 训练集的特征矩阵
features_train2=features_train[,result2$Feature]
save(features_train2,result2,file ='features_train_cor.Rdata')


library(tidygraph)
library(ggraph)
library(dplyr)
library(stringr)


## 可画网络图，内存大
# 1. 数据预处理
net <- result2[,1:2]
net$Score <- 'Geneset'
colnames(net) <- c('to', 'weight', 'from')

# 确保 weight 为数值型
net$weight <- as.numeric(net$weight)

# 节点分组
group <- net
group$Group <- str_sub(group$to, 1, 6)
group$Group <- ifelse(group$Group == 'resnet', 'resnet', 'cellprofiler')
group <- group[, c(1, 4)]
colnames(group) <- c('ID', 'Group')
group <- rbind(group, data.frame(ID = 'Geneset', Group = 'Geneset'))
net$weight <- abs(net$weight)

# 2. 构建图对象
graph <- tbl_graph(edges = net, nodes = group, directed = FALSE)

dev.off()

# 3. 绘制网络图
ggraph(graph, layout = "fr") + 
  geom_edge_link(aes(width = abs(weight)), alpha = 0.6) +
  geom_node_point(aes(color = Group, size = centrality_degree())) +
  geom_node_text(aes(label =ID), repel = TRUE, size = 3) +  # 添加节点文字
  scale_edge_width(range = c(0.2, 3)) +
  scale_size(range = c(3, 10)) +
  theme_void() +
  labs(title = "Network Graph") +
  scale_color_manual(values = c("resnet" = "firebrick", 
                                "cellprofiler" = "steelblue", 
                                "Geneset" = "grey"))





### 可以选择进一步采用弹性网络筛选
library(glmnet)

y <- as.matrix(gsva_train)

# 将 features 转换为矩阵
X <- as.matrix(features_train2)

# 设置弹性网的 alpha 参数
alpha_value <- 0.2  # 可以根据需要调整

# 执行弹性网回归
set.seed(123)  # 设置随机种子以确保可重复性
cv_model <- cv.glmnet(X, y, alpha = alpha_value)

# 根据最佳的 lambda 训练模型
best_lambda <- cv_model$lambda.min
final_model <- glmnet(X, y, alpha = alpha_value, lambda = best_lambda)

# 查看模型系数
coef_values <- coef(final_model)
coef_values <- as.data.frame(as.matrix(coef_values))

# 筛选出非零系数的特征
selected_features <- coef_values[coef_values$s0 != 0, , drop = FALSE]

# 查看选择的特征
print(selected_features)

# 选择的特征名称
selected_feature_names <- rownames(selected_features)
print(selected_feature_names)
selected_feature_names=selected_feature_names[2:length(selected_feature_names)]

result3=result2[result2$Feature %in% selected_feature_names,]


## 可画网络图，内存大
# 1. 数据预处理
net <- result3[,1:2]
net$Score <- 'Geneset'
colnames(net) <- c('to', 'weight', 'from')

# 确保 weight 为数值型
net$weight <- as.numeric(net$weight)

# 节点分组
group <- net
group$Group <- str_sub(group$to, 1, 6)
group$Group <- ifelse(group$Group == 'resnet', 'resnet', 'cellprofiler')
group <- group[, c(1, 4)]
colnames(group) <- c('ID', 'Group')
group <- rbind(group, data.frame(ID = 'Geneset', Group = 'Geneset'))
net$weight <- abs(net$weight)

# 2. 构建图对象
graph <- tbl_graph(edges = net, nodes = group, directed = FALSE)

dev.off()

# 3. 绘制网络图
ggraph(graph, layout = "fr") + 
  geom_edge_link(aes(width = abs(weight)), alpha = 0.6) +
  geom_node_point(aes(color = Group, size = centrality_degree())) +
  geom_node_text(aes(label =ID), repel = TRUE, size = 3) +  # 添加节点文字
  scale_edge_width(range = c(0.2, 3)) +
  scale_size(range = c(3, 10)) +
  theme_void() +
  labs(title = "Network Graph") +
  scale_color_manual(values = c("resnet" = "firebrick", 
                                "cellprofiler" = "steelblue", 
                                "Geneset" = "green"))



### 进一步弹性网络训练集的特征矩阵
features_train3=features_train2[,result3$Feature]
#features_train3=features_train2
#result3=result2
save(features_train3,result3,file ='features_train_enet.Rdata')

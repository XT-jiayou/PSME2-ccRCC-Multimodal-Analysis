rm(list=ls())
# 初始化环境并设置预后模型分析工作目录
setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/55/model/')
#### 预后模型分析依赖包安装 ----
## 检查并安装所需Bioconductor依赖包
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("CoxBoost", quietly = TRUE)) devtools::install_github("binderh/CoxBoost")

if (!requireNamespace("fastAdaboost", quietly = TRUE)) devtools::install_github("souravc83/fastAdaboost")
# 可选依赖包安装命令（当前未执行）
#install.packages('data.tree')
#install.packages('DiagrammeR')
#install.packages('./randomForestSRC_3.2.3.tar.gz')
setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/55/')
# Mime包的本地安装命令（当前未执行）
#devtools::install_local('./Mime-main.zip')
library(Mime1)
library(tibble)
load('./features_train_enet.Rdata')
load('./features_vali.Rdata')
library(data.table)

suv=fread('./TCGA-KIRC.survival.tsv.gz',check.names = F,data.table = F)
colnames(suv)[4]='Sample'

suv=suv[!duplicated(suv$Sample),]
rownames(suv)=suv$Sample

## 按TCGA样本编号匹配训练集病理特征与生存信息
sss=intersect(rownames(features_train3),rownames(suv))
suv_train=suv[sss,]
features_train3=features_train3[sss,]

suv_train=suv_train[,c(2,3)]
colnames(suv_train)=c("OS.time",'OS')

ssss=intersect(rownames(features_vali),rownames(suv))
suv_test=suv[ssss,]
features_vali=features_vali[ssss,colnames(features_train3)]
suv_test=suv_test[,c(2,3)]
colnames(suv_test)=c("OS.time",'OS')



## 构建训练集和内部验证集的生存分析数据

train=cbind(suv_train,features_train3)
# 将训练集转换为Mime分析所需的数据格式
train=rownames_to_column(train,var='ID')
colnames(train) <- gsub("_", ".", colnames(train))
features=colnames(train)


vali=cbind(suv_test,features_vali)
vali=rownames_to_column(vali,var='ID')
colnames(vali) <- gsub("_", ".", colnames(vali))

# 检查训练集与验证集共有的变量名称
ss=intersect(colnames(train),colnames(vali))


## 训练集单因素Cox比例风险回归分析
outTab=data.frame()
pFilter=0.05
sigfeature=c("OS.time","OS")
library(survival)
for(feature in colnames(train[,4:ncol(train)])){
  set.seed(123456)
  cox=coxph(Surv(OS.time, OS) ~ train[,feature], data = train)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigfeature=c(sigfeature,feature)
    outTab=rbind(outTab,
                 cbind(feature=feature,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP))
    print(coxP)
  }
}

outTab$feature=gsub('\\.','_',outTab$feature)
write.csv(outTab,file ='unicox.csv')

uni_feature=outTab$feature
save(uni_feature,file ='unifeature.Rdata')

# 加载森林图绘制所需R包
library(forestplot)

outTab$HR=as.numeric(outTab$HR)
outTab$HR.95L=as.numeric(outTab$HR.95L)
outTab$HR.95H=as.numeric(outTab$HR.95H)

str(outTab)

### 单因素Cox回归结果可视化 -------------------------

# 加载森林图绘制所需R包
library(forestplot)

# outTab为单因素Cox回归结果数据框
# 包含特征名称、风险比、95%置信区间及P值

# 整理森林图所需的HR及95%置信区间数据
forest_data <- data.frame(
  mean = outTab$HR, 
  lower = outTab$HR.95L, 
  upper = outTab$HR.95H, 
  feature = outTab$feature
)

# 生成森林图中显示的特征名称及数值标签
forest_data$feature <- as.character(forest_data$feature)
forest_data$HR <- sprintf("%.2f", forest_data$mean)
forest_data$CI <- paste0("(", sprintf("%.2f", forest_data$lower), ", ", sprintf("%.2f", forest_data$upper), ")")

# 合并HR与95%置信区间标签
forest_data$HR_CI <- paste(forest_data$HR, forest_data$CI, sep = " ")

# 根据HR是否大于1设置森林图点的颜色
forest_data$point_color <- ifelse(forest_data$mean > 1, "red", "blue")

# 绘制单因素Cox回归森林图
forestplot(
  labeltext = list(forest_data$feature, forest_data$HR_CI), # 显示特征名称及HR（95% CI）
  mean = forest_data$mean, # HR估计值
  lower = forest_data$lower, # 95%置信区间下限
  upper = forest_data$upper, # 95%置信区间上限
  new_page = TRUE, 
  xlab = "Hazard Ratio (95% CI)", 
  title = "Forest Plot of Hazard Ratios", 
  boxsize = 0.2, 
  lineheight = unit(8, "mm"), 
  # 可选的森林图配色参数（当前未启用）
  # col = fpColors(box = "lightgray", line = "darkblue", summary = "red"), # 森林图方框、线条及汇总行的可选颜色
  # 根据HR设置点颜色的备选绘图方案（当前未启用）
  graphwidth = unit(6, "cm"), 
  zero = 1, # HR=1的无效应参考线
  
)



### 基于Mime的多算法预后模型构建与比较 -------------------------
list_train_vali_Data=list(train,vali)
names(list_train_vali_Data)=c('Dataset1','Dataset2')
# 查看模型开发函数的帮助文档（当前未执行）
#?ML.Dev.Prog.Sig
res <- Mime1::ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Dataset1,
                              list_train_vali_Data = list_train_vali_Data,
                              unicox.filter.for.candi = T,
                              unicox_p_cutoff = 0.05,
                              candidate_genes =features,
                              mode = 'all', ## 运行全部支持的单算法及组合算法模型
                              nodesize =5,
                              seed = 123 
)



save(res,file ='ML_res.Rdata')

load('ML_res.Rdata')

cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order = names(list_train_vali_Data),
               width = 0.35
)

pdf(
  file = "Cindex_all_models_TCGA.pdf",
  width = 12,
  height = 24
)

cindex_dis_all(
  res,
  validate_set = names(list_train_vali_Data)[-1],
  order = NULL,
  width = 0.35
)

dev.off()




## 查看指定预后模型在训练集和内部验证集中的C-index
setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/55/model/')
cindex_dis_select(res,
                  model="RSF + plsRcox",### 当前代码指定展示的模型
                  order= names(list_train_vali_Data))


survplot <- vector("list",2) 

# 绘制指定模型在两个数据集中的Kaplan–Meier生存曲线
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "RSF + plsRcox",## 当前代码指定展示的模型
                              dataset = names(list_train_vali_Data)[i],
                              # 可选的生存曲线颜色设置（当前未启用）
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Years")
  )
}
aplot::plot_list(gglist=survplot,ncol=2)


all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set = "Dataset1",
            order = names(list_train_vali_Data),
            width = 0.35,
            year = 1)

roc_vis(all.auc.3y,
        model_name = "RSF + plsRcox",
        dataset = "Dataset1",
        order = "Dataset1",
        anno_position = c(0.65, 0.55),
        year = 3)
### 提取指定预后模型的特征信息 -----
res2=res$ml.res$`RSF + plsRcox`  ### 当前代码提取的模型

class(res2)
str(res2, max.level = 2)
final_feature <- colnames(res2$dataX)
final_feature
# 另一模型结果的提取示例（当前未启用）
#res3=res$ml.res$`Lasso + plsRcox`  ### 备选模型
final_feature <- gsub("\\.", "_", final_feature)

save(final_feature, file = "final_feature.Rdata")

write.csv(
  data.frame(feature = final_feature),
  file = "final_feature.csv",
  row.names = FALSE
)




# 基于模型正则化参数提取系数的备选代码（当前未执行）
#coefficients <- coef(res2, s = "lambda.min")
#coefficients_matrix <- as.matrix(coefficients)
#write.csv(coefficients_matrix,file ='coef_matrix.csv')

#index <- which(coefficients  != 0)
#index
#final_feature=row.names(coefficients)[index]
#final_feature<- gsub("\\.", "_", final_feature)
#final_feature

final_feature
save(final_feature,file ='final_feature.Rdata')


### 提取指定预后模型的训练集风险评分

riskscore=res$riskscore$`RSF + plsRcox`$Dataset1  ### 提取Dataset1训练集的风险评分
save(riskscore,file ='riskscore.Rdata')

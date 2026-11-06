#harmony整合，GSE210038,7个肿瘤样本，2个正常样本；GSE156632:5个正常，7个肿瘤样本
rm(list=ls())
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Rcpp)
library(scDblFinder)
library(SeuratData)
library(DT)
###GSE210038
####T1_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415686.data <- read.table(file="T1.txt")
GSM6415686_1 <- CreateSeuratObject(counts = GSM6415686.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415686@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415686_1)
GSM6415686_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415686_1$scDblFinder.class)
# GSM6415686@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415686_2 <-as.Seurat(GSM6415686_1)
GSM6415686 = GSM6415686_2[,GSM6415686_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415686
GSM6415686@assays$RNA@counts[1:5,1:5]
GSM6415686[["percent.mt"]] <- PercentageFeatureSet(GSM6415686, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415686@assays$RNA)) 
HB.genes <- rownames(GSM6415686@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415686[["percent.HB"]]<-PercentageFeatureSet(GSM6415686, features=HB.genes) 
col.num <- length(levels(GSM6415686@active.ident))

violin <- VlnPlot(GSM6415686,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415686_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415686 <- subset(GSM6415686, subset =nFeature_RNA <5000 & nCount_RNA <20000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415686@meta.data$group <- "Tumor"
GSM6415686<- RenameCells(GSM6415686, add.cell.id = "GSM6415686")
GSM6415686@meta.data$orig.ident <- "GSM6415686"
GSM6415686@meta.data$stage <- 4
GSM6415686@meta.data$age <- 69
GSM6415686@meta.data$gender <- 1
GSM6415686@meta.data$sample <- "T1"

####GSM6415687_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415687.data <- read.table(file="T2.txt")
GSM6415687_1 <- CreateSeuratObject(counts = GSM6415687.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415687@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415687_1)
GSM6415687_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415687_1$scDblFinder.class)
# GSM6415687@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415687_2 <-as.Seurat(GSM6415687_1)
GSM6415687 = GSM6415687_2[,GSM6415687_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415687
GSM6415687@assays$RNA@counts[1:5,1:5]
GSM6415687[["percent.mt"]] <- PercentageFeatureSet(GSM6415687, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415687@assays$RNA)) 
HB.genes <- rownames(GSM6415687@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415687[["percent.HB"]]<-PercentageFeatureSet(GSM6415687, features=HB.genes) 
col.num <- length(levels(GSM6415687@active.ident))

violin <- VlnPlot(GSM6415687,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415687_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415687 <- subset(GSM6415687, subset =nFeature_RNA <4000 & nCount_RNA <20000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415687@meta.data$group <- "Tumor"
GSM6415687<- RenameCells(GSM6415687, add.cell.id = "GSM6415687")
GSM6415687@meta.data$orig.ident <- "GSM6415687"
GSM6415687@meta.data$stage <- 3
GSM6415687@meta.data$age <- 84
GSM6415687@meta.data$gender <- 2
GSM6415687@meta.data$sample <- "T2"

####GSM6415688_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415688.data <- read.table(file="T3.txt")
GSM6415688_1 <- CreateSeuratObject(counts = GSM6415688.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415688@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415688_1)
GSM6415688_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415688_1$scDblFinder.class)
# GSM6415688@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415688_2 <-as.Seurat(GSM6415688_1)
GSM6415688 = GSM6415688_2[,GSM6415688_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415688
GSM6415688@assays$RNA@counts[1:5,1:5]
GSM6415688[["percent.mt"]] <- PercentageFeatureSet(GSM6415688, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415688@assays$RNA)) 
HB.genes <- rownames(GSM6415688@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415688[["percent.HB"]]<-PercentageFeatureSet(GSM6415688, features=HB.genes) 
col.num <- length(levels(GSM6415688@active.ident))

violin <- VlnPlot(GSM6415688,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415688_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415688 <- subset(GSM6415688, subset =nFeature_RNA <4000 & nCount_RNA <20000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415688@meta.data$group <- "Tumor"
GSM6415688<- RenameCells(GSM6415688, add.cell.id = "GSM6415688")
GSM6415688@meta.data$orig.ident <- "GSM6415688"
GSM6415688@meta.data$stage <- 4
GSM6415688@meta.data$age <- 58
GSM6415688@meta.data$gender <- 2
GSM6415688@meta.data$sample <- "T3"


####GSM6415689_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415689.data <- read.table(file="T4.txt")
GSM6415689_1 <- CreateSeuratObject(counts = GSM6415689.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415689@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415689_1)
GSM6415689_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415689_1$scDblFinder.class)
# GSM6415689@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415689_2 <-as.Seurat(GSM6415689_1)
GSM6415689 = GSM6415689_2[,GSM6415689_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415689
GSM6415689@assays$RNA@counts[1:5,1:5]
GSM6415689[["percent.mt"]] <- PercentageFeatureSet(GSM6415689, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415689@assays$RNA)) 
HB.genes <- rownames(GSM6415689@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415689[["percent.HB"]]<-PercentageFeatureSet(GSM6415689, features=HB.genes) 
col.num <- length(levels(GSM6415689@active.ident))

violin <- VlnPlot(GSM6415689,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415689_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415689 <- subset(GSM6415689, subset =nFeature_RNA <4000 & nCount_RNA <15000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415689@meta.data$group <- "Tumor"
GSM6415689<- RenameCells(GSM6415689, add.cell.id = "GSM6415689")
GSM6415689@meta.data$orig.ident <- "GSM6415689"
GSM6415689@meta.data$stage <- 3
GSM6415689@meta.data$age <- 70
GSM6415689@meta.data$gender <- 1
GSM6415689@meta.data$sample <- "T4"

####T5_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415690.data <- read.table(file="T5.txt")
GSM6415690_1 <- CreateSeuratObject(counts = GSM6415690.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415690@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415690_1)
GSM6415690_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415690_1$scDblFinder.class)
# GSM6415690@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415690_2 <-as.Seurat(GSM6415690_1)
GSM6415690 = GSM6415690_2[,GSM6415690_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415690
GSM6415690@assays$RNA@counts[1:5,1:5]
GSM6415690[["percent.mt"]] <- PercentageFeatureSet(GSM6415690, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415690@assays$RNA)) 
HB.genes <- rownames(GSM6415690@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415690[["percent.HB"]]<-PercentageFeatureSet(GSM6415690, features=HB.genes) 
col.num <- length(levels(GSM6415690@active.ident))

violin <- VlnPlot(GSM6415690,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415690_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415690 <- subset(GSM6415690, subset =nFeature_RNA <4000 & nCount_RNA <15000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415690@meta.data$group <- "Tumor"
GSM6415690<- RenameCells(GSM6415690, add.cell.id = "GSM6415690")
GSM6415690@meta.data$orig.ident <- "GSM6415690"
GSM6415690@meta.data$stage <- 2
GSM6415690@meta.data$age <- 63
GSM6415690@meta.data$gender <- 1
GSM6415690@meta.data$sample <- "T5"

####GSM6415691_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415691.data <- read.table(file="T6.txt")
GSM6415691_1 <- CreateSeuratObject(counts = GSM6415691.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415691@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415691_1)
GSM6415691_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415691_1$scDblFinder.class)
# GSM6415691@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415691_2 <-as.Seurat(GSM6415691_1)
GSM6415691 = GSM6415691_2[,GSM6415691_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415691
GSM6415691@assays$RNA@counts[1:5,1:5]
GSM6415691[["percent.mt"]] <- PercentageFeatureSet(GSM6415691, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415691@assays$RNA)) 
HB.genes <- rownames(GSM6415691@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415691[["percent.HB"]]<-PercentageFeatureSet(GSM6415691, features=HB.genes) 
col.num <- length(levels(GSM6415691@active.ident))

violin <- VlnPlot(GSM6415691,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415691_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415691 <- subset(GSM6415691, subset =nFeature_RNA <4000 & nCount_RNA <15000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415691@meta.data$group <- "Tumor"
GSM6415691<- RenameCells(GSM6415691, add.cell.id = "GSM6415691")
GSM6415691@meta.data$orig.ident <- "GSM6415691"
GSM6415691@meta.data$stage <- 1
GSM6415691@meta.data$age <- 48
GSM6415691@meta.data$gender <- 1
GSM6415691@meta.data$sample <- "T6"

####GSM6415692_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415692.data <- read.table(file="T7.txt")
GSM6415692_1 <- CreateSeuratObject(counts = GSM6415692.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415692@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415692_1)
GSM6415692_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415692_1$scDblFinder.class)
# GSM6415692@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415692_2 <-as.Seurat(GSM6415692_1)
GSM6415692 = GSM6415692_2[,GSM6415692_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415692
GSM6415692@assays$RNA@counts[1:5,1:5]
GSM6415692[["percent.mt"]] <- PercentageFeatureSet(GSM6415692, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415692@assays$RNA)) 
HB.genes <- rownames(GSM6415692@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415692[["percent.HB"]]<-PercentageFeatureSet(GSM6415692, features=HB.genes) 
col.num <- length(levels(GSM6415692@active.ident))

violin <- VlnPlot(GSM6415692,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415692_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415692 <- subset(GSM6415692, subset =nFeature_RNA <4000 & nCount_RNA <15000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415692@meta.data$group <- "Tumor"
GSM6415692<- RenameCells(GSM6415692, add.cell.id = "GSM6415692")
GSM6415692@meta.data$orig.ident <- "GSM6415692"
GSM6415692@meta.data$stage <- 3
GSM6415692@meta.data$age <- 69
GSM6415692@meta.data$gender <- 2
GSM6415692@meta.data$sample <- "T7"


####GSM6415693_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415693.data <- read.table(file="N1.txt")
GSM6415693_1 <- CreateSeuratObject(counts = GSM6415693.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415693@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415693_1)
GSM6415693_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415693_1$scDblFinder.class)
# GSM6415693@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415693_2 <-as.Seurat(GSM6415693_1)
GSM6415693 = GSM6415693_2[,GSM6415693_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415693
GSM6415693@assays$RNA@counts[1:5,1:5]
GSM6415693[["percent.mt"]] <- PercentageFeatureSet(GSM6415693, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415693@assays$RNA)) 
HB.genes <- rownames(GSM6415693@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415693[["percent.HB"]]<-PercentageFeatureSet(GSM6415693, features=HB.genes) 
col.num <- length(levels(GSM6415693@active.ident))

violin <- VlnPlot(GSM6415693,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415693_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415693 <- subset(GSM6415693, subset =nFeature_RNA <2000 & nCount_RNA <2500 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415693@meta.data$group <- "Normal"
GSM6415693<- RenameCells(GSM6415693, add.cell.id = "GSM6415693")
GSM6415693@meta.data$orig.ident <- "GSM6415693"
GSM6415693@meta.data$sample <- "N1"
####GSM6415694_______________________________________________________________________________
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE210038/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM6415694.data <- read.table(file="N2.txt")
GSM6415694_1 <- CreateSeuratObject(counts = GSM6415694.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM6415694@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM6415694_1)
GSM6415694_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM6415694_1$scDblFinder.class)
# GSM6415694@colData@listData[["scDblFinder.class"]],数据存放位置
GSM6415694_2 <-as.Seurat(GSM6415694_1)
GSM6415694 = GSM6415694_2[,GSM6415694_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM6415694
GSM6415694@assays$RNA@counts[1:5,1:5]
GSM6415694[["percent.mt"]] <- PercentageFeatureSet(GSM6415694, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM6415694@assays$RNA)) 
HB.genes <- rownames(GSM6415694@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM6415694[["percent.HB"]]<-PercentageFeatureSet(GSM6415694, features=HB.genes) 
col.num <- length(levels(GSM6415694@active.ident))

violin <- VlnPlot(GSM6415694,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM6415694_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM6415694 <- subset(GSM6415694, subset =nFeature_RNA <4000 & nCount_RNA <15000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM6415694@meta.data$group <- "Normal"
GSM6415694<- RenameCells(GSM6415694, add.cell.id = "GSM6415694")
GSM6415694@meta.data$orig.ident <- "GSM6415694"
GSM6415694@meta.data$sample <- "N2"


##############GSE156632,5个正常样本，7个肿瘤样本
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735364/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735364.data <- read.csv(file="GSM4735364.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735364.data <- GSM4735364.data[, -1]
GSM4735364.data <- GSM4735364.data[!duplicated(GSM4735364.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735364.data) <- GSM4735364.data$Symbol                # 设置 Symbol 为行名
GSM4735364.data <- GSM4735364.data[, colnames(GSM4735364.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735364.data) <- gsub("_", "-", rownames(GSM4735364.data))###将——改为_
GSM4735364_1 <- CreateSeuratObject(counts = GSM4735364.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735364@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735364_1)
GSM4735364_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735364_1$scDblFinder.class)
# GSM4735364@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735364_2 <-as.Seurat(GSM4735364_1)
GSM4735364 = GSM4735364_2[,GSM4735364_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735364
GSM4735364@assays$RNA@counts[1:5,1:5]
GSM4735364[["percent.mt"]] <- PercentageFeatureSet(GSM4735364, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735364@assays$RNA)) 
HB.genes <- rownames(GSM4735364@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735364[["percent.HB"]]<-PercentageFeatureSet(GSM4735364, features=HB.genes) 
col.num <- length(levels(GSM4735364@active.ident))

violin <- VlnPlot(GSM4735364,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735364_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735364 <- subset(GSM4735364, subset =nFeature_RNA <3000 & nCount_RNA <15000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM4735364@meta.data$group <- "Tumor"
GSM4735364<- RenameCells(GSM4735364, add.cell.id = "GSM4735364")
GSM4735364@meta.data$orig.ident <- "GSM4735364"
GSM4735364@meta.data$stage <- 3
GSM4735364@meta.data$age <- 68
GSM4735364@meta.data$gender <- 1
GSM4735364@meta.data$sample <- "T8"



setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735365/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735365.data <- read.csv(file="GSM4735365.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735365.data <- GSM4735365.data[, -1]
GSM4735365.data <- GSM4735365.data[!duplicated(GSM4735365.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735365.data) <- GSM4735365.data$Symbol                # 设置 Symbol 为行名
GSM4735365.data <- GSM4735365.data[, colnames(GSM4735365.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735365.data) <- gsub("_", "-", rownames(GSM4735365.data))###将——改为_
GSM4735365_1 <- CreateSeuratObject(counts = GSM4735365.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735365@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735365_1)
GSM4735365_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735365_1$scDblFinder.class)
# GSM4735365@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735365_2 <-as.Seurat(GSM4735365_1)
GSM4735365 = GSM4735365_2[,GSM4735365_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735365
GSM4735365@assays$RNA@counts[1:5,1:5]
GSM4735365[["percent.mt"]] <- PercentageFeatureSet(GSM4735365, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735365@assays$RNA)) 
HB.genes <- rownames(GSM4735365@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735365[["percent.HB"]]<-PercentageFeatureSet(GSM4735365, features=HB.genes) 
col.num <- length(levels(GSM4735365@active.ident))

violin <- VlnPlot(GSM4735365,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735365_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735365 <- subset(GSM4735365, subset =nFeature_RNA <3000 & nCount_RNA <10000 & percent.mt < 15& percent.HB < 1 )

#在meta.data中添加group这一项
GSM4735365@meta.data$group <- "Normal"
GSM4735365<- RenameCells(GSM4735365, add.cell.id = "GSM4735365")
GSM4735365@meta.data$orig.ident <- "GSM4735365"
GSM4735365@meta.data$sample <- "N3"



######
setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735366/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735366.data <- read.csv(file="GSM4735366.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735366.data <- GSM4735366.data[, -1]
GSM4735366.data <- GSM4735366.data[!duplicated(GSM4735366.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735366.data) <- GSM4735366.data$Symbol                # 设置 Symbol 为行名
GSM4735366.data <- GSM4735366.data[, colnames(GSM4735366.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735366.data) <- gsub("_", "-", rownames(GSM4735366.data))###将——改为_
GSM4735366_1 <- CreateSeuratObject(counts = GSM4735366.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735366@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735366_1)
GSM4735366_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735366_1$scDblFinder.class)
# GSM4735366@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735366_2 <-as.Seurat(GSM4735366_1)
GSM4735366 = GSM4735366_2[,GSM4735366_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735366
GSM4735366@assays$RNA@counts[1:5,1:5]
GSM4735366[["percent.mt"]] <- PercentageFeatureSet(GSM4735366, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735366@assays$RNA)) 
HB.genes <- rownames(GSM4735366@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735366[["percent.HB"]]<-PercentageFeatureSet(GSM4735366, features=HB.genes) 
col.num <- length(levels(GSM4735366@active.ident))

violin <- VlnPlot(GSM4735366,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735366_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735366 <- subset(GSM4735366, subset =nFeature_RNA <5000 & nCount_RNA <20000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM4735366@meta.data$group <- "Tumor"
GSM4735366<- RenameCells(GSM4735366, add.cell.id = "GSM4735366")
GSM4735366@meta.data$orig.ident <- "GSM4735366"
GSM4735366@meta.data$stage <- 2
GSM4735366@meta.data$age <- 51
GSM4735366@meta.data$gender <- 1
GSM4735366@meta.data$sample <- "T9"

setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735367/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735367.data <- read.csv(file="GSM4735367.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735367.data <- GSM4735367.data[, -1]
GSM4735367.data <- GSM4735367.data[!duplicated(GSM4735367.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735367.data) <- GSM4735367.data$Symbol                # 设置 Symbol 为行名
GSM4735367.data <- GSM4735367.data[, colnames(GSM4735367.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735367.data) <- gsub("_", "-", rownames(GSM4735367.data))###将——改为_
GSM4735367_1 <- CreateSeuratObject(counts = GSM4735367.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735367@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735367_1)
GSM4735367_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735367_1$scDblFinder.class)
# GSM4735367@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735367_2 <-as.Seurat(GSM4735367_1)
GSM4735367 = GSM4735367_2[,GSM4735367_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735367
GSM4735367@assays$RNA@counts[1:5,1:5]
GSM4735367[["percent.mt"]] <- PercentageFeatureSet(GSM4735367, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735367@assays$RNA)) 
HB.genes <- rownames(GSM4735367@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735367[["percent.HB"]]<-PercentageFeatureSet(GSM4735367, features=HB.genes) 
col.num <- length(levels(GSM4735367@active.ident))

violin <- VlnPlot(GSM4735367,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735367_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735367 <- subset(GSM4735367, subset =nFeature_RNA <5000 & nCount_RNA <20000 & percent.mt < 15 & percent.HB < 1)

#在meta.data中添加group这一项
GSM4735367@meta.data$group <- "Normal"
GSM4735367<- RenameCells(GSM4735367, add.cell.id = "GSM4735367")
GSM4735367@meta.data$orig.ident <- "GSM4735367"
GSM4735367@meta.data$sample <- "N4"

setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735368/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735368.data <- read.csv(file="GSM4735368.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735368.data <- GSM4735368.data[, -1]
GSM4735368.data <- GSM4735368.data[!duplicated(GSM4735368.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735368.data) <- GSM4735368.data$Symbol                # 设置 Symbol 为行名
GSM4735368.data <- GSM4735368.data[, colnames(GSM4735368.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735368.data) <- gsub("_", "-", rownames(GSM4735368.data))###将——改为_
GSM4735368_1 <- CreateSeuratObject(counts = GSM4735368.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735368@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735368_1)
GSM4735368_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735368_1$scDblFinder.class)
# GSM4735368@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735368_2 <-as.Seurat(GSM4735368_1)
GSM4735368 = GSM4735368_2[,GSM4735368_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735368
GSM4735368@assays$RNA@counts[1:5,1:5]
GSM4735368[["percent.mt"]] <- PercentageFeatureSet(GSM4735368, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735368@assays$RNA)) 
HB.genes <- rownames(GSM4735368@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735368[["percent.HB"]]<-PercentageFeatureSet(GSM4735368, features=HB.genes) 
col.num <- length(levels(GSM4735368@active.ident))

violin <- VlnPlot(GSM4735368,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735368_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735368 <- subset(GSM4735368, subset =nFeature_RNA <4000 & nCount_RNA <15000 & percent.mt < 10& percent.HB < 1 )

#在meta.data中添加group这一项
GSM4735368@meta.data$group <- "Tumor"
GSM4735368<- RenameCells(GSM4735368, add.cell.id = "GSM4735368")
GSM4735368@meta.data$orig.ident <- "GSM4735368"
GSM4735368@meta.data$stage <- 2
GSM4735368@meta.data$age <- 68
GSM4735368@meta.data$gender <- 1
GSM4735368@meta.data$sample <- "T10"



setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735369/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735369.data <- read.csv(file="GSM4735369.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735369.data <- GSM4735369.data[, -1]
GSM4735369.data <- GSM4735369.data[!duplicated(GSM4735369.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735369.data) <- GSM4735369.data$Symbol                # 设置 Symbol 为行名
GSM4735369.data <- GSM4735369.data[, colnames(GSM4735369.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735369.data) <- gsub("_", "-", rownames(GSM4735369.data))###将——改为_
GSM4735369_1 <- CreateSeuratObject(counts = GSM4735369.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735369@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735369_1)
GSM4735369_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735369_1$scDblFinder.class)
# GSM4735369@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735369_2 <-as.Seurat(GSM4735369_1)
GSM4735369 = GSM4735369_2[,GSM4735369_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735369
GSM4735369@assays$RNA@counts[1:5,1:5]
GSM4735369[["percent.mt"]] <- PercentageFeatureSet(GSM4735369, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735369@assays$RNA)) 
HB.genes <- rownames(GSM4735369@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735369[["percent.HB"]]<-PercentageFeatureSet(GSM4735369, features=HB.genes) 
col.num <- length(levels(GSM4735369@active.ident))

violin <- VlnPlot(GSM4735369,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735369_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735369 <- subset(GSM4735369, subset =nFeature_RNA <2000 & nCount_RNA <10000 & percent.mt < 10& percent.HB < 1 )

#在meta.data中添加group这一项
GSM4735369@meta.data$group <- "Normal"
GSM4735369<- RenameCells(GSM4735369, add.cell.id = "GSM4735369")
GSM4735369@meta.data$orig.ident <- "GSM4735369"
GSM4735369@meta.data$sample <- "N5"


setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735370/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735370.data <- read.csv(file="GSM4735370.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735370.data <- GSM4735370.data[, -1]
GSM4735370.data <- GSM4735370.data[!duplicated(GSM4735370.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735370.data) <- GSM4735370.data$Symbol                # 设置 Symbol 为行名
GSM4735370.data <- GSM4735370.data[, colnames(GSM4735370.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735370.data) <- gsub("_", "-", rownames(GSM4735370.data))###将——改为_
GSM4735370_1 <- CreateSeuratObject(counts = GSM4735370.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735370@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735370_1)
GSM4735370_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735370_1$scDblFinder.class)
# GSM4735370@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735370_2 <-as.Seurat(GSM4735370_1)
GSM4735370 = GSM4735370_2[,GSM4735370_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735370
GSM4735370@assays$RNA@counts[1:5,1:5]
GSM4735370[["percent.mt"]] <- PercentageFeatureSet(GSM4735370, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735370@assays$RNA)) 
HB.genes <- rownames(GSM4735370@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735370[["percent.HB"]]<-PercentageFeatureSet(GSM4735370, features=HB.genes) 
col.num <- length(levels(GSM4735370@active.ident))

violin <- VlnPlot(GSM4735370,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735370_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735370 <- subset(GSM4735370, subset =nFeature_RNA <6000 & nCount_RNA <50000 & percent.mt < 10 & percent.HB < 1)

#在meta.data中添加group这一项
GSM4735370@meta.data$group <- "Tumor"
GSM4735370<- RenameCells(GSM4735370, add.cell.id = "GSM4735370")
GSM4735370@meta.data$orig.ident <- "GSM4735370"
GSM4735370@meta.data$stage <- 3
GSM4735370@meta.data$age <- 59
GSM4735370@meta.data$gender <- 2
GSM4735370@meta.data$sample <- "T11"


setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735371/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735371.data <- read.csv(file="GSM4735371.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735371.data <- GSM4735371.data[, -1]
GSM4735371.data <- GSM4735371.data[!duplicated(GSM4735371.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735371.data) <- GSM4735371.data$Symbol                # 设置 Symbol 为行名
GSM4735371.data <- GSM4735371.data[, colnames(GSM4735371.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735371.data) <- gsub("_", "-", rownames(GSM4735371.data))###将——改为_
GSM4735371_1 <- CreateSeuratObject(counts = GSM4735371.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735371@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735371_1)
GSM4735371_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735371_1$scDblFinder.class)
# GSM4735371@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735371_2 <-as.Seurat(GSM4735371_1)
GSM4735371 = GSM4735371_2[,GSM4735371_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735371
GSM4735371@assays$RNA@counts[1:5,1:5]
GSM4735371[["percent.mt"]] <- PercentageFeatureSet(GSM4735371, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735371@assays$RNA)) 
HB.genes <- rownames(GSM4735371@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735371[["percent.HB"]]<-PercentageFeatureSet(GSM4735371, features=HB.genes) 
col.num <- length(levels(GSM4735371@active.ident))

violin <- VlnPlot(GSM4735371,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735371_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735371 <- subset(GSM4735371, subset =nFeature_RNA <2000 & nCount_RNA <10000 & percent.mt < 10& percent.HB < 1 )

#在meta.data中添加group这一项
GSM4735371@meta.data$group <- "Normal"
GSM4735371<- RenameCells(GSM4735371, add.cell.id = "GSM4735371")
GSM4735371@meta.data$orig.ident <- "GSM4735371"
GSM4735371@meta.data$sample <- "N6"

setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735372/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735372.data <- read.csv(file="GSM4735372.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735372.data <- GSM4735372.data[, -1]
GSM4735372.data <- GSM4735372.data[!duplicated(GSM4735372.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735372.data) <- GSM4735372.data$Symbol                # 设置 Symbol 为行名
GSM4735372.data <- GSM4735372.data[, colnames(GSM4735372.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735372.data) <- gsub("_", "-", rownames(GSM4735372.data))###将——改为_
GSM4735372_1 <- CreateSeuratObject(counts = GSM4735372.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735372@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735372_1)
GSM4735372_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735372_1$scDblFinder.class)
# GSM4735372@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735372_2 <-as.Seurat(GSM4735372_1)
GSM4735372 = GSM4735372_2[,GSM4735372_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735372
GSM4735372@assays$RNA@counts[1:5,1:5]
GSM4735372[["percent.mt"]] <- PercentageFeatureSet(GSM4735372, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735372@assays$RNA)) 
HB.genes <- rownames(GSM4735372@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735372[["percent.HB"]]<-PercentageFeatureSet(GSM4735372, features=HB.genes) 
col.num <- length(levels(GSM4735372@active.ident))

violin <- VlnPlot(GSM4735372,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735372_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735372 <- subset(GSM4735372, subset =nFeature_RNA <4000 & nCount_RNA <15000 & percent.mt < 10 & percent.HB < 1)

#在meta.data中添加group这一项
GSM4735372@meta.data$group <- "Tumor"
GSM4735372<- RenameCells(GSM4735372, add.cell.id = "GSM4735372")
GSM4735372@meta.data$orig.ident <- "GSM4735372"
GSM4735372@meta.data$stage <- 2
GSM4735372@meta.data$age <- 41
GSM4735372@meta.data$gender <- 2
GSM4735372@meta.data$sample <- "T12"

setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735373/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735373.data <- read.csv(file="GSM4735373.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735373.data <- GSM4735373.data[, -1]
GSM4735373.data <- GSM4735373.data[!duplicated(GSM4735373.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735373.data) <- GSM4735373.data$Symbol                # 设置 Symbol 为行名
GSM4735373.data <- GSM4735373.data[, colnames(GSM4735373.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735373.data) <- gsub("_", "-", rownames(GSM4735373.data))###将——改为_
GSM4735373_1 <- CreateSeuratObject(counts = GSM4735373.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735373@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735373_1)
GSM4735373_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735373_1$scDblFinder.class)
# GSM4735373@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735373_2 <-as.Seurat(GSM4735373_1)
GSM4735373 = GSM4735373_2[,GSM4735373_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735373
GSM4735373@assays$RNA@counts[1:5,1:5]
GSM4735373[["percent.mt"]] <- PercentageFeatureSet(GSM4735373, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735373@assays$RNA)) 
HB.genes <- rownames(GSM4735373@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735373[["percent.HB"]]<-PercentageFeatureSet(GSM4735373, features=HB.genes) 
HB.genes <- c("HBA1", "HBA2", "HBB")
col.num <- length(levels(GSM4735373@active.ident))

violin <- VlnPlot(GSM4735373,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735373_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735373 <- subset(GSM4735373, subset =nFeature_RNA <2000 & nCount_RNA <5000 & percent.mt < 10 & percent.HB < 1)

#在meta.data中添加group这一项
GSM4735373@meta.data$group <- "Normal"
GSM4735373<- RenameCells(GSM4735373, add.cell.id = "GSM4735373")
GSM4735373@meta.data$orig.ident <- "GSM4735373"
GSM4735373@meta.data$sample <- "N7"

setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735374/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735374.data <- read.csv(file="GSM4735374.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735374.data <- GSM4735374.data[, -1]
GSM4735374.data <- GSM4735374.data[!duplicated(GSM4735374.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735374.data) <- GSM4735374.data$Symbol                # 设置 Symbol 为行名
GSM4735374.data <- GSM4735374.data[, colnames(GSM4735374.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735374.data) <- gsub("_", "-", rownames(GSM4735374.data))###将——改为_
GSM4735374_1 <- CreateSeuratObject(counts = GSM4735374.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735374@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735374_1)
GSM4735374_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735374_1$scDblFinder.class)
# GSM4735374@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735374_2 <-as.Seurat(GSM4735374_1)
GSM4735374 = GSM4735374_2[,GSM4735374_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735374
GSM4735374@assays$RNA@counts[1:5,1:5]
GSM4735374[["percent.mt"]] <- PercentageFeatureSet(GSM4735374, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735374@assays$RNA)) 
HB.genes <- rownames(GSM4735374@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735374[["percent.HB"]]<-PercentageFeatureSet(GSM4735374, features=HB.genes) 
col.num <- length(levels(GSM4735374@active.ident))

violin <- VlnPlot(GSM4735374,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735374_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735374 <- subset(GSM4735374, subset =nFeature_RNA <4000 & nCount_RNA <20000 & percent.mt < 10& percent.HB < 1 )

#在meta.data中添加group这一项
GSM4735374@meta.data$group <- "Tumor"
GSM4735374<- RenameCells(GSM4735374, add.cell.id = "GSM4735374")
GSM4735374@meta.data$orig.ident <- "GSM4735374"
GSM4735374@meta.data$stage <- 2
GSM4735374@meta.data$age <- 71
GSM4735374@meta.data$gender <- 2
GSM4735374@meta.data$sample <- "T13"

setwd("/home/stu-xt/ccRcc_24.9.2/samples/GSE156632/GSM4735375/")#每一个样本的1.创建seurat,2.质控，3.标准化，寻找高变gene
GSM4735375.data <- read.csv(file="GSM4735375.csv",sep=",",header=T,check.names=F)
# 去掉第一列
GSM4735375.data <- GSM4735375.data[, -1]
GSM4735375.data <- GSM4735375.data[!duplicated(GSM4735375.data$Symbol), ]  # 去除重复 Symbol
rownames(GSM4735375.data) <- GSM4735375.data$Symbol                # 设置 Symbol 为行名
GSM4735375.data <- GSM4735375.data[, colnames(GSM4735375.data) != "Symbol"]  # 删除 Symbol 列
rownames(GSM4735375.data) <- gsub("_", "-", rownames(GSM4735375.data))###将——改为_
GSM4735375_1 <- CreateSeuratObject(counts = GSM4735375.data,min.cells = 3, min.features = 200)
#scDblFinder
#sce <-GSM4735375@assays[["RNA"]]@counts
sce <- as.SingleCellExperiment(GSM4735375_1)
GSM4735375_1 <- scDblFinder(sce, dbr=0.1) 
table(GSM4735375_1$scDblFinder.class)
# GSM4735375@colData@listData[["scDblFinder.class"]],数据存放位置
GSM4735375_2 <-as.Seurat(GSM4735375_1)
GSM4735375 = GSM4735375_2[,GSM4735375_2@meta.data[["scDblFinder.class"]] %in% c("singlet")]
GSM4735375
GSM4735375@assays$RNA@counts[1:5,1:5]
GSM4735375[["percent.mt"]] <- PercentageFeatureSet(GSM4735375, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBM3","ALAS2")
HB_m <- match(HB.genes, rownames(GSM4735375@assays$RNA)) 
HB.genes <- rownames(GSM4735375@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM4735375[["percent.HB"]]<-PercentageFeatureSet(GSM4735375, features=HB.genes) 
col.num <- length(levels(GSM4735375@active.ident))

violin <- VlnPlot(GSM4735375,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#可视化
violin
ggsave("GSM4735375_vlnplot_qc.pdf", plot = violin, width = 8, height = 6) 

#质控
GSM4735375 <- subset(GSM4735375, subset =nFeature_RNA <5000 & nCount_RNA <20000 & percent.mt < 15& percent.HB < 1)

#在meta.data中添加group这一项
GSM4735375@meta.data$group <- "Tumor"
GSM4735375<- RenameCells(GSM4735375, add.cell.id = "GSM4735375")
GSM4735375@meta.data$orig.ident <- "GSM4735375"
GSM4735375@meta.data$stage <- 3
GSM4735375@meta.data$age <- 65
GSM4735375@meta.data$gender <- 2
GSM4735375@meta.data$sample <- "T14"



###开始harmony整合
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/Identify/cluster/")
library(harmony)
scRNA_harmony <- merge(GSM6415686, y=c(GSM6415687,GSM6415688,GSM6415689,GSM6415690,GSM6415691,GSM6415692,GSM6415693,GSM6415694,
                                       GSM4735364,GSM4735365,GSM4735366,GSM4735367,GSM4735368,
                                       GSM4735369,GSM4735370,GSM4735371,GSM4735372,GSM4735373,GSM4735374,GSM4735375))
scRNA_harmony
saveRDS(scRNA_harmony, "scRNA_merge.rds")
scRNA_harmony <- NormalizeData(scRNA_harmony,normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData() %>% RunPCA(pc.genes=VariableFeatures(object = scRNA_harmony,verbose = FALSE))

pdf(file="PCA.pdf",width=6.5,height=6)
ElbowPlot(scRNA_harmony,ndims=50)
dev.off()
saveRDS(scRNA_harmony, "scRNA_pca.rds")


####使用Harmony去批次###########################################################################################
#setwd("/home/stu-xt/ccRcc_24.9.2/harmony_22/xian_10/")
#scRNA_harmony <-readRDS("scRNA_pca.rds")
getOption("warn")
options(warn=1)
getOption("warn")
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident", plot_convergence = TRUE)
scRNA_harmony
saveRDS(scRNA_harmony, "harmony_integrate.rds")



############################################harmony循环的整个代码################################################
#rm(list=ls())
#setwd("/home/stu-xt/ccRcc_24.9.2/harmony_51/big_cluster/")
#scRNA_harmony <-readRDS("harmony_integrate.rds")
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/xunhuan/")
for (i in c(19,22,25,30,35,40,45,50)){
  print(i)
  pc.num = 1:i
  print(pc.num)
  ##细胞聚类
  scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = pc.num) 
  scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.5)
  ##非线性降维
  #tSNE
  scRNA_harmony = RunTSNE(scRNA_harmony, reduction = "harmony", dims = pc.num)
  plot1 = DimPlot(scRNA_harmony, pt.size = 0.1, label = TRUE,  reduction = "tsne") 
  #UMAP
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num)
  plot2 = DimPlot(scRNA_harmony, pt.size = 0.1, label = TRUE,  reduction = "umap") 
  #合并tSNE与UMAP
  plotc <- plot1+plot2+ plot_layout(guides = 'collect')
  tsne_umap_pdf <- paste(i, "_tSNE_UMAP.pdf", sep = "")
  ggsave(tsne_umap_pdf, plot = plotc, width = 10, height = 5)
}


#########选择不同的dim###########################################################################################
rm(list=ls())
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/Identify/cluster/")
scRNA_harmony <-readRDS("harmony_integrate.rds")
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/daqun/35_0.5/cluster/")
pcSelect =35
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:pcSelect) %>% FindClusters(resolution =0.5)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:pcSelect)                     
pdf(file="umap_35_0.5.pdf",width=6.5,height=6)
UMAPPlot(object = scRNA_harmony, pt.size = 0.1, label = TRUE,raster=FALSE)    
dev.off()

scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:pcSelect)         
pdf(file="TSNE_35_0.5.pdf",width=6.5,height=6)
TSNEPlot(object = scRNA_harmony, pt.size = 0.1, label = TRUE,raster=FALSE)    
dev.off()


pdf(file="sample1.pdf",width=6.5,height=6)
DimPlot(scRNA_harmony, reduction= "umap", label= F, group.by='orig.ident',repel = TRUE,raster=FALSE)
dev.off()



pdf(file="group1.pdf",width=6.5,height=6)
DimPlot(scRNA_harmony, reduction= "umap", label= F, group.by='group',repel = TRUE,raster=FALSE)
dev.off()

pdf(file="sample2.pdf",width=6.5,height=6)
DimPlot(scRNA_harmony, reduction= "tsne", label= F, group.by='orig.ident',repel = TRUE,raster=FALSE)
dev.off()



pdf(file="group2.pdf",width=6.5,height=6)
DimPlot(scRNA_harmony, reduction= "tsne", label= F, group.by='group',repel = TRUE,raster=FALSE)
dev.off()

saveRDS(scRNA_harmony, "scRNA_harmony.rds")



####展示marker
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/daqun/35_0.5/zhushi/")

#小提琴图
gene <- c("KLRD1","GNLY","KLRF1","NKG7","GNLY","NCAM1","KLRD1")
pdf(file="NK_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()


gene <- c("CD3D","CD3E","CD2","CD3G")
pdf(file="T_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()


gene <- c("CD163","CD68")
pdf(file="Macrophage_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <- c("CLEC9A","CPNE3","CD1C","FCER1A","CCR7","BIRC3")
pdf(file="DC_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <- c("CSF1R","CSF3R","CD68","MARCO",'FCGR3A',"FCN1")
pdf(file="myeloid _vio.pdf",width=20,height=10)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <- c("MZB1","CD79A","CD79B")
pdf(file="B_PLA_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()


gene <- c("EPCAM","KRT19","CLDN4","KRT18")
pdf(file="EPI_vio.pdf",width=20,height=10)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <- c('KIT',"CPA3","TPSAB1","GATA2")
pdf(file="MAST_vio.pdf",width=20,height=10)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()


gene <- c("MKI67","STMN1","PCNA")
pdf(file="zengzhi_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <- c('KDR',"VWF","CD34","FLT1","RAMP2")
pdf(file="endo_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <- c('LUM',"DCN","COL1A1","ACTA2","SPARC","COL1A2","THY1")
pdf(file="FIB_vio.pdf",width=20,height=10)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()



gene <- c("CD163","CD68","FCN1","APOBEC3A","THBS1","HLA-DQB2","HLA-DPB1","BIRC3")
pdf(file="Myeloid1_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <- c('C1QA',	'C1QB',	'C1QC',"FCN1","CD1C","LYN")
pdf(file="C1QA_vio.pdf",width=20,height=5)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()




gene <- c("C1QA","C1QB","C1QC","SELENOP","RNASE1","DAB2","LGMN","PLTP","MAF","SLCO2B1")
pdf(file="Mac_vio.pdf",width=50,height=20)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()


gene <- c("VCAN","FCN1","CD300E","S100A12","EREG","APOBEC3A","STXBP2","ASGR1","CCR2","NRG1")
pdf(file="Mono_vio.pdf",width=50,height=20)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <- c("FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2" )
pdf(file="neutrophils_vio.pdf",width=50,height=20)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()


gene <-  c("MT1G","MT1H","ALDOB","GPX3","MIOX")
pdf(file="n_epi_vio.pdf",width=50,height=20)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()

gene <-  c("CRYAB","CD24","NDUFA4L2","NNMT")
pdf(file="T_epi_vio.pdf",width=50,height=20)
VlnPlot(object =scRNA_harmony, features = gene, pt.size = 0, slot = "data", assay = "RNA",raster=FALSE)
dev.off()


####去除双细胞；
#19群：成纤维和内皮；26：巨噬和内皮；27：上皮和成纤维；28：成纤维和巨噬；29：巨噬和上皮；30：成纤维和上皮
rm(list=ls())
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/daqun/35_0.5/cluster/")
scRNA_harmony<- readRDS("scRNA_harmony.rds")
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/daqun/35_0.5/cluster/quchu/")
scRNA_harmony<- subset(scRNA_harmony, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25))
pdf(file="umap_35_0.5.pdf",width=6.5,height=6)
UMAPPlot(object = scRNA_harmony, pt.size = 0.1, label = TRUE,raster=FALSE)    
dev.off()

       
pdf(file="TSNE_35_0.5.pdf",width=6.5,height=6)
TSNEPlot(object = scRNA_harmony, pt.size = 0.1, label = TRUE,raster=FALSE)    
dev.off()


pdf(file="sample1.pdf",width=6.5,height=6)
DimPlot(scRNA_harmony, reduction= "umap", label= F, group.by='orig.ident',repel = TRUE,raster=FALSE)
dev.off()



pdf(file="group1.pdf",width=6.5,height=6)
DimPlot(scRNA_harmony, reduction= "umap", label= F, group.by='group',repel = TRUE,raster=FALSE)
dev.off()

pdf(file="sample2.pdf",width=6.5,height=6)
DimPlot(scRNA_harmony, reduction= "tsne", label= F, group.by='orig.ident',repel = TRUE,raster=FALSE)
dev.off()



pdf(file="group2.pdf",width=6.5,height=6)
DimPlot(scRNA_harmony, reduction= "tsne", label= F, group.by='group',repel = TRUE,raster=FALSE)
dev.off()

saveRDS(scRNA_harmony, "scRNA_harmony.rds")
#######每个群的marker
#top10marker
logFCfilter=0.5
adjPvalFilter=0.05
m.markers <- FindAllMarkers(object =scRNA_harmony,
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.markers=m.markers[(abs(as.numeric(as.vector(m.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(m.markers$p_val_adj))<adjPvalFilter),]
top10 <- m.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(top10,file="top20_cluster.txt",sep="\t",row.names=F,quote=F)



####大群手动注释
rm(list=ls())
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/daqun/35_0.5/cluster/quchu/")
scRNA_harmony<- readRDS("scRNA_harmony.rds")
setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/PSME2/daqun/Rename/")
new.cluster.ids <- c("Epithelial","Endothelial","T","NK","Fibroblast","Macrophages","T","Epithelial",
                  "Monocytes","Epithelial","Endothelial","Epithelial","Epithelial","Endothelial", "B","Epithelial",
                    "Proliferation","Fibroblast","Mast","Epithelial","Epithelial","Endothelial","Epithelial",
                     "Epithelial","Epithelial")
names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony<- RenameIdents(scRNA_harmony, new.cluster.ids)
pdf(file="umap_rename.pdf",width=7,height=5.5)
DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, pt.size = 0.15,raster=FALSE)
dev.off()

pdf(file="tsne_rename.pdf",width=7,height=5.5)
DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE, pt.size = 0.15,raster=FALSE)
dev.off()

saveRDS(scRNA_harmony, "harmony_Rename.rds")



setwd("/home/stu-xt/ccRcc_24.9.2/harmony_21/daqun/35_0.5/Rename/")
scRNA_harmony<- readRDS("harmony_Rename.rds")



pdf(file="geneMarker.pdf",width=20,height=8)
DotPlot(scRNA_harmony,features = c("SPP1","APOE","TREM2","FABP5","C1QA","C1QB","C1QC","CD9","FCGRT","MRC1","CD163","MAFB",
                         "CD14","VCAN","FCGR3A","CDKN1C","LILRB2","S100A8","S100A9","CSF3R","PLAC8"
                         ))+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")+scale_color_gradientn(colours = rev(c("red3","grey")))
dev.off()

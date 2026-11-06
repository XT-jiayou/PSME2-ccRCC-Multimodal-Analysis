rm(list = ls())
# ============================================================
# 1. 设置工作目录并读取数据
# ============================================================
setwd("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/Epi/cellchat")
library(Matrix)
library(CellChat)
library(qs)
library(Seurat)
library(patchwork)
library(BiocParallel)
library(ggplot2)
register(MulticoreParam(workers = 8, progressbar = TRUE))
scRNA <- readRDS(
  "/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/Epi/cellchat/scRNA_tumor_PSME2_annotated.rds")
# ============================================================
# 2. 检查细胞注释并提取表达矩阵
# ============================================================
Idents(scRNA) <- "cellchat_annotation"
table(Idents(scRNA))
phe <- scRNA@meta.data
table(phe$cellchat_annotation)
p1 <- DimPlot(
  scRNA,
  reduction = "umap",
  pt.size = 0.1,
  group.by = "cellchat_annotation",
  label = TRUE,
  repel = TRUE)
ggsave( "CellChat_UMAP.pdf", plot = p1, width = 12, height = 9)
data.input <- GetAssayData(
  scRNA,
  assay = "RNA",
  layer = "data")
dim(data.input)
# ============================================================
# 3. 创建CellChat对象
# ============================================================
meta <- scRNA@meta.data[, c("orig.ident", "cellchat_annotation")]
colnames(meta) <- c("group", "labels")
table(meta$labels)
stopifnot(identical(rownames(meta), colnames(data.input)))
cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "labels")
levels(cellchat@idents)
table(cellchat@idents)
# ============================================================
# 4. 设置配体-受体相互作用数据库
# ============================================================
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
dim(CellChatDB.use$interaction)
head(CellChatDB.use$interaction$ligand)
cellchat@DB <- CellChatDB.use
# ============================================================
# 5. 预处理细胞-细胞通讯分析的表达数据
# ============================================================
cellchat <- subsetData(cellchat)
dim(cellchat@data.signaling)
future::plan("sequential")
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
saveRDS(cellchat, "cellchat_step1.rds")
# ============================================================
# 6. 计算细胞间通讯概率
# ============================================================
cellchat <- computeCommunProb(
  cellchat,
  type = "triMean",
  raw.use = FALSE)
# 过滤细胞数量少于10的通讯
cellchat <- filterCommunication(
  cellchat,
  min.cells = 10)
# 保存结果
saveRDS(cellchat, "cellchat_step2.rds")

###7.在信号通路水平上推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)
# 提取全部11个细胞群的通讯结果
df.net <- subsetCommunication(cellchat)
# 保存CellChat对象
saveRDS(cellchat, "cellchat_step3.rds")
save(df.net, file = "Tumor_PSME2_df.net.Rdata")
write.csv( df.net, "Tumor_PSME2_df.net.csv", row.names = FALSE)

#8.可视化流程
# 聚合全部11个细胞群的通讯网络
cellchat <- aggregateNet(cellchat)
# 获取各细胞群的细胞数量
groupSize <- as.numeric(table(cellchat@idents))
# 保存通讯网络图
pdf("Tumor_PSME2_CellChat_circle.pdf", width = 16, height = 8)
par(mfrow = c(1, 2), xpd = TRUE)
# 通讯数量
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions")
# 通讯强度
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength")
dev.off()

# 保存完整CellChat对象
saveRDS(cellchat, "cellchat_step4.rds")

setwd("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/Epi/cellchat")
cellchat <- readRDS("cellchat_step4.rds")
setwd("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/Epi/cellchat/1/")
##9.heatmap展示互作数据
pdf("Tumor_PSME2_CellChat_count_heatmap.pdf", width = 12, height = 10)
pheatmap::pheatmap(
  cellchat@net$count,
  border_color = "black",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize = 10,
  display_numbers = TRUE,
  number_color = "black",
  number_format = "%.0f")
dev.off()

###10.贝克图
# ============================================================
# 10. 绘制11个细胞群分别作为发送者的通讯网络
# ============================================================

# 设置11个细胞群的显示顺序
celltype_order <- c( "PSME2- Epi", "PSME2+ Epi","Macrophage", "Monocyte", "T",
  "NK", "B","Endothelial", "Fibroblast", "Mast", "Proliferation")
# 提取通讯强度矩阵并排序
mat <- cellchat@net$weight
mat <- mat[celltype_order, celltype_order, drop = FALSE]
# 设置细胞群颜色
color.use <- rainbow(nrow(mat))
names(color.use) <- rownames(mat)
# 设置细胞数量
groupSize <- as.numeric(table(cellchat@idents)[celltype_order])
# 保存PDF
pdf("Tumor_PSME2_CellChat_each_sender.pdf", width = 16, height = 12)
par(mfrow = c(3, 4), xpd = TRUE, mar = c(1, 1, 2, 1))
for (i in seq_len(nrow(mat))) {
  mat2 <- matrix(
    0,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat) )
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    arrow.size = 0.05,
    arrow.width = 1,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i],
    color.use = color.use )}
dev.off()

#11.层次结构图
# ============================================================
library(CellChat)
setwd("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/Epi/cellchat")
cellchat <- readRDS("cellchat_step4.rds")
# 查看信号通路及细胞群
cellchat@netP$pathways
levels(cellchat@idents)
# 选择信号通路
pathways.show <- "MIF"
# ============================================================
# 1. Hierarchy plot
# ============================================================
# 左侧设置PSME2- Epi、PSME2+ Epi
vertex.receiver <- match(
  c("PSME2- Epi", "PSME2+ Epi"),
  levels(cellchat@idents))
pdf("MIF_hierarchy.pdf", width = 12, height = 8)
netVisual_aggregate(
  cellchat,
  signaling = pathways.show,
  vertex.receiver = vertex.receiver,
  layout = "hierarchy")
dev.off()
# ============================================================
# 2. Circle plot
# ============================================================

pdf("MIF_circle.pdf", width = 10, height = 10)

netVisual_aggregate(
  cellchat,
  signaling = pathways.show,
  layout = "circle"
)

dev.off()

# ============================================================
# 3. Chord diagram
# ============================================================

pdf("MIF_chord.pdf", width = 12, height = 12)

netVisual_aggregate(
  cellchat,
  signaling = pathways.show,
  layout = "chord"
)

dev.off()

# ============================================================
# 4. 分组Chord diagram
# ============================================================

group.cellType <- c(
  "PSME2- Epi" = "Epithelial",
  "PSME2+ Epi" = "Epithelial",
  "Macrophage" = "Myeloid",
  "Monocyte" = "Myeloid",
  "T" = "Lymphoid",
  "NK" = "Lymphoid",
  "B" = "Lymphoid",
  "Endothelial" = "Stromal",
  "Fibroblast" = "Stromal",
  "Mast" = "Myeloid",
  "Proliferation" = "Proliferation"
)

group.cellType <- group.cellType[levels(cellchat@idents)]

pdf("MIF_chord_group.pdf", width = 12, height = 12)

netVisual_chord_cell(
  cellchat,
  signaling = pathways.show,
  group = group.cellType,
  title.name = paste0(pathways.show, " signaling network")
)

dev.off()

# ============================================================
# 5. Heatmap
# ============================================================

pdf("MIF_heatmap.pdf", width = 12, height = 10)

netVisual_heatmap(
  cellchat,
  signaling = pathways.show,
  color.heatmap = "Reds"
)

dev.off()


# ============================================================
# 12. MIF通路配体-受体贡献分析及可视化
# ============================================================
library(CellChat)
setwd("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/Epi/cellchat")
#cellchat <- readRDS("cellchat_step4.rds")
pathways.show <- "MIF"
# 1. 计算配体-受体对MIF通路的贡献
pdf("MIF_LR_contribution.pdf", width = 8, height = 6)
netAnalysis_contribution(
  cellchat,
  signaling = pathways.show)
dev.off()
# 2. 提取MIF通路的配体-受体对
pairLR <- extractEnrichedLR(
  cellchat,
  signaling = pathways.show,
  geneLR.return = FALSE)
print(pairLR)

# 选择MIF-CD74/CXCR4
LR.show <- pairLR[
  pairLR$interaction_name == "MIF_CD74_CXCR4",
  ,
  drop = FALSE
]

# 如果数据库中没有该组合，则使用第一个配体-受体对
if (nrow(LR.show) == 0) {
  LR.show <- pairLR[1, , drop = FALSE]
}

print(LR.show)

# 3. Hierarchy plot
vertex.receiver <- match(
  c("PSME2- Epi", "PSME2+ Epi"),
  levels(cellchat@idents)
)

pdf("MIF_LR_hierarchy.pdf", width = 12, height = 8)

netVisual_individual(
  cellchat,
  signaling = pathways.show,
  pairLR.use = LR.show,
  vertex.receiver = vertex.receiver,
  layout = "hierarchy"
)

dev.off()

# 4. Circle plot
pdf("MIF_LR_circle.pdf", width = 10, height = 10)

netVisual_individual(
  cellchat,
  signaling = pathways.show,
  pairLR.use = LR.show,
  layout = "circle"
)

dev.off()


# ============================================================
# 13. PSME2+/- Epi → Macrophage 信号通路信息流比较
# ============================================================

rm(list = ls())

library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/Epi/cellchat")

cellchat <- readRDS("cellchat_step4.rds")

# 1. 提取通路水平通讯概率
prob <- cellchat@netP$prob

groups <- c("PSME2- Epi", "PSME2+ Epi")

stopifnot( all(groups %in% dimnames(prob)[[1]]), "Macrophage" %in% dimnames(prob)[[2]] )

# 2. 提取两组上皮细胞到巨噬细胞的通讯强度
df <- data.frame(
  Pathway = dimnames(prob)[[3]],
  Negative = as.numeric(prob["PSME2- Epi", "Macrophage", ]),
  Positive = as.numeric(prob["PSME2+ Epi", "Macrophage", ])
)

# 去除两组均为0的通路
df <- df %>%
  filter(
    is.finite(Negative),
    is.finite(Positive),
    Negative + Positive > 0
  ) %>%
  mutate(
    Total = Negative + Positive,
    Negative_relative = Negative / Total,
    Positive_relative = Positive / Total
  ) %>%
  arrange(desc(Total))

# 保存原始比较结果
write.csv( df, "PSME2_Epi_to_Macrophage_pathway_comparison.csv", row.names = FALSE )

# 3. 整理绘图数据
pathway_order <- rev(df$Pathway)

df.absolute <- df %>%
  select(Pathway, Negative, Positive) %>%
  pivot_longer(
    cols = c(Negative, Positive),
    names_to = "Group",
    values_to = "Information_flow"
  )

df.relative <- df %>%
  select(Pathway, Negative_relative, Positive_relative) %>%
  pivot_longer(
    cols = c(Negative_relative, Positive_relative),
    names_to = "Group",
    values_to = "Relative_flow"
  ) %>%
  mutate( Group = ifelse( Group == "Negative_relative", "Negative", "Positive" ) )

df.absolute$Pathway <- factor(
  df.absolute$Pathway,
  levels = pathway_order
)

df.relative$Pathway <- factor(
  df.relative$Pathway,
  levels = pathway_order
)

# 4. 设置颜色
group_colors <- c(
  "Negative" = "#F8766D",
  "Positive" = "#00BFC4"
)

# 5. 左图：相对信息流
p1 <- ggplot(
  df.relative,
  aes(x = Relative_flow, y = Pathway, fill = Group)
) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous( limits = c(0, 1), breaks = seq(0, 1, 0.25) ) +
  labs( x = "Relative information flow", y = NULL, fill = NULL ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(size = 8)
  )

# 6. 右图：通讯强度
p2 <- ggplot(
  df.absolute,
  aes(x = Information_flow, y = Pathway, fill = Group)
) +
  geom_col( position = position_dodge(width = 0.8), width = 0.7 ) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous( expand = expansion(mult = c(0, 0.08)) ) +
  labs( x = "Information flow", y = NULL, fill = NULL ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(size = 8)
  )

# 7. 合并图片
p <- p1 + p2 +
  plot_layout(widths = c(1, 1))

# 8. 保存PDF
ggsave(
  "PSME2_Epi_Macrophage_Information_flow.pdf",
  plot = p,
  width = 14,
  height = max(7, nrow(df) * 0.23 + 2)
)

# 9. 查看MIF、SPP1和MHC-II的具体数值
df %>%
  filter(Pathway %in% c("MIF", "SPP1", "MHC-II")) %>%
  print()


# ============================================================
# PSME2+/- Epi → Macrophage
# Top 5 signaling pathways information flow
# ============================================================

rm(list = ls())

library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd("/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/Epi/cellchat")

# 1. 读取CellChat对象
cellchat <- readRDS("cellchat_step4.rds")

# 2. 提取信号通路通讯概率
prob <- cellchat@netP$prob

df.all <- data.frame(
  Pathway = dimnames(prob)[[3]],
  Negative = as.numeric(prob["PSME2- Epi", "Macrophage", ]),
  Positive = as.numeric(prob["PSME2+ Epi", "Macrophage", ])
)

# 3. 计算累计通讯概率并筛选前5条
df.all <- df.all %>%
  filter(
    is.finite(Negative),
    is.finite(Positive),
    Negative + Positive > 0
  ) %>%
  mutate(
    Total = Negative + Positive,
    Negative_relative = Negative / Total,
    Positive_relative = Positive / Total
  ) %>%
  arrange(desc(Total))

# 保存全部通路结果
write.csv( df.all, "PSME2_Epi_Macrophage_All_pathways.csv", row.names = FALSE )

# 只取前5条
df <- df.all %>%
  slice_head(n = 5)

print(df)

# 4. 整理绘图数据
pathway_order <- rev(df$Pathway)

df.relative <- df %>%
  select(Pathway, Negative_relative, Positive_relative) %>%
  pivot_longer(
    cols = c(Negative_relative, Positive_relative),
    names_to = "Group",
    values_to = "Relative_flow"
  ) %>%
  mutate(
    Group = ifelse(
      Group == "Negative_relative",
      "PSME2- Epi",
      "PSME2+ Epi"
    ),
    Pathway = factor(Pathway, levels = pathway_order)
  )

df.absolute <- df %>%
  select(Pathway, Negative, Positive) %>%
  pivot_longer(
    cols = c(Negative, Positive),
    names_to = "Group",
    values_to = "Information_flow"
  ) %>%
  mutate(
    Group = ifelse(
      Group == "Negative",
      "PSME2- Epi",
      "PSME2+ Epi"
    ),
    Pathway = factor(Pathway, levels = pathway_order)
  )

# 5. 设置颜色
group_colors <- c(
  "PSME2- Epi" = "#F8766D",
  "PSME2+ Epi" = "#00BFC4"
)

# 6. 左图：相对信息流
p1 <- ggplot(
  df.relative,
  aes( x = Relative_flow, y = Pathway, fill = Group )
) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous( limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0) ) +
  labs( x = "Relative information flow", y = NULL, fill = NULL ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom"
  )

# 7. 右图：通讯强度
p2 <- ggplot(
  df.absolute,
  aes( x = Information_flow, y = Pathway, fill = Group )
) +
  geom_col( position = position_dodge(width = 0.8), width = 0.7 ) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous( expand = expansion(mult = c(0, 0.1)) ) +
  labs( x = "Information flow", y = NULL, fill = NULL ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom"
  )

# 8. 合并两张图
p <- p1 + p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# 9. 保存PDF
ggsave( "PSME2_Epi_Macrophage_Top5_Information_flow.pdf", plot = p, width = 12, height = 5.5 )

# 10. 保存前5条通路数据
write.csv( df, "PSME2_Epi_Macrophage_Top5_pathways.csv", row.names = FALSE )


####最终结果出图代码
###排序图
## ============================================================
## CellChat Macrophage-associated Signaling Pathways
##
## 自动计算全部显著信号通路
## 根据累计通讯概率选取Top 5
##
## 四个通讯方向：
## PSME2+ Epi -> Macrophage
## PSME2- Epi -> Macrophage
## Macrophage -> PSME2+ Epi
## Macrophage -> PSME2- Epi
##
## P < 0.05
##
## Output: PDF + PNG + CSV
## Folder: cellchat/final/
## ============================================================

library(CellChat)
library(ggplot2)
library(dplyr)
library(scales)

# ============================================================
# 1. 设置输入和输出路径
# ============================================================

indir <- paste0(
  "/home/stu-xt/ccRcc_24.9.2/ZC/Epi/",
  "scissor/WSI/Epi/cellchat"
)

outdir <- file.path(indir, "final")

dir.create( outdir, recursive = TRUE, showWarnings = FALSE )

# ============================================================
# 2. 读取CellChat对象
# ============================================================

cellchat <- readRDS(
  file.path(indir, "cellchat_step4.rds")
)

cat("\nCellChat对象类型：\n")
print(class(cellchat))

cat("\n细胞群名称：\n")
print(levels(cellchat@idents))

# ============================================================
# 3. 检查细胞群
# ============================================================

required_groups <- c(
  "PSME2+ Epi",
  "PSME2- Epi",
  "Macrophage"
)

actual_groups <- levels(cellchat@idents)

if (!all(required_groups %in% actual_groups)) {
  stop("CellChat细胞群名称不匹配。")
}

# ============================================================
# 4. 提取全部通讯数据
# ============================================================

df_all <- subsetCommunication(
  cellchat,
  thresh = 1
)

cat( "\n全部通讯记录：", nrow(df_all), "\n" )

# ============================================================
# 5. 提取四个通讯方向
# ============================================================

df_core <- df_all %>%
  filter(
    (
      source %in% c(
        "PSME2+ Epi",
        "PSME2- Epi"
      ) &
      target == "Macrophage"
    ) |
    (
      source == "Macrophage" &
      target %in% c(
        "PSME2+ Epi",
        "PSME2- Epi"
      )
    )
  ) %>%
  mutate( Group = paste0( source, " -> ", target ) )

cat( "\n四个方向通讯记录：", nrow(df_core), "\n" )

# ============================================================
# 6. 筛选显著通讯
# ============================================================

df_sig <- df_core %>%
  filter(
    !is.na(pathway_name),
    !is.na(prob),
    !is.na(pval),
    pval < 0.05
  )

cat( "\n显著通讯记录：", nrow(df_sig), "\n" )

if (nrow(df_sig) == 0) {
  stop("没有找到P < 0.05的显著通讯。")
}

# ============================================================
# 7. 计算全部通路累计通讯概率
# ============================================================

# 每条通路：
# 四个通讯方向中所有显著LR通讯概率求和
#
# 不预先指定通路名称

pathway_summary_all <- df_sig %>%
  group_by(pathway_name) %>%
  summarise(
    Cumulative_Probability = sum(
      prob,
      na.rm = TRUE
    ),

    Interaction_Count = n(),

    .groups = "drop"
  ) %>%
  arrange( desc(Cumulative_Probability) )

# ============================================================
# 8. 仅选取累计通讯概率最高的前5条通路
# ============================================================

# 本次唯一的分析筛选修改

pathway_summary <- pathway_summary_all %>%
  slice_head(n = 5)

n_pathways <- nrow(pathway_summary)

cat("\n====================================\n")

cat( "全部显著信号通路数量：", nrow(pathway_summary_all), "\n" )

cat( "实际展示Top通路数量：", n_pathways, "\n" )

cat("\nTop 5信号通路：\n")

print( as.data.frame(pathway_summary), row.names = FALSE )

cat("====================================\n")

# ============================================================
# 9. 设置纵坐标顺序
# ============================================================

# ggplot2的Y轴从下向上排列
# 使用rev()使最高值位于最上方

pathway_order <- as.character(
  pathway_summary$pathway_name
)

pathway_summary$pathway_name <- factor(
  pathway_summary$pathway_name,
  levels = rev(pathway_order)
)

# ============================================================
# 10. 设置数值标签
# ============================================================

pathway_summary <- pathway_summary %>%
  mutate( Label = sprintf( "%.3f", Cumulative_Probability ) )

# ============================================================
# 11. 自动设置横坐标范围
# ============================================================

max_prob <- max(
  pathway_summary$Cumulative_Probability,
  na.rm = TRUE
)

x_max <- max_prob * 1.15

# ============================================================
# 12. 设置图片尺寸
# ============================================================

plot_width <- 7

plot_height <- max(
  4,
  n_pathways * 0.42 + 1.8
)

cat( "\n图片尺寸：", plot_width, "x", plot_height, "英寸\n" )

# ============================================================
# 13. 绘制横向柱形图
# ============================================================

p <- ggplot(
  pathway_summary,
  aes( x = Cumulative_Probability, y = pathway_name )
) +

  # 蓝色柱形
  geom_col(
    width = 0.56,
    fill = "#6BAED6"
  ) +

  # 柱形末端数值
  geom_text( aes( label = Label ), hjust = -0.08, size = 3.2, color = "black" ) +

  # 横坐标
  scale_x_continuous(
    limits = c(0, x_max),
    breaks = pretty(
      c(0, max_prob),
      n = 6
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +

  # 纵坐标
  scale_y_discrete( limits = rev(pathway_order), drop = FALSE ) +

  # 标题和坐标轴
  labs(
    title = paste0(
      "Macrophage-associated signaling pathways ",
      "between PSME2 states"
    ),

    subtitle = "Ranked by cumulative interaction probability",

    x = "Cumulative Probability",

    y = "Pathway Name"
  ) +

  theme_classic( base_size = 11 ) +

  theme(

    # 标题
    plot.title = element_text(
      size = 11,
      face = "bold",
      color = "black",
      hjust = 0
    ),

    # 副标题
    plot.subtitle = element_text(
      size = 10,
      color = "grey40",
      hjust = 0
    ),

    # X轴标题
    axis.title.x = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),

    # Y轴标题
    axis.title.y = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),

    # X轴刻度
    axis.text.x = element_text(
      size = 9,
      color = "black"
    ),

    # Y轴通路名称
    axis.text.y = element_text(
      size = 10,
      face = "bold",
      color = "#B22222"
    ),

    # 坐标轴线
    axis.line = element_line(
      linewidth = 0.5,
      color = "black"
    ),

    axis.ticks = element_line(
      color = "black"
    ),

    # 纵向网格线
    panel.grid.major.x = element_line(
      color = "#D9D9D9",
      linewidth = 0.4
    ),

    panel.grid.major.y = element_blank(),

    panel.grid.minor = element_blank(),

    plot.margin = margin(
      10, 15, 10, 10
    )
  )

# ============================================================
# 14. 显示图片
# ============================================================

print(p)

# ============================================================
# 15. 保存PDF
# ============================================================

ggsave(
  filename = file.path(
    outdir,
    "CellChat_Top5_Pathways_Cumulative_Probability.pdf"
  ),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

# ============================================================
# 16. 保存PNG
# ============================================================

ggsave(
  filename = file.path(
    outdir,
    "CellChat_Top5_Pathways_Cumulative_Probability.png"
  ),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

# ============================================================
# 17. 保存Top 5通路计算结果
# ============================================================

write.csv(
  as.data.frame(pathway_summary),
  file = file.path(
    outdir,
    "CellChat_Top5_Pathways_Cumulative_Probability.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 18. 保存全部通路计算结果
# ============================================================

# 保留全部通路数据，方便核对排序

write.csv(
  as.data.frame(pathway_summary_all),
  file = file.path(
    outdir,
    "CellChat_All_Pathways_Cumulative_Probability.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 19. 保存显著通讯原始数据
# ============================================================

write.csv(
  as.data.frame(df_sig),
  file = file.path(
    outdir,
    "CellChat_All_Pathways_Significant_Interactions.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 20. 最终检查
# ============================================================

cat("\n====================================\n")

cat("CellChat Top 5信号通路排序图完成！\n")

cat( "全部显著信号通路数量：", nrow(pathway_summary_all), "\n" )

cat( "展示信号通路数量：", n_pathways, "\n" )

cat( "显著通讯记录数量：", nrow(df_sig), "\n" )

cat( "累计通讯概率最高的通路：", pathway_order[1], "\n" )

cat( "最高累计通讯概率：", max_prob, "\n" )

cat( "输出文件夹：", outdir, "\n" )

cat("====================================\n")

####总的图

## ============================================================
## CellChat Overall Communication Network
##
## Left: Number of interactions
## Right: Interaction weights/strength
##
## Input: cellchat_step4.rds
## Output: PDF + PNG
## Folder: cellchat/final/
## ============================================================

library(CellChat)
library(ggplot2)

# ============================================================
# 1. 设置输入和输出路径
# ============================================================

indir <- paste0(
  "/home/stu-xt/ccRcc_24.9.2/ZC/Epi/",
  "scissor/WSI/Epi/cellchat"
)

outdir <- file.path(indir, "final")

dir.create( outdir, recursive = TRUE, showWarnings = FALSE )

# ============================================================
# 2. 读取CellChat对象
# ============================================================

cellchat <- readRDS(
  file.path(indir, "cellchat_step4.rds")
)

cat("\nCellChat对象类型：\n")
print(class(cellchat))

cat("\n细胞群名称：\n")
print(levels(cellchat@idents))

# ============================================================
# 3. 检查总体通讯矩阵
# ============================================================

if (
  is.null(cellchat@net$count) ||
  is.null(cellchat@net$weight)
) {
  stop( "CellChat对象缺少count或weight矩阵，请检查分析步骤。" )
}

mat_count <- cellchat@net$count
mat_weight <- cellchat@net$weight

cat("\n通讯数量矩阵维度：\n")
print(dim(mat_count))

cat("\n通讯强度矩阵维度：\n")
print(dim(mat_weight))

# ============================================================
# 4. 设置细胞群大小
# ============================================================

group_size <- as.numeric(
  table(cellchat@idents)[
    rownames(mat_count)
  ]
)

if (anyNA(group_size)) {
  stop("细胞群数量与通讯矩阵名称不匹配。")
}

# ============================================================
# 5. 设置细胞群颜色
# ============================================================

group_names <- rownames(mat_count)

group_colors <- grDevices::hcl.colors(
  length(group_names),
  palette = "Dynamic"
)

names(group_colors) <- group_names

# ============================================================
# 6. 定义绘图函数
# ============================================================

plot_circle_network <- function() {

  par(
    mfrow = c(1, 2),
    mar = c(2, 2, 3, 2),
    xpd = TRUE
  )

  # ----------------------------------------------------------
  # 左图：通讯数量
  # ----------------------------------------------------------

  netVisual_circle(
    mat_count,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = group_colors,
    vertex.label.cex = 0.9,
    vertex.size.max = 18,
    edge.width.max = 8,
    title.name = "Number of interactions"
  )

  # ----------------------------------------------------------
  # 右图：通讯强度
  # ----------------------------------------------------------

  netVisual_circle(
    mat_weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = group_colors,
    vertex.label.cex = 0.9,
    vertex.size.max = 18,
    edge.width.max = 8,
    title.name = "Interaction weights/strength"
  )
}

# ============================================================
# 7. 在RStudio中显示
# ============================================================

plot_circle_network()

# ============================================================
# 8. 保存PDF
# ============================================================

pdf(
  file = file.path(
    outdir,
    "CellChat_Overall_Communication_Network.pdf"
  ),
  width = 12,
  height = 6.5,
  useDingbats = FALSE
)

plot_circle_network()

dev.off()

# ============================================================
# 9. 保存PNG
# ============================================================

png(
  filename = file.path(
    outdir,
    "CellChat_Overall_Communication_Network.png"
  ),
  width = 12,
  height = 6.5,
  units = "in",
  res = 300,
  bg = "white"
)

plot_circle_network()

dev.off()

# ============================================================
# 10. 保存通讯矩阵
# ============================================================

write.csv( mat_count, file = file.path( outdir, "CellChat_Interaction_Count_Matrix.csv" ) )

write.csv( mat_weight, file = file.path( outdir, "CellChat_Interaction_Weight_Matrix.csv" ) )

# ============================================================
# 11. 输出完成信息
# ============================================================

cat("\n====================================\n")
cat("CellChat总体通讯网络图绘制完成！\n")

cat( "细胞群数量：", length(group_names), "\n" )

cat( "输出文件夹：", outdir, "\n" )

cat("====================================\n")


####气泡图

## ============================================================
## CellChat Ligand-Receptor Bubble Plot
##
## Six signaling pathways:
## IL6 / MIF / SPP1 / MHC-II / APP / CCL
##
## Four communication directions
##
## 本次修改：
## 1. MIF显示在最上面
## 2. SPP1紧接MIF
## 3. 其他配体-受体对保持字母排序
## 4. complete()后重新设置因子顺序
## 5. scale_y_discrete()明确锁定纵坐标顺序
##
## Color: Communication Probability
## Size: -log10(p-value)
##
## Output: PDF + PNG + CSV
## Folder: cellchat/final/
## ============================================================

library(CellChat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(grid)

# ============================================================
# 1. 设置输入和输出文件夹
# ============================================================

indir <- paste0(
  "/home/stu-xt/ccRcc_24.9.2/ZC/Epi/",
  "scissor/WSI/Epi/cellchat"
)

outdir <- file.path(indir, "final")

dir.create( outdir, recursive = TRUE, showWarnings = FALSE )

# ============================================================
# 2. 读取CellChat对象
# ============================================================

cellchat <- readRDS(
  file.path(indir, "cellchat_step4.rds")
)

cat("\nCellChat对象类型：\n")
print(class(cellchat))

cat("\n细胞群名称：\n")
print(levels(cellchat@idents))

# ============================================================
# 3. 设置四个通讯方向
# ============================================================

group_order <- c(
  "PSME2+ Epi -> Macrophage",
  "PSME2- Epi -> Macrophage",
  "Macrophage -> PSME2+ Epi",
  "Macrophage -> PSME2- Epi"
)

required_groups <- c(
  "PSME2+ Epi",
  "PSME2- Epi",
  "Macrophage"
)

actual_groups <- levels(cellchat@idents)

if (!all(required_groups %in% actual_groups)) {
  stop( "CellChat对象中的细胞群名称不匹配。" )
}

# ============================================================
# 4. 提取全部通讯结果
# ============================================================

df_all <- subsetCommunication(
  cellchat,
  thresh = 1
)

# ============================================================
# 5. 提取四个通讯方向
# ============================================================

df_core <- df_all %>%
  filter(
    (
      source %in% c(
        "PSME2+ Epi",
        "PSME2- Epi"
      ) &
      target == "Macrophage"
    ) |
    (
      source == "Macrophage" &
      target %in% c(
        "PSME2+ Epi",
        "PSME2- Epi"
      )
    )
  ) %>%
  mutate(
    Group = paste0(
      source,
      " -> ",
      target
    ),
    interaction_name = as.character(
      interaction_name
    ),
    interaction = as.character(
      interaction_name_2
    ),
    ligand = as.character(ligand),
    pathway_name = as.character(pathway_name)
  ) %>%
  filter(
    !is.na(interaction_name),
    !is.na(interaction),
    !is.na(prob),
    !is.na(pval)
  )

df_core$Group <- factor(
  df_core$Group,
  levels = group_order
)

cat( "\n四组全部通讯记录：", nrow(df_core), "\n" )

# ============================================================
# 6. 只筛选6条信号通路
# ============================================================

selected_pathways <- c(
  "IL6",
  "MIF",
  "SPP1",
  "MHC-II",
  "APP",
  "CCL"
)

df_core <- df_core %>%
  filter(
    pathway_name %in% selected_pathways
  )

cat( "\n六条通路全部通讯记录：", nrow(df_core), "\n" )

cat("\n六条通路实际匹配情况：\n")

print(
  as.data.frame(
    df_core %>%
      count(pathway_name)
  )
)

missing_pathways <- setdiff(
  selected_pathways,
  unique(df_core$pathway_name)
)

cat("\n未匹配的通路：\n")
print(missing_pathways)

# ============================================================
# 7. 提取显著通讯
# ============================================================

df_sig <- df_core %>%
  filter(
    pval < 0.05
  )

cat( "\n六条通路显著通讯记录：", nrow(df_sig), "\n" )

# ============================================================
# 8. 显著配体-受体对去重
# ============================================================

LR_sig <- df_sig %>%
  distinct(
    interaction_name,
    .keep_all = TRUE
  ) %>%
  select( interaction_name, interaction, ligand, pathway_name )

n_LR <- nrow(LR_sig)

cat( "\n去重后显著配体-受体对：", n_LR, "\n" )

if (n_LR == 0) {
  stop( "六条通路中没有找到显著配体-受体对。" )
}

# ============================================================
# 9. 保留所有显著配体-受体对
# ============================================================

selected_LR <- as.character(
  LR_sig$interaction_name
)

df <- df_core %>%
  filter(
    interaction_name %in% selected_LR
  )

# ============================================================
# 10. 计算-log10(p-value)
# ============================================================

df <- df %>%
  mutate(
    log10_p = pmin(
      -log10(
        pmax(pval, 1e-10)
      ),
      10
    ),
    significant = ifelse(
      pval < 0.05,
      "yes",
      "no"
    )
  )

# ============================================================
# 11. 建立唯一配体-受体标签
#     MIF第一，SPP1第二
#     其他保持字母排序
# ============================================================

LR_order_df <- LR_sig %>%
  mutate(
    interaction_name = as.character(
      interaction_name
    ),
    interaction = as.character(
      interaction
    ),
    ligand = as.character(
      ligand
    ),

    sort_group = case_when(
      ligand == "MIF"  ~ 1,
      ligand == "SPP1" ~ 2,
      TRUE            ~ 3
    )
  ) %>%
  arrange(
    sort_group,
    ligand,
    interaction
  ) %>%
  select( -sort_group )

LR_order_df$plot_id <- make.unique(
  LR_order_df$interaction,
  sep = "_"
)

# ============================================================
# 12. 使用match()精确匹配
# ============================================================

df$interaction_name <- as.character(
  df$interaction_name
)

df$plot_id <- LR_order_df$plot_id[
  match( df$interaction_name, LR_order_df$interaction_name )
]

n_missing <- sum(
  is.na(df$plot_id)
)

cat( "\n未匹配记录数量：", n_missing, "\n" )

if (n_missing > 0) {
  stop( "存在配体-受体标签匹配失败，请检查interaction_name。" )
}

# ============================================================
# 13. 设置配体-受体顺序
#     MIF最上面，SPP1第二
# ============================================================

lr_order <- as.character(
  LR_order_df$plot_id
)

# ggplot2的Y轴从下向上排列
# 因此需要反转顺序

y_levels <- rev(lr_order)

df$plot_id <- factor(
  as.character(df$plot_id),
  levels = y_levels
)

df$Group <- factor(
  as.character(df$Group),
  levels = group_order
)

# ============================================================
# 14. 补齐四组网格
# ============================================================

plot_df <- df %>%
  tidyr::complete(
    plot_id = factor(
      y_levels,
      levels = y_levels
    ),
    Group = factor(
      group_order,
      levels = group_order
    )
  )

# 关键修正：
# complete()之后再次设置因子顺序

plot_df$plot_id <- factor(
  as.character(plot_df$plot_id),
  levels = y_levels
)

plot_df$Group <- factor(
  as.character(plot_df$Group),
  levels = group_order
)

# ============================================================
# 15. 提取显著通讯气泡
# ============================================================

plot_points <- plot_df %>%
  filter(
    !is.na(prob),
    !is.na(pval),
    pval < 0.05
  )

plot_points$plot_id <- factor(
  as.character(plot_points$plot_id),
  levels = y_levels
)

# ============================================================
# 15.1 检查最终纵坐标顺序
# ============================================================

cat( "\n========== 最终纵坐标顺序 ==========\n" )

cat("\n最上面的10种配体-受体对：\n")

print( head( rev(levels(plot_df$plot_id)), 10 ) )

cat("\n最上面的配体-受体所属通路：\n")

print( head( LR_order_df$pathway_name, 10 ) )

# ============================================================
# 16. 检查绘图行数
# ============================================================

n_rows_with_points <- plot_points %>%
  distinct(plot_id) %>%
  nrow()

cat("\n========== 绘图检查 ==========\n")

cat( "应绘制配体-受体行数：", n_LR, "\n" )

cat( "实际有显著气泡的行数：", n_rows_with_points, "\n" )

cat( "实际显著气泡数量：", nrow(plot_points), "\n" )

if (n_rows_with_points != n_LR) {
  stop( "错误：存在整行没有显著气泡的配体-受体对。" )
}

cat( "检查通过：每行至少有一个显著气泡。\n" )

# ============================================================
# 17. 设置颜色
# ============================================================

bubble_colors <- c(
  "#4575B4",
  "#91BFDB",
  "#E0F3F8",
  "#FFFFBF",
  "#FEE090",
  "#FDAE61",
  "#F46D43",
  "#D73027",
  "#A50026"
)

# ============================================================
# 18. 设置通讯概率范围
# ============================================================

prob_range <- range(
  plot_points$prob,
  na.rm = TRUE
)

cat("\n通讯概率范围：\n")
print(prob_range)

# ============================================================
# 19. 绘制气泡图
# ============================================================

p <- ggplot(
  plot_df,
  aes( x = Group, y = plot_id )
) +

  # 网格保持紧凑
  geom_tile(
    width = 0.92,
    height = 0.92,
    fill = "white",
    color = "#E8E8E8",
    linewidth = 0.18
  ) +

  # 显著通讯气泡
  geom_point( data = plot_points, aes( color = prob, size = log10_p ), shape = 16, na.rm = TRUE ) +

  # 通讯概率颜色
  scale_color_gradientn( colours = bubble_colors, limits = prob_range, name = "Commun. Prob." ) +

  # 气泡大小
  scale_size_continuous(
    range = c(0.5, 3.2),
    limits = c(0, 10),
    breaks = c(2, 4, 6, 8, 10),
    name = expression(
      -log[10](p-value)
    )
  ) +

  # 横坐标保持不变
  scale_x_discrete( drop = FALSE, expand = expansion(add = 0.5) ) +

  # 关键修改：明确锁定Y轴顺序
  scale_y_discrete( limits = y_levels, drop = FALSE, expand = expansion(add = 0.5) ) +

  labs( x = NULL, y = NULL ) +

  theme_bw( base_size = 9 ) +

  theme(

    panel.grid.major = element_blank(),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),

    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 8,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 6,
      color = "black"
    ),

    axis.ticks = element_blank(),

    axis.title = element_blank(),

    legend.position = "right",

    legend.box = "vertical",

    legend.title = element_text(
      size = 9,
      color = "black"
    ),

    legend.text = element_text(
      size = 8,
      color = "black"
    ),

    legend.spacing.y = unit(
      0.4,
      "cm"
    ),

    plot.margin = margin(
      10, 10, 10, 10
    )
  ) +

  guides(

    color = guide_colorbar(
      order = 1,
      barheight = unit(3, "cm"),
      barwidth = unit(0.45, "cm")
    ),

    size = guide_legend(
      order = 2,
      override.aes = list(
        color = "black"
      )
    )
  )

# ============================================================
# 20. 显示图片
# ============================================================

print(p)

# ============================================================
# 21. 自动设置图片尺寸
# ============================================================

plot_height <- max(
  6,
  n_LR * 0.15 + 2.5
)

plot_width <- 4.6

cat( "\n图片尺寸：", plot_width, "x", plot_height, "英寸\n" )

# ============================================================
# 22. 保存PDF
# ============================================================

ggsave(
  filename = file.path(
    outdir,
    "PSME2_CellChat_6Pathways_LR_bubble.pdf"
  ),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

# ============================================================
# 23. 保存PNG
# ============================================================

ggsave(
  filename = file.path(
    outdir,
    "PSME2_CellChat_6Pathways_LR_bubble.png"
  ),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

# ============================================================
# 24. 保存全部绘图数据
# ============================================================

write.csv(
  as.data.frame(df),
  file = file.path(
    outdir,
    "PSME2_CellChat_6Pathways_LR_data.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 25. 保存去重显著配体-受体列表
# ============================================================

write.csv(
  as.data.frame(LR_order_df),
  file = file.path(
    outdir,
    "PSME2_CellChat_6Pathways_significant_LR_list.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 26. 输出各通路显著配体-受体数量
# ============================================================

cat("\n========== 各通路统计 ==========\n")

pathway_summary <- LR_order_df %>%
  count(
    pathway_name,
    name = "LR_count"
  )

print( as.data.frame(pathway_summary) )

# ============================================================
# 27. 最终输出信息
# ============================================================

cat("\n====================================\n")

cat( "六条通路全部通讯记录：", nrow(df_core), "\n" )

cat( "显著通讯记录：", nrow(df_sig), "\n" )

cat( "去重后显著配体-受体对：", n_LR, "\n" )

cat( "实际有气泡的行数：", n_rows_with_points, "\n" )

cat( "实际显著气泡数量：", nrow(plot_points), "\n" )

cat( "结果保存位置：", outdir, "\n" )

cat("====================================\n")

## ============================================================
## PSME2 scTenifoldKnk Virtual Knockout
## 保留原分析参数
## 结果保存至新的 scTenifoldKnk 文件夹
## ============================================================
rm(list=ls())
setwd('/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/')
library(Seurat)
library(scTenifoldKnk)
library(Matrix)
library(dplyr)

# ============================================================
# 1. 设置输出目录
# ============================================================

outdir <- "/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/fig4/scTenifoldKnk/4"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. 读取原始单细胞对象
# ============================================================

Epi <- readRDS(
  "/home/stu-xt/ccRcc_24.9.2/harmony_21/PSME2/Epithelial/30_1/Epi.rds"
)

DefaultAssay(Epi) <- "RNA"

# 查看细胞和基因数量
print(dim(Epi))

# 查看细胞分组
print(table(Epi$group, useNA = "ifany"))

# ============================================================
# 3. 设置虚拟敲除基因
# ============================================================

gKO <- "PSME2"

# 原代码中的机制相关基因
mechanistic_genes <- c(
  "SOCS1", "SOCS3", "CISH",
  "STAT3", "STAT6", "JAK1", "JAK2",
  "IL4", "IL10", "IL13",
  
  "MIF", "CD74",
  "CXCR2", "CXCR4",
  "ERK1", "ERK2",
  "AKT1", "HIF1A", "VEGFA",
  
  "SPP1", "CD44",
  "ITGAV", "ITGB1", "ITGB3", "ITGB5",
  "MMP9", "MMP2",
  "TWIST1", "SNAI1", "SNAI2"
)

# ============================================================
# 4. 提取3000个高变基因
# ============================================================

Epi <- FindVariableFeatures(
  Epi,
  selection.method = "vst",
  nfeatures = 3000
)

highVarGeneNames <- VariableFeatures(Epi)

# 合并高变基因、PSME2及机制相关基因
selected_genes <- unique(c(
  highVarGeneNames,
  gKO,
  mechanistic_genes
))

selected_genes <- intersect(
  selected_genes,
  rownames(Epi)
)

cat("纳入分析的基因数：", length(selected_genes), "\n")

# 保存纳入分析的基因列表
write.csv(
  data.frame(gene = selected_genes),
  file = file.path(outdir, "PSME2_selected_genes.csv"),
  row.names = FALSE
)

# ============================================================
# 5. 提取原始counts矩阵
# ============================================================

# Seurat v5
countMatrix <- GetAssayData(
  Epi,
  assay = "RNA",
  layer = "counts"
)

countMatrix_filtered <- countMatrix[
  selected_genes,
  ,
  drop = FALSE
]

cat("输入矩阵维度：\n")
print(dim(countMatrix_filtered))

stopifnot(gKO %in% rownames(countMatrix_filtered))

# ============================================================
# 6. 运行scTenifoldKnk
# ============================================================

set.seed(123456)

result <- scTenifoldKnk(
  countMatrix = countMatrix_filtered,
  gKO = gKO,
  qc = FALSE,
  nc_nNet = 10,
  nc_nCells = 500,
  nc_nComp = 3
)

# ============================================================
# 7. 保存完整虚拟敲除结果
# ============================================================

saveRDS(
  result,
  file = file.path(
    outdir,
    "PSME2_scTenifoldKnk_result.rds"
  )
)

# ============================================================
# 8. 导出差异调控统计结果
# ============================================================

dr <- result$diffRegulation

write.csv(
  dr,
  file = file.path(
    outdir,
    "PSME2_diffRegulation_all.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 9. 提取FDR显著扰动基因
# ============================================================

dr_sig <- dr %>%
  filter(
    gene != gKO,
    !is.na(p.adj),
    p.adj < 0.05
  ) %>%
  arrange(p.adj)

write.csv(
  dr_sig,
  file = file.path(
    outdir,
    "PSME2_diffRegulation_FDR005.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10. 单独检查IL6、MIF、SPP1
# ============================================================

target_results <- dr %>%
  filter(gene %in% c("IL6", "MIF", "SPP1"))

print(target_results)

write.csv(
  target_results,
  file = file.path(
    outdir,
    "PSME2_IL6_MIF_SPP1_results.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 11. 保存运行信息
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(outdir, "sessionInfo.txt")
)

cat("\n分析完成！\n")
cat("FDR < 0.05的扰动基因数量：", nrow(dr_sig), "\n")
cat("结果保存位置：", outdir, "\n")








## ============================================================
## PSME2 scTenifoldKnk 三张图
## 直接读取已经完成的虚拟敲除结果
## ============================================================

library(scTenifoldKnk)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)

# ============================================================
# 1. 读取你的虚拟敲除结果
# ============================================================

outdir <- "/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/fig4/scTenifoldKnk/2"

result <- readRDS(
  file.path(outdir, "PSME2_scTenifoldKnk_result.rds")
)


# ============================================================
# Figure 1. PSME2虚拟敲除调控网络图
# 采用推文中的plotKO代码
# ============================================================

plotKO(
  result,
  gKO = "PSME2",
  annotate = FALSE,
  q = 0.99
)


# ============================================================
# Figure 2. Top 20 Differentially Regulated Genes
# 采用推文中的柱状图代码
# ============================================================

df_KO <- result[["diffRegulation"]]

# 推文使用FC计算log2FC
# 先确认结果中存在FC列
print(colnames(df_KO))

if ("FC" %in% colnames(df_KO)) {
  
  df_KO$log2FC <- log2(df_KO$FC)
  
  top20 <- df_KO %>%
    filter(p.value < 0.05 & log2FC > 0) %>%
    arrange(desc(log2FC)) %>%
    slice_head(n = 20)
  
  p2 <- top20 %>%
    mutate(gene = fct_reorder(gene, log2FC)) %>%
    ggplot(aes(x = gene, y = log2FC, fill = log2FC)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    coord_flip() +
    geom_text(
      aes(label = round(log2FC, 1)),
      hjust = -0.15,
      size = 3.6,
      color = "grey20"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.08))
    ) +
    scale_fill_gradient(
      low = "#9ecae1",
      high = "#08519c"
    ) +
    labs(
      title = "Top 20 Differentially Regulated Genes",
      subtitle = paste0(
        "Ranked by log2FC (n = ",
        nrow(df_KO),
        ")"
      ),
      x = NULL,
      y = "log2FC"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      plot.subtitle = element_text(color = "grey40")
    )
  
  print(p2)
  
} else {
  
  message(
    "当前结果中没有FC列，无法直接复现推文的log2FC柱状图。"
  )
  
}


# ============================================================
# Figure 3. Differential Regulation Volcano Plot
# 采用推文中的火山图代码
# ============================================================

df <- result$diffRegulation

df$log_pval <- -log10(df$p.value)

label_genes <- subset(
  df,
  abs(Z) > 1 & p.value < 0.05
)

p3 <- ggplot(df, aes(x = Z, y = log_pval)) +
  geom_point(alpha = 0.5) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = c(2),
    linetype = "dashed",
    color = "blue"
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 50
  ) +
  labs(
    title = "Differential Regulation Volcano Plot",
    x = "Z-score",
    y = "-log10(p-value)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 14,
      color = "black",
      angle = 0,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black",
      hjust = 1
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    legend.text = element_text(size = 14),
    legend.title = element_text(
      size = 14,
      face = "bold"
    )
  )

print(p3)

















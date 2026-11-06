
rm(list = ls())
library(Seurat)
library(CytoTRACE)
library(qs)
library(ggplot2)
library(dplyr)

# 1. 设置路径
indir <- "/home/stu-xt/project/work/ccRcc_mac_psme2/fig5"
outdir <- file.path(indir, "CytoTRACE")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 2. 读取巨噬细胞对象
Mac <- qs::qread(file.path(indir, "Mac_UCell.qs"))
stopifnot(ncol(Mac) == 7380)

# 3. 设置巨噬细胞亚群
stopifnot("mac_subtype" %in% colnames(Mac@meta.data))
Mac$mac_subtype <- factor(
  Mac$mac_subtype,
  levels = paste0("C", 0:12)
)
Idents(Mac) <- "mac_subtype"

# 4. 提取RNA counts
exp1 <- GetAssayData(Mac, assay = "RNA", layer = "counts")
exp1 <- exp1[Matrix::rowSums(exp1 > 0) >= 5, , drop = FALSE]
exp1 <- as.matrix(exp1)

# 5. CytoTRACE计算
set.seed(123)
results <- CytoTRACE(exp1, ncores = 20)

# 6. 保存原始计算结果
save(results, file = file.path(outdir, "CytoTRACE_results.RData"))
saveRDS(results, file.path(outdir, "CytoTRACE_results.rds"))

# 7. 提取CytoTRACE评分
scores <- results$CytoTRACE
stopifnot(length(scores) == ncol(Mac))
if (!is.null(names(scores))) {
  scores <- scores[colnames(Mac)]
} else {
  names(scores) <- colnames(Mac)
}
stopifnot(!anyNA(scores))

Mac$CytoTRACE <- as.numeric(scores)
qs::qsave(Mac, file.path(outdir, "Mac_CytoTRACE.qs"))

# 8. 导出评分表
score_df <- data.frame(
  Cell = colnames(Mac),
  Cluster = as.character(Mac$mac_subtype),
  CytoTRACE = Mac$CytoTRACE
)
write.csv(
  score_df,
  file.path(outdir, "Mac_CytoTRACE_scores.csv"),
  row.names = FALSE
)

# 9. 绘制各亚群CytoTRACE评分
p1 <- ggplot(
  score_df,
  aes(x = Cluster, y = CytoTRACE, fill = Cluster)
) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA,
               fill = "white") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = "Macrophage cluster", y = "CytoTRACE score")

ggsave(
  file.path(outdir, "Mac_CytoTRACE_Violin.pdf"),
  p1, width = 9, height = 5
)

# 10. 绘制CytoTRACE降维图
if ("umap" %in% Reductions(Mac)) {
  p2 <- FeaturePlot(
    Mac,
    features = "CytoTRACE",
    reduction = "umap",
    cols = c("#DCEAF5", "#08519C"),
    order = TRUE
  ) + ggtitle("CytoTRACE")

  ggsave(
    file.path(outdir, "Mac_CytoTRACE_UMAP.pdf"),
    p2, width = 6, height = 5
  )
}

# 11. CytoTRACE原始绘图函数
phenot <- as.character(Mac$mac_subtype)
names(phenot) <- colnames(Mac)

if ("tsne" %in% Reductions(Mac)) {
  emb <- Embeddings(Mac, reduction = "tsne")
  plotCytoTRACE(
    results,
    phenotype = phenot,
    emb = emb,
    outputDir = file.path(outdir, "CytoTRACE_original")
  )
}

plotCytoGenes(
  results,
  numOfGenes = 30,
  outputDir = outdir
)

# 12. 保存软件版本
writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)

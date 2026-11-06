
rm(list = ls())
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(RColorBrewer)

# 1. 设置工作目录
outdir <- "/home/stu-xt/project/work/slingshot_env"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

# 2. 读取之前保存的PCA和聚类数据
x <- readRDS(file.path(outdir, "slingshot_input.rds"))
pca_mat <- as.matrix(x$pca)
clusters <- as.character(x$clusters)
names(clusters) <- rownames(pca_mat)

# 3. 保留C0-C12共13个巨噬细胞亚群
keep <- clusters %in% paste0("C", 0:12)
pca_7380 <- pca_mat[keep, 1:30, drop = FALSE]
clusters_7380 <- clusters[keep]
clusters_7380 <- factor(clusters_7380, levels = paste0("C", 0:12))

# 4. 检查细胞数量及PCA维度
print(dim(pca_7380))
print(table(clusters_7380))
stopifnot(nrow(pca_7380) == 7380)
stopifnot(ncol(pca_7380) == 30)
stopifnot(!anyNA(pca_7380), !anyNA(clusters_7380))

# 5. Slingshot轨迹分析（PCA前30维）
set.seed(42)
slingshot_7380 <- slingshot(
  pca_7380,
  clusterLabels = clusters_7380
)

# 6. 检查轨迹数量及各轨迹经过的亚群
lineages <- slingLineages(slingshot_7380)
print(length(lineages))
print(lineages)
stopifnot(length(lineages) == 3)

# 7. 保存Slingshot分析对象
saveRDS(
  slingshot_7380,
  file.path(outdir, "slingshot_PCA_7380_3lineages.rds")
)

# 8. 提取三条轨迹的拟时序
pt_7380 <- slingPseudotime(slingshot_7380)
colnames(pt_7380) <- paste0("Lineage", 1:ncol(pt_7380))
rownames(pt_7380) <- rownames(pca_7380)
print(dim(pt_7380))

# 9. 保存拟时序结果
saveRDS(pt_7380, file.path(outdir, "slingshot_pseudotime_7380.rds"))
write.csv(pt_7380, file.path(outdir, "slingshot_pseudotime_7380.csv"))

# 10. 绘制PCA空间中的三条轨迹
cluster_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(13),
  paste0("C", 0:12)
)

pdf(file.path(outdir, "Mac_Slingshot_PCA_3lineages.pdf"), width = 8, height = 6)
plot(
  pca_7380[, 1:2],
  col = cluster_cols[as.character(clusters_7380)],
  pch = 16, cex = 0.35,
  xlab = "PC1", ylab = "PC2",
  main = "Macrophage Slingshot Trajectories"
)
curves <- slingCurves(slingshot_7380)
curve_cols <- c("#E64B35", "#4DBBD5", "#00A087")
for (i in seq_along(curves)) {
  lines(curves[[i]]$s[, 1:2], col = curve_cols[i], lwd = 2.5)
}
legend("topright", legend = paste0("Lineage ", 1:3),
       col = curve_cols, lwd = 2.5, bty = "n")
dev.off()

# 11. 保存
save.image(file.path(outdir, "Mac_Slingshot_Analysis.RData"))
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "Slingshot_sessionInfo.txt"))

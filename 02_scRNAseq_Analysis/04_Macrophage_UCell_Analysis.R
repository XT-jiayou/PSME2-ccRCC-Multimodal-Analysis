rm(list = ls())
library(Seurat); library(UCell); library(qs); library(ggplot2)
indir <- "/home/stu-xt/project/work/ccRcc_mac_psme2/fig5"
outdir <- "/home/stu-xt/ccRcc_24.9.2/ZC/Epi/scissor/WSI/macrophage/UCell2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
Mac <- qs::qread(file.path(indir, "Mac_subtype.qs"))
load(file.path(indir, "m1m2.Rdata"))
print(names(m1m2_pws))
M1 <- m1m2_pws[["GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP"]]
M2 <- m1m2_pws[["GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN"]]
counts <- Seurat::GetAssayData(Mac, assay = "RNA", layer = "counts")
signatures <- list(M1 = intersect(M1, rownames(counts)), M2 = intersect(M2, rownames(counts)))
print(lengths(signatures))
scores <- UCell::ScoreSignatures_UCell(matrix = counts, features = signatures, chunk.size = 500, ncores = 1)
Mac <- Seurat::AddMetaData(Mac, metadata = scores)
Mac$M2_minus_M1_UCell <- Mac$M2_UCell - Mac$M1_UCell
qs::qsave(Mac, file.path(outdir, "Mac_UCell.qs"))
score_df <- data.frame(Cell = colnames(Mac), Mac@meta.data[, c("M1_UCell", "M2_UCell", "M2_minus_M1_UCell"), drop = FALSE], check.names = FALSE)
write.csv(score_df, file.path(outdir, "Mac_UCell_scores.csv"), row.names = FALSE)
p <- FeaturePlot(Mac, features = c("M1_UCell", "M2_UCell", "M2_minus_M1_UCell"), reduction = "tsne", cols = c("#08306B", "#FFFFFF", "#A50026"), order = TRUE, ncol = 3, combine = TRUE)
ggsave(file.path(outdir, "Figure5B_UCell_wide.pdf"), plot = p, width = 12, height = 6)

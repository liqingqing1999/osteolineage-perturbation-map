# vc_65a2 (S2 part 2): R UCell 打分 —— 与 vc_55 逐字同函数
# ScoreSignatures_UCell(M, features = gs, maxRank = 1500)，M = 原始 counts（秩不变）
# 预注册 30_ §一.1；产出 T65/export/ucell_scores.csv
suppressPackageStartupMessages({ library(Matrix); library(UCell) })

EXP <- Sys.getenv("T65_EXPORT")
if (EXP == "") stop("Set T65_EXPORT to the exported-matrix directory (see README)")
t0 <- Sys.time()

done <- file.path(EXP, "ucell_scores.csv")
if (file.exists(done)) {
  cat("[vc_65a2] scores 已存在，幂等跳过\n"); quit(save = "no")
}

M <- readMM(file.path(EXP, "matrix.mtx.gz"))
genes <- readLines(file.path(EXP, "genes.tsv"))
bars  <- readLines(file.path(EXP, "barcodes.tsv"))
M <- as(M, "dgCMatrix")
M <- Matrix::t(M)   # scipy 导出为 细胞×基因；UCell 需要 基因×细胞（同 vc_55）
rownames(M) <- genes; colnames(M) <- bars
stopifnot(nrow(M) == length(genes), ncol(M) == length(bars))
cat("[vc_65a2] M:", dim(M)[1], "x", dim(M)[2], "read", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

gs <- list(
  OS_identity = c("RUNX2","SP7","ALPL","BGLAP","SPP1","IBSP","DMP1",
                  "COL1A1","SOST","PHEX","BMP2","BMP4"),
  AD_drift    = c("PPARG","CEBPA","CEBPB","ADIPOQ","LEP","FABP4",
                  "PLIN1","LPL","GPD1","RETN","ACSL1")
)
missing <- unlist(lapply(gs, function(g) setdiff(g, genes)))
stopifnot(length(missing) == 0)
gs <- lapply(gs, function(g) intersect(g, genes))
cat("[vc_65a2] 基因集全命中：AD", length(gs$AD_drift), "/ OS", length(gs$OS_identity), "\n")

sc <- ScoreSignatures_UCell(M, features = gs, maxRank = 1500)
out <- data.frame(cell = bars,
                  AD_UCell = as.numeric(sc[, "AD_drift_UCell"]),
                  OS_UCell = as.numeric(sc[, "OS_identity_UCell"]))
write.csv(out, done, row.names = FALSE)
cat("[vc_65a2] 完成 n =", nrow(out), "，耗时", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

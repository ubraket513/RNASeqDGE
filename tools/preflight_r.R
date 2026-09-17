# Offline pinned environment smoke check. Called only by tools/toolchain.R.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]), getRversion() == "4.5.3")
pins <- c(DESeq2 = "1.50.2", apeglm = "1.32.0", BiocParallel = "1.44.0", pheatmap = "1.0.13")
for (name in names(pins)) {
  stopifnot(as.character(packageVersion(name)) == pins[[name]])
  suppressPackageStartupMessages(library(name, character.only = TRUE))
}
set.seed(17092026)
cd <- data.frame(condition = factor(rep(c("control", "treated"), each = 3),
  levels = c("control", "treated")), row.names = paste0("sample_", 1:6))
x <- matrix(rnbinom(6000, mu = rep(seq(40, 440, length.out = 1000), 6), size = 10),
  nrow = 1000, dimnames = list(paste0("gene_", 1:1000), rownames(cd)))
x[1:100, 4:6] <- x[1:100, 4:6] * 3L
dds <- DESeq(DESeqDataSetFromMatrix(x, cd, ~ condition), quiet = TRUE, parallel = FALSE)
raw <- results(dds, contrast = c("condition", "treated", "control"), alpha = 0.05)
shrunk <- lfcShrink(dds, coef = "condition_treated_vs_control", res = raw,
                    type = "apeglm", quiet = TRUE, parallel = FALSE)
stopifnot(nrow(raw) == 1000, identical(rownames(raw), rownames(shrunk)),
          all(is.finite(shrunk$log2FoldChange)), median(raw$log2FoldChange[1:100]) > 0)
png(file.path(args[1], "MA.png"), width = 800, height = 600)
plotMA(shrunk, alpha = 0.05)
dev.off()
stopifnot(file.info(file.path(args[1], "MA.png"))$size > 0)
writeLines(capture.output(sessionInfo()), file.path(args[1], "sessionInfo.txt"))
cat("PASS: pinned R packages, DESeq2, apeglm, PNG\n")

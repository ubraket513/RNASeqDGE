# Offline backend smoke test, not statistical parity with the original study.
packages <- c("DESeq2", "apeglm", "GEOquery", "optparse", "tidyverse",
              "EnhancedVolcano", "RColorBrewer", "pheatmap", "ggrepel",
              "org.Hs.eg.db", "BiocParallel")
for (package in packages) {
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}
out <- "tests/output/p0/r"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(17092026)
sample_data <- data.frame(
  condition = factor(rep(c("control", "treated"), each = 3),
                     levels = c("control", "treated")),
  row.names = paste0("sample_", seq_len(6)))
means <- rep(seq(40, 440, length.out = 1000), 6)
counts <- matrix(rnbinom(6000, mu = means, size = 10), nrow = 1000,
                 dimnames = list(paste0("synthetic_", seq_len(1000)), rownames(sample_data)))
counts[1:100, 4:6] <- counts[1:100, 4:6] * 3L
dds <- DESeqDataSetFromMatrix(counts, sample_data, ~ condition)
dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)
raw <- results(dds, contrast = c("condition", "treated", "control"), alpha = 0.05)
coefficient <- "condition_treated_vs_control"
stopifnot(coefficient %in% resultsNames(dds))
shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm", quiet = TRUE)
stopifnot(nrow(raw) == 1000, nrow(shrunk) == 1000,
          identical(rownames(raw), rownames(shrunk)),
          all(is.finite(shrunk$log2FoldChange)),
          median(raw$log2FoldChange[1:100]) > 0)
png(file.path(out, "MA.png"), width = 800, height = 600)
plotMA(raw, main = "Synthetic treated vs control")
dev.off()
stopifnot(file.info(file.path(out, "MA.png"))$size > 0)
write.csv(as.data.frame(raw), file.path(out, "unshrunken.csv"))
write.csv(as.data.frame(shrunk), file.path(out, "shrunken.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("PASS: required libraries, DESeq2, explicit contrast, apeglm, PNG device\n")

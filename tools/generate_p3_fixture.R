# Deterministic synthetic fixture; no network or external reference data.
set.seed(17092026)
out <- "tests/fixtures/p3"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = TRUE, row.names = FALSE,
  na = "", qmethod = "double")
samples <- data.frame(sample_id = paste0("sample_", 1:8), condition = rep(c("control", "treated"), each = 4),
  batch = rep(c("A", "B"), 4))
counts <- matrix(rnbinom(8000, mu = rep(seq(40, 440, length.out = 1000), 8), size = 20), 1000, 8)
counts[1:100, 5:8] <- counts[1:100, 5:8] * 5L
counts[101:200, 1:4] <- counts[101:200, 1:4] * 5L
counts[, c(2, 4, 6, 8)] <- counts[, c(2, 4, 6, 8)] * 2L
counts[201:300, 5:8] <- counts[201:300, 5:8] * 2L
counts[501:999, ] <- 0L
for (i in 501:999) counts[i, 1L + (i%%8L)] <- 1L
counts[1000, ] <- 0L
colnames(counts) <- samples$sample_id
write_tsv(data.frame(gene_id = sprintf("gene_%04d", 1:1000), counts, check.names = FALSE), "counts.tsv")
write_tsv(samples, "samples.tsv")
write_tsv(data.frame(key = c("version", "design_terms", "alpha", "filter", "shrinkage", "type.batch",
  "reference.batch"), value = c("1", "batch,condition", "0.05", "zero_total", "apeglm", "categorical",
  "A")), "analysis.tsv")
write_tsv(data.frame(contrast_id = c("treated_vs_control", "control_vs_treated"), factor = "condition",
  numerator = c("treated", "control"), denominator = c("control", "treated")), "contrasts.tsv")
write_tsv(data.frame(gene_id = sprintf("gene_%04d", 1:500), gene_symbol = c("SHARED", "SHARED", paste0("S",
  3:500)), biotype = "synthetic", chromosome = "fixture"), "annotation.tsv")

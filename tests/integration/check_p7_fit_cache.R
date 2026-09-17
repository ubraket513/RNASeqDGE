# Focused public-CLI regression: compatible models reuse fits, relevelled models do not.
# The oracle below independently fits every requested contrast with literal settings.
rscript <- normalizePath(file.path(R.home("bin"), "Rscript"))
work <- tempfile("p7-cache-", tmpdir = normalizePath("tests/output"))
dir.create(work)
fixture <- file.path(work, "input"); dir.create(fixture)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA", qmethod = "double")
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
set.seed(7017)
samples <- data.frame(sample_id = paste0("s", 1:12), condition = rep(c("A", "B", "C"), each = 4))
x <- matrix(rnbinom(600 * 12, mu = rep(seq(30, 500, length.out = 600), 12), size = 12),
  nrow = 600, dimnames = list(sprintf("gene_%04d", 1:600), samples$sample_id))
x[1:100, 5:8] <- x[1:100, 5:8] * 3L
x[101:200, 9:12] <- x[101:200, 9:12] * 4L
x[600, ] <- 0L
write_tsv(data.frame(gene_id = rownames(x), x, check.names = FALSE), file.path(fixture, "counts.tsv"))
write_tsv(samples, file.path(fixture, "samples.tsv"))
write_tsv(data.frame(key = c("version", "design_terms", "alpha", "filter", "shrinkage"),
  value = c("1", "condition", "0.05", "zero_total", "apeglm")), file.path(fixture, "analysis.tsv"))
contrasts <- data.frame(contrast_id = c("B_vs_A", "C_vs_A", "A_vs_B"),
  factor = "condition", numerator = c("B", "C", "A"), denominator = c("A", "A", "B"))
write_tsv(contrasts, file.path(fixture, "contrasts.tsv"))
out <- file.path(work, "out")
args <- c("--vanilla", "run_deg_analysis_offline.R", "--counts", file.path(fixture, "counts.tsv"),
  "--samples", file.path(fixture, "samples.tsv"), "--analysis", file.path(fixture, "analysis.tsv"),
  "--contrasts", file.path(fixture, "contrasts.tsv"), "--out", out, "--workers", "1")
status <- system2(rscript, vapply(args, shQuote, ""), stdout = file.path(work, "cli.log"), stderr = file.path(work, "cli.log"))
stopifnot(status == 0L)
cat("Artifacts:", work, "\n")
stopifnot(file.exists(file.path(out, "timings.tsv")))
timings <- read_tsv(file.path(out, "timings.tsv"))
stopifnot(all(c("stage", "contrast_id", "model_id", "status", "workers", "wall_seconds", "cpu_seconds", "self_peak_rss_kb") %in% names(timings)),
          all(timings$wall_seconds >= 0), all(timings$cpu_seconds >= 0))
fit <- timings[timings$stage == "fit", ]
stopifnot(nrow(fit) == 3L, identical(fit$status, c("executed", "reused", "executed")),
          fit$model_id[1] == fit$model_id[2], fit$model_id[1] != fit$model_id[3])
for (stage in c("normalized_counts", "vst", "pca", "sample_distances")) {
  rows <- timings[timings$stage == stage, ]
  stopifnot(nrow(rows) == 3L, identical(rows$status, c("executed", "reused", "executed")))
}
cat("PASS: compatible models reuse fit/transforms; inverse baseline remains independent\n")
if (Sys.getenv("P7_TIMINGS_ONLY") == "1") quit(status = 0L)
suppressPackageStartupMessages(library(DESeq2))
x <- x[1:599, , drop = FALSE]
close_num <- function(actual, reference, label) {
  actual <- as.numeric(actual); reference <- as.numeric(reference)
  stopifnot(length(actual) == length(reference), identical(is.na(actual), is.na(reference)))
  keep <- !is.na(reference)
  if (any(abs(actual[keep] - reference[keep]) > 1e-10 + 1e-7 * abs(reference[keep]))) stop("numeric mismatch: ", label)
}
for (i in 1:3) {
  row <- contrasts[i, ]
  levels <- if (i == 3L) c("B", "A", "C") else c("A", "B", "C")
  cd <- data.frame(condition = factor(samples$condition, levels = levels), row.names = samples$sample_id)
  dds <- DESeq(DESeqDataSetFromMatrix(x, cd, ~ condition), quiet = TRUE, parallel = FALSE)
  raw <- results(dds, contrast = c("condition", row$numerator, row$denominator), alpha = 0.05)
  set.seed(1)
  shrunk <- lfcShrink(dds, coef = paste0("condition_", row$numerator, "_vs_", row$denominator),
    res = raw, type = "apeglm", quiet = TRUE, parallel = FALSE)
  vst <- SummarizedExperiment::assay(varianceStabilizingTransformation(dds, blind = FALSE))
  root <- file.path(out, "contrasts", row$contrast_id)
  artifact <- function(name) read_tsv(file.path(root, paste0(name, ".tsv")))
  for (pair in list(list("results", raw), list("shrunken", shrunk))) {
    got <- artifact(pair[[1]]); expected <- as.data.frame(pair[[2]])
    stopifnot(identical(got$gene_id, rownames(expected)))
    for (name in names(expected)) close_num(got[[name]], expected[[name]], paste(row$contrast_id, pair[[1]], name))
  }
  close_num(as.matrix(artifact("normalized_counts")[-1]), counts(dds, normalized = TRUE), "normalization")
  stopifnot(identical(artifact("model_matrix")$sample_id, samples$sample_id))
  close_num(as.matrix(artifact("model_matrix")[-1]), model.matrix(~ condition, cd), "model matrix")
  significant <- !is.na(raw$padj) & raw$padj < 0.05
  volcano <- artifact("volcano_data")
  stopifnot(identical(volcano$significant, significant))
  close_num(volcano$log2FoldChange, shrunk$log2FoldChange, "volcano LFC")
  close_num(volcano$padj, raw$padj, "volcano padj")
  close_num(volcano$neg_log10_padj, -log10(pmax(raw$padj, .Machine$double.xmin)), "volcano y")
  pca <- prcomp(t(vst), center = TRUE, scale. = FALSE)
  for (k in 1:2) {
    got <- artifact("pca_data")[[paste0("PC", k)]]
    sign <- if (sum(got * pca$x[, k]) < 0) -1 else 1
    close_num(sign * got, pca$x[, k], "PCA coordinates")
    close_num(artifact("pca_data")[[paste0("PC", k, "_variance")]],
      rep(pca$sdev[k]^2 / sum(pca$sdev^2), nrow(samples)), "PCA variance")
  }
  close_num(as.matrix(artifact("sample_distances")[-1]), as.matrix(dist(t(vst))), "sample distances")
  ids <- which(significant); ids <- head(ids[order(raw$padj[ids], rownames(raw)[ids], method = "radix")], 50L)
  stopifnot(identical(artifact("heatmap_data")$gene_id, rownames(x)[ids]))
  close_num(as.matrix(artifact("heatmap_data")[-1]), vst[ids, , drop = FALSE], "heatmap VST")
  stopifnot(length(list.files(root, pattern = "\\.png$")) == 7L)
  cat("PASS independent per-contrast oracle:", row$contrast_id, "\n")
  rm(dds, raw, shrunk, vst); gc()
}
writeLines("PASS: strict numeric/NA/identity/plot-data parity; two fits for three contrasts", file.path(work, "VERIFIED"))
cat("PASS P7 focused cache oracle:", work, "\n")

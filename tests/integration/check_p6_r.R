# Independent real-data oracle: never sources the production analysis script.
# Usage: Rscript --vanilla tests/integration/check_p6_r.R INPUT_DIR OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(DESeq2))
samples <- read_tsv(file.path(args[1], "samples.tsv"))
input <- read_tsv(file.path(args[1], "counts.tsv"))
stopifnot(identical(names(input)[-1], samples$sample_id), nrow(input) == 47886L)
x <- as.matrix(input[-1]); rownames(x) <- input$gene_id
x <- x[rowSums(x) > 0, , drop = FALSE]
cd <- data.frame(condition = factor(samples$condition, levels = c("control", "bipolar")),
                 row.names = samples$sample_id)
dds <- DESeq(DESeqDataSetFromMatrix(x, cd, ~ condition), quiet = TRUE, parallel = FALSE)
raw <- results(dds, contrast = c("condition", "bipolar", "control"), alpha = 0.05)
set.seed(1)
shrunk <- lfcShrink(dds, coef = "condition_bipolar_vs_control", res = raw,
                   type = "apeglm", quiet = TRUE, parallel = FALSE)
out <- file.path(args[2], "contrasts", "bipolar_vs_control")
compare <- function(actual, expected, label) {
  # apeglm vectors can carry gene names; identity is checked separately above.
  # Compare numeric masks/values without incidental vector-name attributes.
  actual <- as.numeric(actual); expected <- as.numeric(expected)
  stopifnot(length(actual) == length(expected), identical(is.na(actual), is.na(expected)))
  keep <- !is.na(expected)
  if (any(abs(actual[keep] - expected[keep]) > 1e-10 + 1e-7 * abs(expected[keep])))
    stop("numeric mismatch: ", label)
  if (label %in% c("pvalue", "padj")) {
    stopifnot(identical(actual[keep] == 0, expected[keep] == 0))
    positive <- keep & expected > 0
    stopifnot(all(abs(log(actual[positive]) - log(expected[positive])) <= 1e-7))
  }
}
for (name in c("results", "shrunken")) {
  actual <- read_tsv(file.path(out, paste0(name, ".tsv")))
  expected <- if (name == "results") raw else shrunk
  stopifnot(identical(actual$gene_id, rownames(expected)))
  for (column in names(as.data.frame(expected))) compare(actual[[column]], expected[[column]], column)
  stopifnot(identical(actual$gene_id[!is.na(actual$padj) & actual$padj < 0.05],
                      rownames(expected)[!is.na(expected$padj) & expected$padj < 0.05]))
}
norm <- read_tsv(file.path(out, "normalized_counts.tsv"))
stopifnot(identical(norm$gene_id, rownames(x)), identical(names(norm)[-1], samples$sample_id))
compare(as.numeric(as.matrix(norm[-1])), as.numeric(counts(dds, normalized = TRUE)), "normalized counts")
pngs <- list.files(out, pattern = "\\.png$", full.names = TRUE)
stopifnot(length(pngs) == 7L)
for (path in pngs) {
  con <- file(path, "rb"); signature <- readBin(con, "raw", n = 8); close(con)
  stopifnot(identical(as.integer(signature), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
}
cat("PASS: real-data independent raw/shrunken/normalization oracle, gene/sample identity, NA/zero masks, log-p values and DEG membership; seven PNG signatures.\n")

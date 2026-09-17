# Integration expectations deliberately run the public CLI, never source it.
rscript <- normalizePath(file.path(R.home("bin"), "Rscript"))
fixture <- "tests/fixtures/p3"
dir.create("tests/output/p3", recursive = TRUE, showWarnings = FALSE)
work <- tempfile("p3-integration-", tmpdir = normalizePath("tests/output/p3"))
dir.create(work)
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
                                         row.names = FALSE, na = "NA", qmethod = "double")
invoke <- function(out, workers = 1, dir = fixture, extra = character(), annotation = TRUE) {
  args <- c("--vanilla", "run_deg_analysis_offline.R",
    "--counts", file.path(dir, "counts.tsv"), "--samples", file.path(dir, "samples.tsv"),
    "--analysis", file.path(dir, "analysis.tsv"), "--contrasts", file.path(dir, "contrasts.tsv"),
    if (annotation) c("--annotation", file.path(dir, "annotation.tsv")), "--out", out, "--workers", workers, extra)
  log <- suppressWarnings(system2(rscript, vapply(args, shQuote, ""), stdout = TRUE, stderr = TRUE))
  writeLines(log, file.path(work, paste0(basename(out), "-workers-", workers, ".log")))
  status <- attr(log, "status")
  list(status = if (is.null(status)) 0L else status, log = log)
}
# Output schemas must retain unique column names for every valid input.
# Run this small boundary suite alone with P3_SCHEMA_ONLY=1.
check_output_schema_names <- function() {
  cases <- c("PC1", "PC2", "PC1_variance", "PC2_variance", "sample_id", "model_sample_id")
  for (reserved in cases) {
    dir <- file.path(work, paste0("reserved-", reserved)); dir.create(dir)
    invisible(file.copy(list.files(fixture, pattern = "\\.tsv$", full.names = TRUE), dir))
    samples <- read_tsv(file.path(dir, "samples.tsv"))
    analysis <- read_tsv(file.path(dir, "analysis.tsv"))
    counts <- read_tsv(file.path(dir, "counts.tsv"))
    if (reserved == "sample_id") {
      samples$sample_id[1] <- "sample_id"
      names(counts)[2] <- "sample_id"
    } else if (reserved == "model_sample_id") {
      samples$sample <- rep(c("_id", "_id", "other", "other"), 2)
      analysis$value[analysis$key == "design_terms"] <- "batch,condition,sample"
      analysis <- rbind(analysis, data.frame(key = c("type.sample", "reference.sample"),
                                            value = c("categorical", "other")))
    } else {
      samples[[reserved]] <- seq_len(nrow(samples))
      analysis <- rbind(analysis, data.frame(key = paste0("type.", reserved), value = "numeric"))
    }
    write_tsv(samples, file.path(dir, "samples.tsv"))
    write_tsv(analysis, file.path(dir, "analysis.tsv"))
    write_tsv(counts, file.path(dir, "counts.tsv"))
    target <- file.path(work, paste0("reserved-", reserved, "-out"))
    run <- invoke(target, dir = dir)
    if (run$status == 0L || dir.exists(target))
      stop("reserved output-column regression: accepted ", reserved, " and published ambiguous TSV headers")
    stopifnot(any(grepl("samples.tsv.*collid|samples.tsv.*reserved", run$log)),
              any(grepl(if (reserved == "model_sample_id") "sample_id" else reserved, run$log, fixed = TRUE)))
    cat("PASS: preflight rejects output-column collision", reserved, "\n")
  }
  stopifnot(length(list.files(work, pattern = "^\\.p3-stage", all.files = TRUE)) == 0L)
}
check_output_schema_names()
if (identical(Sys.getenv("P3_SCHEMA_ONLY"), "1")) {
  positive <- file.path(work, "schema-positive")
  run <- invoke(positive)
  if (run$status != 0L) stop("schema regression positive control failed:\n", paste(run$log, collapse = "\n"))
  tables <- list.files(positive, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  stopifnot(length(tables) == 18L,
            all(vapply(tables, function(path) !anyDuplicated(names(read_tsv(path))), TRUE)))
  cat("PASS: ordinary manifest publishes 18 tables with unique headers\n")
  cat("PASS: focused output-schema regressions; artifacts:", work, "\n")
  quit(status = 0L)
}
# RED: absent implementation must fail this success assertion, with R's missing-file diagnostic.
out <- file.path(work, "workers-1")
run <- invoke(out)
if (run$status != 0L) stop("positive offline invocation failed:\n", paste(run$log, collapse = "\n"))
stopifnot(file.exists(file.path(out, "contrasts", "treated_vs_control", "results.tsv")))
cat("PASS: offline public command succeeds\n")
close_num <- function(actual, expected, label) {
  actual <- as.numeric(actual); expected <- as.numeric(expected)
  if (length(actual) != length(expected) || !identical(is.na(actual), is.na(expected)) ||
      any(abs(actual[!is.na(actual)] - expected[!is.na(expected)]) >
          1e-10 + 1e-7 * abs(expected[!is.na(expected)]))) stop("numeric mismatch: ", label)
}
# Independent oracle fixes the model, levels, filter and comparison literally.
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))
samples <- read_tsv(file.path(fixture, "samples.tsv"))
input <- read_tsv(file.path(fixture, "counts.tsv"))
x <- as.matrix(input[-1]); rownames(x) <- input$gene_id
x <- x[1:999, , drop = FALSE]
cd <- data.frame(batch = factor(samples$batch, levels = c("A", "B")),
                 condition = factor(samples$condition, levels = c("control", "treated")),
                 row.names = samples$sample_id)
dds <- DESeqDataSetFromMatrix(x, cd, ~ batch + condition)
dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)
raw <- results(dds, contrast = c("condition", "treated", "control"), alpha = 0.05)
set.seed(1)
shrunk <- lfcShrink(dds, coef = "condition_treated_vs_control", res = raw,
                    type = "apeglm", quiet = TRUE, parallel = FALSE)
v <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
# PCA uses all retained genes, making the exact assay auditable.
pca <- prcomp(t(v), center = TRUE, scale. = FALSE)
expected_genes <- sprintf("gene_%04d", 1:999)
expected_sig <- expected_genes[!is.na(raw$padj) & raw$padj < 0.05]
stopifnot(length(expected_sig) > 100L, median(raw$log2FoldChange[1:100]) > 1.5,
          anyNA(raw$padj), any(!is.na(raw$pvalue) & is.na(raw$padj)))
root1 <- file.path(out, "contrasts", "treated_vs_control")
read_artifact <- function(root, name) read_tsv(file.path(root, paste0(name, ".tsv")))
result <- read_artifact(root1, "results")
shrink <- read_artifact(root1, "shrunken")
for (tab in list(result, shrink)) {
  stopifnot(identical(tab$gene_id, expected_genes), nrow(tab) == 999L,
            identical(tab$gene_symbol[1:2], c("SHARED", "SHARED")),
            all(is.na(tab$gene_symbol[501:999])))
}
for (name in names(as.data.frame(raw))) close_num(result[[name]], raw[[name]], paste("raw", name))
for (name in names(as.data.frame(shrunk))) close_num(shrink[[name]], shrunk[[name]], paste("shrink", name))
check_artifacts <- function(root) {
  expected <- c("results.tsv", "shrunken.tsv", "normalized_counts.tsv", "model_matrix.tsv",
    "summary.tsv", "volcano_data.tsv", "pca_data.tsv", "sample_distances.tsv", "heatmap_data.tsv",
    "MA_raw.png", "MA_shrunken.png", "dispersion.png", "volcano.png", "pca.png",
    "sample_distances.png", "heatmap.png", "sessionInfo.txt")
  stopifnot(setequal(list.files(root), expected), all(file.info(file.path(root, expected))$size > 0))
  for (name in expected[grepl("\\.tsv$", expected)])
    stopifnot(!anyDuplicated(names(read_tsv(file.path(root, name)))))
  for (name in expected[grepl("\\.png$", expected)]) {
    con <- file(file.path(root, name), "rb"); signature <- readBin(con, "raw", n = 8); close(con)
    stopifnot(identical(as.integer(signature), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
  }
  stopifnot(any(grepl("DESeq2", readLines(file.path(root, "sessionInfo.txt")))))
}
for (id in c("treated_vs_control", "control_vs_treated")) check_artifacts(file.path(out, "contrasts", id))
norm <- read_artifact(root1, "normalized_counts")
stopifnot(identical(norm$gene_id, expected_genes), identical(names(norm)[-1], samples$sample_id))
close_num(as.matrix(norm[-1]), counts(dds, normalized = TRUE), "normalized counts")
model <- read_artifact(root1, "model_matrix")
stopifnot(identical(model$sample_id, samples$sample_id),
          identical(names(model)[-1], c("(Intercept)", "batchB", "conditiontreated")))
close_num(as.matrix(model[-1]), model.matrix(~ batch + condition, cd), "model matrix")
volcano <- read_artifact(root1, "volcano_data")
stopifnot(identical(volcano$gene_id, expected_genes),
          identical(volcano$gene_id[volcano$significant], expected_sig))
close_num(volcano$log2FoldChange, shrunk$log2FoldChange, "volcano x")
close_num(volcano$padj, raw$padj, "volcano FDR")
close_num(volcano$neg_log10_padj, -log10(pmax(raw$padj, .Machine$double.xmin)), "volcano y")
summary <- read_artifact(root1, "summary")
get_stat <- function(name) summary$value[match(name, summary$key)]
stopifnot(get_stat("genes_input") == "1000", get_stat("genes_retained") == "999",
          get_stat("genes_zero_total") == "1", get_stat("significant") == as.character(length(expected_sig)),
          get_stat("factor") == "condition", get_stat("numerator") == "treated",
          get_stat("denominator") == "control", get_stat("coefficient") == "condition_treated_vs_control")
pca_data <- read_artifact(root1, "pca_data")
stopifnot(identical(pca_data$sample_id, samples$sample_id))
# PC sign is arbitrary; compare after independent sign alignment.
for (k in 1:2) {
  got <- pca_data[[paste0("PC", k)]]
  sign <- if (sum(got * pca$x[, k]) < 0) -1 else 1
  close_num(sign * got, pca$x[, k], paste("PCA", k))
}
close_num(pca_data$PC1_variance, rep(pca$sdev[1]^2 / sum(pca$sdev^2), 8), "PC1 variance")
distance <- read_artifact(root1, "sample_distances")
stopifnot(identical(distance$sample_id, samples$sample_id), identical(names(distance)[-1], samples$sample_id))
close_num(as.matrix(distance[-1]), as.matrix(dist(t(v))), "sample distances")
heat <- read_artifact(root1, "heatmap_data")
heat_idx <- which(!is.na(raw$padj) & raw$padj < 0.05)
heat_idx <- head(heat_idx[order(raw$padj[heat_idx], expected_genes[heat_idx])], 50)
stopifnot(identical(heat$gene_id, expected_genes[heat_idx]), identical(names(heat)[-1], samples$sample_id))
close_num(as.matrix(heat[-1]), v[heat_idx, , drop = FALSE], "heatmap VST")
reverse <- read_artifact(file.path(out, "contrasts", "control_vs_treated"), "results")
# Check the reverse fit independently too. Relevelled optimizations differ slightly
# (observed max inversion error 2.58e-5); the 1e-4 inversion bound is separate
# from the strict per-model oracle and workers 1/2/4 numeric tolerance.
reverse_cd <- cd
reverse_cd$condition <- factor(samples$condition, levels = c("treated", "control"))
reverse_dds <- DESeq(DESeqDataSetFromMatrix(x, reverse_cd, ~ batch + condition), quiet = TRUE, parallel = FALSE)
reverse_raw <- results(reverse_dds, contrast = c("condition", "control", "treated"), alpha = 0.05)
for (name in names(as.data.frame(reverse_raw)))
  close_num(reverse[[name]], reverse_raw[[name]], paste("reverse oracle", name))
stopifnot(identical(reverse$gene_id, expected_genes), identical(is.na(reverse$padj), is.na(result$padj)),
          max(abs(reverse$log2FoldChange + result$log2FoldChange), na.rm = TRUE) < 1e-4,
          identical(reverse$gene_id[!is.na(reverse$padj) & reverse$padj < 0.05], expected_sig))
cat("PASS: independent DESeq2 oracle, identities, inversion, FDR, annotation and all plot semantics\n")
# Reordering, dropping, or changing worker-dependent values must be detected.
for (workers in c(2, 4)) {
  parallel_out <- file.path(work, paste0("workers-", workers))
  run <- invoke(parallel_out, workers)
  if (run$status != 0L) stop(paste(run$log, collapse = "\n"))
  for (id in c("treated_vs_control", "control_vs_treated")) {
    a <- file.path(out, "contrasts", id); b <- file.path(parallel_out, "contrasts", id)
    check_artifacts(b)
    for (name in c("results", "shrunken", "normalized_counts", "model_matrix", "volcano_data",
                   "pca_data", "sample_distances", "heatmap_data")) {
      reference <- read_artifact(a, name); actual <- read_artifact(b, name)
      stopifnot(identical(names(reference), names(actual)), identical(dim(reference), dim(actual)))
      for (col in names(reference)) {
        if (is.numeric(reference[[col]])) close_num(actual[[col]], reference[[col]], paste(workers, name, col))
        else stopifnot(identical(actual[[col]], reference[[col]]))
      }
    }
  }
  cat("PASS: workers", workers, "parity (1e-10 + 1e-7 * abs(reference))\n")
}
# Invalid input never publishes, and existing output is retained unchanged.
mutate_case <- function(name, mutate, pattern) {
  dir <- file.path(work, name); dir.create(dir)
  file.copy(list.files(fixture, pattern = "\\.tsv$", full.names = TRUE), dir)
  mutate(dir)
  target <- file.path(work, paste0(name, "-out"))
  run <- invoke(target, dir = dir)
  if (run$status == 0L || dir.exists(target) || !any(grepl(pattern, run$log, ignore.case = TRUE)))
    stop("negative case failed: ", name, "\n", paste(run$log, collapse = "\n"))
  cat("PASS: reject", name, "\n")
}
edit <- function(dir, file, fn) {
  path <- file.path(dir, paste0(file, ".tsv")); write_tsv(fn(read_tsv(path)), path)
}
# These rows were accepted by generic read.delim but violate the v1 dialect.
mutate_case("quoted-tab", function(d) edit(d, "annotation", function(x) {x$gene_symbol[1] <- "bad\tcell"; x}), "TSV|cell")
mutate_case("quoted-newline", function(d) edit(d, "annotation", function(x) {x$gene_symbol[1] <- "bad\ncell"; x}), "TSV|cell")
mutate_case("ragged-record", function(d) {
  p <- file.path(d, "annotation.tsv"); lines <- readLines(p); lines[2] <- paste0(lines[2], "\textra"); writeLines(lines, p)
}, "TSV|field|schema")
mutate_case("blank-record", function(d) cat("\n", file = file.path(d, "annotation.tsv"), append = TRUE), "TSV|record")
mutate_case("fitting-failure", function(d) edit(d, "counts", function(x) {
  for (i in seq_len(nrow(x) - 1L)) x[i, 2L + (i %% 8L)] <- 0L
  x
}), "every gene contains|geometric mean")
mutate_case("reordered-samples", function(d) edit(d, "samples", function(x) x[8:1, ]), "sample.*order")
mutate_case("missing-sample", function(d) edit(d, "samples", function(x) x[-1, ]), "sample.*order")
mutate_case("count-overflow", function(d) edit(d, "counts", function(x) {x[1, 2] <- "2147483648"; x}), "count")
mutate_case("count-fraction", function(d) edit(d, "counts", function(x) {x[1, 2] <- "1.5"; x}), "count")
mutate_case("count-negative", function(d) edit(d, "counts", function(x) {x[1, 2] <- "-1"; x}), "count")
mutate_case("undeclared-covariate", function(d) edit(d, "samples", function(x) {x$age <- 1:8; x}), "type declaration")
mutate_case("absent-covariate", function(d) edit(d, "samples", function(x) {x$batch <- NULL; x}), "absent")
for (value in c("NA", "NaN", "Inf", "1e999")) local({
  bad <- value
  mutate_case(paste0("numeric-", bad), function(d) {
    edit(d, "analysis", function(x) rbind(x, data.frame(key = "type.age", value = "numeric")))
    edit(d, "samples", function(x) {x$age <- as.character(1:8); x$age[1] <- bad; x})
  }, "finite decimal")
})
mutate_case("absent-level", function(d) edit(d, "contrasts", function(x) {x$numerator[1] <- "absent"; x}), "observed")
mutate_case("rank-deficient", function(d) edit(d, "samples", function(x) {x$batch <- rep(c("A", "B"), each = 4); x}), "rank")
mutate_case("zero-residual-df", function(d) {
  edit(d, "samples", function(x) {x$batch <- c("A", "B", "C", "D", "E", "F", "G", "G"); x})
}, "residual|rank")
mutate_case("unknown-annotation-gene", function(d) edit(d, "annotation", function(x) {x$gene_id[1] <- "unknown"; x}), "annotation")
mutate_case("duplicate-gene", function(d) edit(d, "counts", function(x) {x$gene_id[2] <- x$gene_id[1]; x}), "gene")
mutate_case("reference-conflict", function(d) edit(d, "analysis", function(x)
  rbind(x, data.frame(key = "reference.condition", value = "control"))), "reference")
mutate_case("unknown-setting", function(d) edit(d, "analysis", function(x)
  rbind(x, data.frame(key = "unexpected", value = "value"))), "unknown")
mutate_case("all-zero", function(d) edit(d, "counts", function(x) {x[-1] <- 0L; x}), "zero|retained")
mutate_case("unsafe-contrast-id", function(d) edit(d, "contrasts", function(x) {x$contrast_id[1] <- "../escape"; x}), "contrast.*ID")
# A full-rank 8x8 design (seven batches, split the shared batch across conditions).
mutate_case("saturated-model", function(d) edit(d, "samples", function(x) {
  x$batch <- c("A", "B", "C", "D", "D", "E", "F", "G"); x
}), "residual")
snapshot <- tools::md5sum(list.files(out, recursive = TRUE, full.names = TRUE))
run <- invoke(out)
stopifnot(run$status != 0L, identical(snapshot, tools::md5sum(names(snapshot))))
for (alias in c(fixture, file.path(fixture, "counts.tsv"), ".")) {
  run <- invoke(alias)
  stopifnot(run$status != 0L, any(grepl("alias|existing|exists", run$log)))
}
link <- file.path(work, "input-link")
stopifnot(file.symlink(normalizePath(fixture), link))
run <- invoke(link)
stopifnot(run$status != 0L)
cat("PASS: existing output and input/output aliases preserved\n")
# Identical distributions in each condition produce no significant genes.
no_signal <- file.path(work, "no-signal"); dir.create(no_signal)
invisible(file.copy(list.files(fixture, pattern = "\\.tsv$", full.names = TRUE), no_signal))
edit(no_signal, "counts", function(x) {x <- x[c(1:500, 1000), ]; x[6:9] <- x[2:5]; x})
empty_out <- file.path(work, "no-signal-out")
edit(no_signal, "analysis", function(x) {
  x$value[x$key == "type.batch"] <- "numeric"; x[x$key != "reference.batch", ]
})
edit(no_signal, "samples", function(x) {x$batch <- as.integer(x$batch == "B"); x})
edit(no_signal, "counts", function(x) {x$gene_id[1:2] <- c("#gene_one", 'gene"two'); x})
run <- invoke(empty_out, dir = no_signal, annotation = FALSE)
if (run$status != 0L) stop(paste(run$log, collapse = "\n"))
for (id in c("treated_vs_control", "control_vs_treated")) {
  root <- file.path(empty_out, "contrasts", id); check_artifacts(root)
  tab <- read_artifact(root, "results"); heat <- read_artifact(root, "heatmap_data")
  stopifnot(!any(!is.na(tab$padj) & tab$padj < 0.05), nrow(heat) == 0L,
            identical(names(heat), c("gene_id", samples$sample_id)),
            identical(tab$gene_id[1:2], c("#gene_one", 'gene"two')),
            !"gene_symbol" %in% names(tab))
}
stopifnot(length(list.files(work, pattern = "^\\.p3-stage", all.files = TRUE)) == 0L)
cat("PASS: zero-significant path; no leaked staged directories\n")
cat("PASS: primary integration cases; artifacts:", work, "\n")
# Exactly one significant row exercises the heatmap's no-row-clustering branch.
single <- file.path(work, "single-significant"); dir.create(single)
invisible(file.copy(list.files(fixture, pattern = "\\.tsv$", full.names = TRUE), single))
edit(single, "analysis", function(x) {x$value[x$key == "alpha"] <- "1e-46"; x})
edit(single, "contrasts", function(x) x[1, , drop = FALSE])
single_out <- file.path(work, "single-significant-out")
run <- invoke(single_out, dir = single)
if (run$status != 0L) stop(paste(run$log, collapse = "\n"))
single_root <- file.path(single_out, "contrasts", "treated_vs_control")
check_artifacts(single_root)
single_raw <- read_artifact(single_root, "results")
single_heat <- read_artifact(single_root, "heatmap_data")
stopifnot(sum(!is.na(single_raw$padj) & single_raw$padj < 1e-46) == 1L,
          identical(single_heat$gene_id, "gene_0170"))
close_num(as.matrix(single_heat[-1]), v["gene_0170", , drop = FALSE], "single-row heatmap VST")
cat("PASS: one-significant-gene heatmap; P3 suite complete\n")

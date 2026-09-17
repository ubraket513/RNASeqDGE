#!/usr/bin/env Rscript
# P3 offline trust boundary. All inputs and models are validated before staging.
fail <- function(...) stop(..., call. = FALSE)
require_true <- function(ok, ...) if (!isTRUE(ok)) fail(...)
is_symlink <- function(path) {
  target <- Sys.readlink(path)
  !is.na(target) && nzchar(target)
}
parse_cli <- function(args) {
  allowed <- c("counts", "samples", "analysis", "contrasts", "annotation", "out", "workers")
  if (identical(args, "--help")) {
    cat("Usage: Rscript --vanilla run_deg_analysis_offline.R --counts FILE --samples FILE\n",
      "  --analysis FILE --contrasts FILE [--annotation FILE] --out NEW_DIR [--workers N]\n",
      "OUT must not exist; its parent directory must exist. No network access is used.\n", sep = "")
    return(NULL)
  }
  require_true(length(args)%%2L == 0L, "options require explicit values")
  values <- list(workers = "1")
  seen <- character()
  for (i in seq.int(1L, length(args), by = 2L)) {
    name <- sub("^--", "", args[i])
    require_true(startsWith(args[i], "--") && name %in% allowed && !name %in% seen, "unknown or duplicate option: ",
      args[i])
    require_true(nzchar(args[i + 1L]), "empty option: ", name)
    values[[name]] <- args[i + 1L]
    seen <- c(seen, name)
  }
  require_true(all(c("counts", "samples", "analysis", "contrasts", "out") %in% seen), "required options: --counts --samples --analysis --contrasts --out")
  require_true(grepl("^[1-9][0-9]*$", values$workers) && is.finite(as.numeric(values$workers)) && as.numeric(values$workers) <=
    .Machine$integer.max, "workers must be a positive integer")
  values$workers <- as.integer(values$workers)
  values
}
read_tsv <- function(path, schema = NULL) {
  require_true(file_test("-f", path) && file.access(path, 4) == 0, "input must be a readable regular file: ",
    path)
  # Match the restricted v1 TSV dialect before read.delim performs decoding.
  # read.delim alone accepts embedded tabs/newlines and malformed quoting.
  lines <- withCallingHandlers(readLines(path, encoding = "UTF-8", warn = FALSE), warning = function(w) fail(path,
    ": ", conditionMessage(w)))
  require_true(length(lines) > 0L, "empty TSV: ", path)
  lines[1] <- sub("^﻿", "", lines[1])
  cell <- "(?:\"(?:[^\"\\t\\r\\n]|\"\")*\"|[^\"\\t\\r\\n]*)"
  syntax <- paste0("^", cell, "(\\t", cell, ")*$")
  require_true(all(validUTF8(lines)) && all(nzchar(lines)) && !any(startsWith(lines, "#")) && !any(grepl("﻿",
    lines, fixed = TRUE)) && all(grepl(syntax, lines, perl = TRUE)), "invalid TSV record or cell: ",
    path)
  connection <- textConnection(lines)
  on.exit(close(connection))
  widths <- count.fields(connection, sep = "\t", quote = "\"", comment.char = "", blank.lines.skip = FALSE)
  require_true(all(!is.na(widths)) && all(widths == widths[1]), "inconsistent TSV field count: ", path)
  close(connection)
  connection <- textConnection(lines)
  x <- withCallingHandlers(read.delim(connection, header = TRUE, sep = "\t", quote = "\"", comment.char = "",
    colClasses = "character", na.strings = character(), check.names = FALSE, stringsAsFactors = FALSE,
    fill = FALSE, blank.lines.skip = FALSE, row.names = NULL), warning = function(w) fail(path, ": ",
    conditionMessage(w)))
  require_true(length(names(x)) > 0L && all(nzchar(names(x))) && !anyDuplicated(names(x)), "empty or duplicate columns: ",
    path)
  require_true(all(vapply(x, function(v) all(!is.na(v) & validUTF8(v)), TRUE)), "missing cells or invalid UTF-8: ",
    path)
  if (!is.null(schema))
    require_true(setequal(names(x), schema), "invalid schema: ", path)
  x
}
valid_ids <- function(x) length(x) > 0L && all(grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", x)) && !anyDuplicated(x)
decimal <- function(x, context) {
  syntax <- "^[+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?$"
  values <- suppressWarnings(as.numeric(x))
  require_true(all(grepl(syntax, x) & is.finite(values)), context, " must be a finite decimal")
  values
}
validate_inputs <- function(opt) {
  # Existing directories, symlinks (including dangling links), files and ancestors
  # cannot be destinations. Exclusive publication also prevents accidental replacement.
  require_true(!file.exists(opt$out) && !dir.exists(opt$out) && !is_symlink(opt$out), "output exists or aliases an existing path: ",
    opt$out)
  parent <- normalizePath(dirname(opt$out), mustWork = TRUE)
  require_true(dir.exists(parent) && file.access(parent, 2) == 0, "output parent must be writable")
  out <- file.path(parent, basename(opt$out))
  samples <- read_tsv(opt$samples)
  require_true(all(c("sample_id", "condition") %in% names(samples)), "samples require sample_id and condition")
  require_true(valid_ids(samples$sample_id), opt$samples, ": invalid or duplicate sample ID")
  reserved_ids <- intersect(samples$sample_id, c("gene_id", "sample_id"))
  require_true(length(reserved_ids) == 0L, opt$samples, ": sample_id value(s) ", paste(reserved_ids,
    collapse = ", "), " are reserved for output identity columns")
  reserved_covariates <- intersect(names(samples), c("PC1", "PC2", "PC1_variance", "PC2_variance"))
  require_true(length(reserved_covariates) == 0L, opt$samples, ": covariate column(s) ", paste(reserved_covariates,
    collapse = ", "), " collide with reserved pca_data.tsv columns")
  settings <- read_tsv(opt$analysis, c("key", "value"))
  require_true(!anyDuplicated(settings$key) && all(nzchar(settings$key)), "duplicate or empty analysis key")
  settings <- setNames(settings$value, settings$key)
  required <- c("version", "design_terms", "alpha", "filter", "shrinkage")
  require_true(all(required %in% names(settings)), "required analysis setting absent")
  require_true(settings[["version"]] == "1" && settings[["filter"]] == "zero_total" && settings[["shrinkage"]] ==
    "apeglm", "analysis requires version=1, filter=zero_total, shrinkage=apeglm")
  alpha <- decimal(settings[["alpha"]], "alpha")
  require_true(alpha > 0 && alpha < 1, "alpha must lie strictly between 0 and 1")
  design_text <- settings[["design_terms"]]
  require_true(grepl("^[A-Za-z][A-Za-z0-9_]*(,[A-Za-z][A-Za-z0-9_]*)*$", design_text), "design_terms must be comma-separated declared identifiers")
  terms <- strsplit(design_text, ",", fixed = TRUE)[[1]]
  require_true(!anyDuplicated(terms) && all(terms %in% setdiff(names(samples), "sample_id")), "duplicate or absent design covariate")
  types <- c(condition = "categorical")
  refs <- character()
  for (key in setdiff(names(settings), required)) {
    if (startsWith(key, "type.")) {
      term <- substring(key, 6L)
      require_true(term %in% setdiff(names(samples), "sample_id"), "type declaration names absent covariate")
      require_true(settings[[key]] %in% c("categorical", "numeric"), "invalid covariate type")
      types[term] <- settings[[key]]
    } else if (startsWith(key, "reference.")) {
      term <- substring(key, 11L)
      require_true(term %in% setdiff(names(samples), "sample_id"), "reference names absent covariate")
      refs[term] <- settings[[key]]
    } else fail("unknown analysis setting: ", key)
  }
  require_true(types[["condition"]] == "categorical", "condition must be categorical")
  for (term in setdiff(names(samples), "sample_id")) {
    require_true(term %in% names(types), "covariate requires type declaration: ", term)
    if (types[[term]] == "numeric")
      samples[[term]] <- decimal(samples[[term]], paste("numeric covariate", term)) else require_true(all(nzchar(samples[[term]])), "empty categorical covariate: ", term)
  }
  for (term in names(refs)) require_true(types[[term]] == "categorical" && refs[[term]] %in% samples[[term]],
    "reference must be an observed categorical level: ", term)
  contrasts <- read_tsv(opt$contrasts, c("contrast_id", "factor", "numerator", "denominator"))
  require_true(valid_ids(contrasts$contrast_id), "invalid or duplicate contrast ID")
  for (i in seq_len(nrow(contrasts))) {
    row <- contrasts[i, ]
    term <- row$factor
    require_true(term %in% terms && types[[term]] == "categorical", "contrast factor must be categorical design term")
    require_true(row$numerator != row$denominator && all(c(row$numerator, row$denominator) %in% samples[[term]]),
      "contrast levels must be distinct observed levels")
    if (term %in% names(refs))
      require_true(refs[[term]] == row$denominator, "explicit reference disagrees with contrast denominator")
  }
  for (term in terms) if (types[[term]] == "categorical" && !term %in% contrasts$factor)
    require_true(term %in% names(refs), "categorical design term requires reference: ", term)
  input <- read_tsv(opt$counts)
  require_true(names(input)[1] == "gene_id" && identical(names(input)[-1], samples$sample_id), "count/sample identity and order must match exactly")
  require_true(nrow(input) > 0L && all(nzchar(input$gene_id)) && !anyDuplicated(input$gene_id), "empty or duplicate gene ID")
  count_text <- as.matrix(input[-1])
  values <- suppressWarnings(as.numeric(count_text))
  require_true(all(grepl("^[0-9]+$", count_text)) && all(is.finite(values)) && all(values <= .Machine$integer.max),
    "counts must be decimal integers in [0, 2147483647]")
  counts <- matrix(as.integer(values), nrow(input), ncol(input) - 1L, dimnames = list(input$gene_id,
    samples$sample_id))
  annotation <- NULL
  if (!is.null(opt$annotation)) {
    annotation <- read_tsv(opt$annotation, c("gene_id", "gene_symbol", "biotype", "chromosome"))
    require_true(all(nzchar(annotation$gene_id)) && !anyDuplicated(annotation$gene_id) && all(annotation$gene_id %in%
      input$gene_id), "annotation has duplicate, empty or unknown gene IDs")
  }
  keep <- rowSums(counts) > 0
  require_true(any(keep), "no genes retained after zero_total filtering")
  counts <- counts[keep, , drop = FALSE]
  # reformulate quotes reserved R identifiers, and receives only validated names.
  formula <- reformulate(vapply(terms, function(x) paste0("`", x, "`"), ""))
  models <- vector("list", nrow(contrasts))
  for (i in seq_len(nrow(contrasts))) {
    contrast <- contrasts[i, ]
    cd <- samples
    rownames(cd) <- cd$sample_id
    cd$sample_id <- NULL
    for (term in names(types)[types == "categorical"]) {
      observed <- sort(unique(cd[[term]]), method = "radix")
      reference <- if (term == contrast$factor)
        contrast$denominator else if (term %in% names(refs))
        refs[[term]] else observed[1]
      cd[[term]] <- factor(cd[[term]], levels = c(reference, setdiff(observed, reference)))
      if (term %in% terms) {
        require_true(nlevels(cd[[term]]) > 1L, "categorical model term needs at least two levels")
        require_true(!anyDuplicated(make.names(levels(cd[[term]]))), "categorical levels collide as coefficient names")
      }
    }
    mm <- model.matrix(formula, cd)
    require_true(all(is.finite(mm)) && qr(mm)$rank == ncol(mm), "model matrix is rank deficient")
    require_true(nrow(mm) > ncol(mm), "model requires positive residual degrees of freedom")
    require_true(!anyDuplicated(make.names(colnames(mm))), "model coefficient names are ambiguous")
    require_true(!"sample_id" %in% colnames(mm), opt$samples, ": model column 'sample_id' collides with reserved model_matrix.tsv identity column")
    coefficient <- make.names(paste0(contrast$factor, "_", contrast$numerator, "_vs_", contrast$denominator))
    models[[i]] <- list(data = cd, matrix = mm, coefficient = coefficient)
  }
  list(out = out, counts = counts, samples = samples, annotation = annotation, alpha = alpha, contrasts = contrasts,
    models = models, formula = formula, genes_input = nrow(input), genes_zero_total = sum(!keep))
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA",
  qmethod = "double")
write_matrix <- function(x, name, path) {
  ids <- data.frame(id = as.character(rownames(x)), stringsAsFactors = FALSE)
  names(ids) <- name
  write_tsv(cbind(ids, as.data.frame(x, check.names = FALSE)), path)
}
plot_png <- function(path, draw) {
  png(path, width = 1000, height = 800, res = 120)
  on.exit(dev.off())
  draw()
}
assert_same <- function(x, y, label) {
  require_true(identical(is.na(x), is.na(y)), "coefficient/contrast NA mismatch: ", label)
  keep <- !is.na(x)
  require_true(all(abs(x[keep] - y[keep]) <= 1e-10 + 1e-07 * abs(y[keep])), "coefficient does not match requested contrast: ",
    label)
}
peak_rss_kb <- function() {
  # Linux VmHWM is the parent R process high-water mark, not aggregate worker RSS.
  if (!file.exists("/proc/self/status"))
    return(NA_real_)
  line <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
  if (length(line) != 1L)
    return(NA_real_)
  as.numeric(strsplit(trimws(sub("^VmHWM:", "", line)), "[[:space:]]+")[[1]][1])
}
record_timing <- function(profile, stage, contrast_id, model_id, workers, status, wall_seconds = 0, cpu_seconds = 0) {
  profile$rows[[length(profile$rows) + 1L]] <- data.frame(stage = stage, contrast_id = contrast_id,
    model_id = model_id, status = status, workers = workers, wall_seconds = max(0, wall_seconds),
    cpu_seconds = max(0, cpu_seconds), self_peak_rss_kb = peak_rss_kb())
}
profile_call <- function(profile, stage, contrast_id, model_id, workers, expression) {
  start <- proc.time()
  value <- force(expression)
  elapsed <- proc.time() - start
  record_timing(profile, stage, contrast_id, model_id, workers, "executed", unname(elapsed[["elapsed"]]),
    sum(elapsed[c("user.self", "sys.self", "user.child", "sys.child")]))
  value
}
model_groups <- function(models) {
  # Counts and formula are one immutable validated object for this invocation.
  # Compare the full data frame too: identical matrices alone do not bind factor
  # levels or unused covariates needed by scientific output and coefficients.
  representatives <- integer()
  groups <- integer(length(models))
  for (i in seq_along(models)) {
    compatible <- vapply(representatives, function(j) identical(models[[i]]$matrix, models[[j]]$matrix) &&
      identical(models[[i]]$data, models[[j]]$data), TRUE)
    if (any(compatible))
      groups[i] <- which(compatible)[1] else {
      representatives <- c(representatives, i)
      groups[i] <- length(representatives)
    }
  }
  groups
}
run_contrast <- function(data, i, opt, stage, bp, fitted, profile, model_id, reused) {
  contrast <- data$contrasts[i, ]
  model <- data$models[[i]]
  dds <- fitted$dds
  timed <- function(name, expression) profile_call(profile, name, contrast$contrast_id, model_id, opt$workers,
    expression)
  cached <- function(name, expression) {
    if (exists(name, envir = fitted, inherits = FALSE)) {
      record_timing(profile, name, contrast$contrast_id, model_id, opt$workers, "reused")
      get(name, envir = fitted, inherits = FALSE)
    } else {
      value <- timed(name, expression)
      assign(name, value, envir = fitted)
      value
    }
  }
  if (reused)
    record_timing(profile, "fit", contrast$contrast_id, model_id, opt$workers, "reused")
  require_true(sum(DESeq2::resultsNames(dds) == model$coefficient) == 1L, "unsupported apeglm coefficient")
  raw <- timed("results", {
    result <- DESeq2::results(dds, contrast = c(contrast$factor, contrast$numerator, contrast$denominator),
      alpha = data$alpha, parallel = opt$workers > 1L, BPPARAM = bp)
    coef_result <- DESeq2::results(dds, name = model$coefficient, alpha = data$alpha)
    require_true(identical(rownames(result), rownames(coef_result)), "coefficient gene identity mismatch")
    for (column in names(result)) assert_same(result[[column]], coef_result[[column]], column)
    result
  })
  # Default apeglm nbinomCR is deterministic; seed also fixes future stochastic internals.
  set.seed(1)
  shrunk <- timed("shrinkage", DESeq2::lfcShrink(dds, coef = model$coefficient, res = raw, type = "apeglm",
    quiet = TRUE, parallel = opt$workers > 1L, BPPARAM = bp))
  require_true(identical(rownames(raw), rownames(shrunk)), "shrinkage gene identity mismatch")
  vst <- cached("vst", SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)))
  normalized <- cached("normalized_counts", DESeq2::counts(dds, normalized = TRUE))
  pc <- cached("pca", prcomp(t(vst), center = TRUE, scale. = FALSE))
  distance <- cached("sample_distances", as.matrix(dist(t(vst))))
  table_start <- proc.time()
  significant <- !is.na(raw$padj) & raw$padj < data$alpha
  root <- file.path(stage, "contrasts", contrast$contrast_id)
  require_true(dir.create(root, recursive = TRUE), "cannot create contrast output")
  result_table <- function(result) {
    tab <- data.frame(gene_id = rownames(result), as.data.frame(result), check.names = FALSE)
    if (!is.null(data$annotation)) {
      ann <- data$annotation[match(tab$gene_id, data$annotation$gene_id), c("gene_symbol", "biotype",
        "chromosome"), drop = FALSE]
      for (column in names(ann)) ann[[column]][!is.na(ann[[column]]) & ann[[column]] == ""] <- NA_character_
      rownames(ann) <- NULL
      tab <- cbind(tab, ann)
    }
    tab
  }
  write_tsv(result_table(raw), file.path(root, "results.tsv"))
  write_tsv(result_table(shrunk), file.path(root, "shrunken.tsv"))
  write_matrix(normalized, "gene_id", file.path(root, "normalized_counts.tsv"))
  write_matrix(model$matrix, "sample_id", file.path(root, "model_matrix.tsv"))
  stats <- c(contrast_id = contrast$contrast_id, factor = contrast$factor, numerator = contrast$numerator,
    denominator = contrast$denominator, coefficient = model$coefficient, alpha = data$alpha, genes_input = data$genes_input,
    genes_zero_total = data$genes_zero_total, genes_retained = nrow(raw), significant = sum(significant),
    padj_na = sum(is.na(raw$padj)), workers = opt$workers)
  write_tsv(data.frame(key = names(stats), value = unname(stats)), file.path(root, "summary.tsv"))
  volcano <- data.frame(gene_id = rownames(raw), log2FoldChange = shrunk$log2FoldChange, padj = raw$padj,
    neg_log10_padj = -log10(pmax(raw$padj, .Machine$double.xmin)), significant = significant)
  write_tsv(volcano, file.path(root, "volcano_data.tsv"))
  # Single retained gene still has meaningful PC1, with a zero PC2 placeholder.
  pc2 <- if (ncol(pc$x) >= 2L)
    pc$x[, 2] else rep(0, nrow(pc$x))
  total_variance <- sum(pc$sdev^2)
  variance <- if (total_variance > 0)
    pc$sdev^2/total_variance else rep(0, length(pc$sdev))
  pca <- data.frame(sample_id = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc2, PC1_variance = variance[1],
    PC2_variance = if (length(variance) >= 2L)
      variance[2] else 0, model$data, check.names = FALSE)
  write_tsv(pca, file.path(root, "pca_data.tsv"))
  write_matrix(distance, "sample_id", file.path(root, "sample_distances.tsv"))
  idx <- which(significant)
  idx <- head(idx[order(raw$padj[idx], rownames(raw)[idx], method = "radix")], 50L)
  heat <- vst[idx, , drop = FALSE]
  write_matrix(heat, "gene_id", file.path(root, "heatmap_data.tsv"))
  table_elapsed <- proc.time() - table_start
  record_timing(profile, "tables_and_plot_data", contrast$contrast_id, model_id, opt$workers, "executed",
    unname(table_elapsed[["elapsed"]]), sum(table_elapsed[c("user.self", "sys.self", "user.child",
      "sys.child")]))
  plot_start <- proc.time()
  title <- paste(contrast$numerator, "vs", contrast$denominator)
  plot_png(file.path(root, "MA_raw.png"), function() DESeq2::plotMA(raw, alpha = data$alpha, main = title))
  plot_png(file.path(root, "MA_shrunken.png"), function() DESeq2::plotMA(shrunk, alpha = data$alpha,
    main = title))
  plot_png(file.path(root, "dispersion.png"), function() DESeq2::plotDispEsts(dds))
  plot_png(file.path(root, "volcano.png"), function() {
    finite <- is.finite(volcano$log2FoldChange) & is.finite(volcano$neg_log10_padj)
    if (!any(finite)) {
      plot.new()
      text(0.5, 0.5, "No finite adjusted p-values")
    } else {
      plot(volcano$log2FoldChange[finite], volcano$neg_log10_padj[finite], pch = 16, cex = 0.5,
        col = ifelse(significant[finite], "firebrick", "grey60"), main = title, xlab = "apeglm log2 fold change",
        ylab = "-log10 adjusted p-value")
      abline(h = -log10(data$alpha), lty = 2)
    }
  })
  plot_png(file.path(root, "pca.png"), function() {
    plot(pca$PC1, pca$PC2, pch = 19, col = as.integer(model$data[[contrast$factor]]), xlab = sprintf("PC1 (%.1f%%)",
      100 * pca$PC1_variance[1]), ylab = sprintf("PC2 (%.1f%%)", 100 * pca$PC2_variance[1]), main = "VST sample PCA")
    text(pca$PC1, pca$PC2, labels = pca$sample_id, pos = 3, cex = 0.7)
  })
  plot_png(file.path(root, "sample_distances.png"), function() pheatmap::pheatmap(distance, main = "Euclidean distances of VST samples",
    silent = FALSE))
  plot_png(file.path(root, "heatmap.png"), function() {
    if (nrow(heat) == 0L) {
      plot.new()
      text(0.5, 0.5, "No significant genes (adjusted p-value < alpha)")
    } else pheatmap::pheatmap(heat, scale = "none", cluster_rows = nrow(heat) > 1L, main = "Significant genes: VST expression",
      fontsize_row = 6, silent = FALSE)
  })
  writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
  plot_elapsed <- proc.time() - plot_start
  record_timing(profile, "plots_and_session", contrast$contrast_id, model_id, opt$workers, "executed",
    unname(plot_elapsed[["elapsed"]]), sum(plot_elapsed[c("user.self", "sys.self", "user.child",
      "sys.child")]))
}
main <- function() {
  total_start <- proc.time()
  opt <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (is.null(opt))
    return(invisible(NULL))
  profile <- new.env(parent = emptyenv())
  profile$rows <- list()
  data <- profile_call(profile, "validation", "", "", opt$workers, validate_inputs(opt))
  profile_call(profile, "package_load", "", "", opt$workers, {
    for (package in c("DESeq2", "apeglm", "BiocParallel", "pheatmap")) require_true(requireNamespace(package,
      quietly = TRUE), "required package unavailable: ", package)
  })
  bp <- if (opt$workers == 1L)
    BiocParallel::SerialParam(RNGseed = 1L) else BiocParallel::MulticoreParam(workers = opt$workers, RNGseed = 1L, progressbar = FALSE)
  stage <- tempfile(".p3-stage-", tmpdir = dirname(data$out))
  require_true(dir.create(stage, mode = "0700"), "cannot create exclusive staging directory")
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  groups <- model_groups(data$models)
  for (group in unique(groups)) {
    indices <- which(groups == group)
    first <- indices[1]
    model_id <- paste0("model_", group)
    # Retain only one compatible fitted model at a time, bounding cache memory.
    fitted <- new.env(parent = emptyenv())
    fitted$dds <- profile_call(profile, "fit", data$contrasts$contrast_id[first], model_id, opt$workers,
      {
        dds <- DESeq2::DESeqDataSetFromMatrix(data$counts, data$models[[first]]$data, data$formula)
        DESeq2::DESeq(dds, quiet = TRUE, parallel = opt$workers > 1L, BPPARAM = bp)
      })
    for (i in indices) run_contrast(data, i, opt, stage, bp, fitted, profile, model_id, i != first)
    rm(fitted, dds)
    invisible(gc())
  }
  total_elapsed <- proc.time() - total_start
  record_timing(profile, "total", "", "", opt$workers, "executed", unname(total_elapsed[["elapsed"]]),
    sum(total_elapsed[c("user.self", "sys.self", "user.child", "sys.child")]))
  write_tsv(do.call(rbind, profile$rows), file.path(stage, "timings.tsv"))
  require_true(!file.exists(data$out) && !dir.exists(data$out) && !is_symlink(data$out), "output appeared during analysis; refusing replacement")
  require_true(file.rename(stage, data$out), "cannot publish staged output")
  cat("Published offline analysis:", data$out, "\n")
}
tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(save = "no", status = 1L)
})

#!/usr/bin/env Rscript
# Compare scientific TSVs, excluding additive performance metadata and worker label.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, all(dir.exists(args)))
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
files <- list.files(file.path(args[1], "contrasts"), pattern = "\\.tsv$", recursive = TRUE)
stopifnot(length(files) > 0L,
  identical(files, list.files(file.path(args[2], "contrasts"), pattern = "\\.tsv$", recursive = TRUE)))
for (file in files) {
  reference <- read_tsv(file.path(args[1], "contrasts", file))
  actual <- read_tsv(file.path(args[2], "contrasts", file))
  stopifnot(identical(names(actual), names(reference)), identical(dim(actual), dim(reference)))
  if (basename(file) == "summary.tsv") {
    reference <- reference[reference$key != "workers", , drop = FALSE]
    actual <- actual[actual$key != "workers", , drop = FALSE]
  }
  for (name in names(reference)) {
    a <- actual[[name]]; b <- reference[[name]]
    if (!identical(is.na(a), is.na(b))) stop("NA mask mismatch: ", file, " / ", name)
    if (is.numeric(b)) {
      keep <- !is.na(b)
      if (any(abs(a[keep] - b[keep]) > 1e-10 + 1e-7 * abs(b[keep])))
        stop("numeric mismatch: ", file, " / ", name)
    } else if (!identical(a, b)) stop("identity mismatch: ", file, " / ", name)
  }
}
cat("PASS", length(files), "scientific tables; exact identities/NA masks; tolerance 1e-10 + 1e-7 * abs(reference)\n")

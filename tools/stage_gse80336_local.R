#!/usr/bin/env Rscript
# Restricted-reference smoke preparation; never a whole-genome method benchmark.
source(file.path(
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
  if (basename(dirname(sub(
    "^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]
  ))) == "integration")
    "../../tools/lib_helpers.R" else "lib_helpers.R"
))
ACCESSION <- "NC_000022.11"
ASSEMBLY <- "GCF_000001405.40_GRCh38.p14"
validate_read_stage <- function(stage) {
  assert(is.character(stage$sampling) && length(stage$sampling) == 1L && startsWith(stage$sampling,
    "first N complete FASTQ records in archive order;"), "read-stage must contain contiguous archive prefixes, not subsampled selections")
  assert(length(stage$reads) > 0L && all(vapply(stage$reads, function(read) is.null(read$selection),
    TRUE)), "read-stage contains a non-prefix selection")
  invisible(TRUE)
}
copy_fastq_prefix <- function(source, destination, records, check_deadline = function() {
}) {
  for (record in seq_len(records)) {
    check_deadline()
    lines <- withCallingHandlers(readLines(source, n = 4L, warn = TRUE), warning = function(w) stop("incomplete FASTQ line"))
    assert(length(lines) == 4L, paste("incomplete FASTQ record", record))
    assert(startsWith(lines[1], "@") && startsWith(lines[3], "+"), "invalid FASTQ delimiters")
    assert(nchar(lines[2], type = "bytes") > 0L && nchar(lines[2], type = "bytes") == nchar(lines[4],
      type = "bytes"), "FASTQ sequence/quality length mismatch")
    assert(grepl("^[ACGTNacgtn]+$", lines[2]), "invalid FASTQ sequence")
    qualities <- as.integer(charToRaw(lines[4]))
    assert(all(qualities >= 33L & qualities <= 126L), "invalid FASTQ quality")
    writeLines(lines, destination, sep = "\n", useBytes = TRUE)
  }
}
main <- function() {
  args <- parse_cli(list(output = file.path(ROOT, "tests/output/p6-local-data"), samples = file.path(ROOT,
    "tests/output/p6-inputs/subset6/samples.tsv"), mapping = file.path(ROOT, "tests/output/p6-inputs/sample_run_mapping.tsv"),
    sources = file.path(ROOT, "config/gse80336-sources.json"), records = "50000", deadline = NULL,
    timeout = "30", `ena-scheme` = "https", `reuse-reference-downloads` = FALSE, `read-stage` = NULL),
    required = "deadline", flags = "reuse-reference-downloads")
  records <- as.integer(args$records)
  timeout <- as.integer(args$timeout)
  assert(!is.na(records) && records > 0L && !is.na(timeout) && timeout > 0L, "records/timeout must be positive")
  assert(args[["ena-scheme"]] %in% c("https", "http", "ftp"), "invalid ENA transport")
  deadline <- as.POSIXct(args$deadline, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  assert(!is.na(deadline), "deadline must be UTC ISO8601")
  start <- Sys.time()
  check_deadline <- function() assert(Sys.time() < deadline, "staging deadline exceeded")
  remaining <- function() {
    check_deadline()
    max(1L, as.integer(difftime(deadline, Sys.time(), units = "secs")))
  }
  curl_args <- function(url) c("curl", "--fail", "--location", "--silent", "--show-error", "--connect-timeout",
    timeout, "--max-time", remaining(), url)
  out <- absolute(args$output)
  downloads <- file.path(out, "downloads")
  reads_dir <- file.path(out, "reads")
  dir.create(downloads, recursive = TRUE, showWarnings = FALSE)
  dir.create(reads_dir, showWarnings = FALSE)
  assert(!file.exists(file.path(out, "STAGED")), "output already staged")
  sources <- read_json(args$sources)
  assert(sources$reference$assembly == ASSEMBLY, "assembly mismatch")
  samples <- table_read(args$samples)
  mapping <- table_read(args$mapping)
  assert(nrow(samples) == 6L && !anyDuplicated(samples$sample_id), "six distinct samples required")
  assert(!anyDuplicated(mapping$sample_id), "duplicate sample/run mapping")
  selected <- mapping[match(samples$sample_id, mapping$sample_id), , drop = FALSE]
  assert(!anyNA(selected) && all(selected$layout == "single" & selected$strandedness == "reverse"),
    "expected mapped reverse single-end samples")
  provenance <- list(study = "GSE80336", assembly = ASSEMBLY, accession = ACCESSION, started_utc = format(start,
    "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), deadline_utc = args$deadline, records_per_run = records, scope = "restricted-reference smoke check, not a full method benchmark",
    sampling = "first N complete FASTQ records in archive order; nonrandom prefix sampling may be biased",
    source_reads_full_md5_verified = FALSE, references = list(), reads = list())
  for (name in names(sources$reference$files)) {
    check_deadline()
    info <- sources$reference$files[[name]]
    path <- file.path(downloads, name)
    if (!isTRUE(args[["reuse-reference-downloads"]])) {
      partial <- paste0(path, ".partial")
      run(c(curl_args(info$url), "-o", partial))
      assert(file.rename(partial, path), "reference publish failed")
    }
    actual <- unname(tools::md5sum(path))
    assert(!is.na(actual) && actual == info$md5, paste("full reference MD5 mismatch:", name))
    provenance$references[[length(provenance$references) + 1L]] <- list(file = name, url = info$url,
      md5 = actual, md5_verified = TRUE, compressed_bytes = file.info(path)$size)
  }
  old <- if (!is.null(args[["read-stage"]]))
    read_json(file.path(args[["read-stage"]], "provenance.json")) else NULL
  if (!is.null(old))
    validate_read_stage(old)
  for (row in rows_list(selected)) {
    check_deadline()
    url <- paste0(args[["ena-scheme"]], "://", row$fastq_ftp)
    assert(!grepl(";", url, fixed = TRUE), "expected one single-end source")
    path <- file.path(reads_dir, paste0(row$run_id, ".fastq"))
    partial <- paste0(path, ".partial")
    begin <- Sys.time()
    if (is.null(old)) {
      # Intentional bounded stream: upstream curl/gzip see SIGPIPE when N records
      # have been consumed. No full source MD5 or gzip-CRC claim is made.
      command <- paste(paste(shQuote(curl_args(url)), collapse = " "), "| gzip -dc")
      incoming <- pipe(command, "r")
    } else {
      matches <- Filter(function(x) x$run_id == row$run_id, old$reads)
      assert(length(matches) == 1L, "missing source prefix")
      prior <- matches[[1]]
      original <- file.path(args[["read-stage"]], prior$path)
      assert(sha256(original) == prior$sha256 && prior$prefix_records >= records, "source prefix hash/count mismatch")
      incoming <- file(original, "r")
      url <- prior$url
    }
    outgoing <- file(partial, "w")
    tryCatch(copy_fastq_prefix(incoming, outgoing, records, check_deadline), finally = {
      close(incoming)
      close(outgoing)
    })
    assert(file.rename(partial, path), "FASTQ prefix publish failed")
    result <- list(sample_id = row$sample_id, run_id = row$run_id, url = url, source_full_md5 = row$fastq_md5,
      source_full_md5_verified = FALSE, source_full_gzip_crc_verified = FALSE, source_total_records = as.numeric(row$read_count),
      prefix_records = records, path = file.path("reads", basename(path)), sha256 = sha256(path),
      bytes = file.info(path)$size, seconds = as.numeric(difftime(Sys.time(), begin, units = "secs")))
    if (!is.null(old)) {
      result$local_verified_parent_stage <- absolute(args[["read-stage"]])
      result$source_prefix_sha256 <- prior$sha256
      result$source_prefix_records <- prior$prefix_records
    }
    write_json(result, file.path(reads_dir, paste0(row$run_id, ".provenance.json")))
    provenance$reads[[length(provenance$reads) + 1L]] <- result
    cat("STAGED", row$run_id, records, "records\n")
  }
  # Stream 10k lines at a time; never materialize the full human reference in RAM.
  incoming <- gzfile(file.path(downloads, paste0(ASSEMBLY, "_genomic.fna.gz")), "rt")
  target <- file(file.path(out, "reference.fa"), "w")
  bases <- 0
  found <- 0L
  active <- FALSE
  tryCatch(repeat {
    lines <- readLines(incoming, n = 10000L)
    if (!length(lines))
      break
    check_deadline()
    first <- 1L
    for (header in which(startsWith(lines, ">"))) {
      if (active && header > first) {
        chunk <- lines[seq.int(first, header - 1L)]
        bases <- bases + sum(nchar(trimws(chunk)))
        writeLines(chunk, target)
      }
      active <- strsplit(substring(lines[header], 2), " ")[[1]][1] == ACCESSION
      if (active) {
        found <- found + 1L
        writeLines(lines[header], target)
      }
      first <- header + 1L
    }
    if (active && first <= length(lines)) {
      chunk <- lines[seq.int(first, length(lines))]
      bases <- bases + sum(nchar(trimws(chunk)))
      writeLines(chunk, target)
    }
  }, finally = {
    close(incoming)
    close(target)
  })
  assert(found == 1L && bases == 50818468, "NC_000022.11 FASTA accession/length mismatch")
  incoming <- gzfile(file.path(downloads, paste0(ASSEMBLY, "_genomic.gtf.gz")), "rt")
  target <- file(file.path(out, "reference.gtf"), "w")
  genes <- list()
  exon_rows <- 0L
  tryCatch(repeat {
    lines <- readLines(incoming, n = 10000L)
    if (!length(lines))
      break
    check_deadline()
    lines <- lines[startsWith(lines, "#") | startsWith(lines, paste0(ACCESSION, "\t"))]
    for (line in lines) {
      if (startsWith(line, "#")) {
        writeLines(line, target)
        next
      }
      fields <- strsplit(line, "\t")[[1]]
      assert(length(fields) == 9L, "invalid GTF row")
      if (fields[1] != ACCESSION)
        next
      assert(as.numeric(fields[4]) >= 1 && as.numeric(fields[5]) <= bases, "GTF coordinates outside accession")
      writeLines(line, target)
      if (fields[3] == "exon") {
        attribute <- function(key) {
          pattern <- paste0("(^|; )", key, " \"([^\"]*)\";")
          m <- regexec(pattern, fields[9])
          match <- regmatches(fields[9], m)[[1]]
          if (length(match))
          match[3] else ""
        }
        gene <- attribute("gene_id")
        assert(nzchar(gene), "missing gene_id")
        genes[[gene]] <- data.frame(gene_id = gene, gene_symbol = attribute("gene"), biotype = attribute("gene_biotype"),
          chromosome = ACCESSION)
        exon_rows <- exon_rows + 1L
      }
    }
  }, finally = {
    close(incoming)
    close(target)
  })
  assert(exon_rows > 0L && length(genes) > 0L, "no annotated exons")
  copy_files(args$samples, file.path(out, "samples.tsv"))
  table_write(data.frame(run_id = selected$run_id, sample_id = selected$sample_id, fastq_1 = paste0("reads/",
    selected$run_id, ".fastq"), fastq_2 = "", layout = "single", strandedness = "reverse"), file.path(out,
    "runs.tsv"))
  names <- c("reference.fa", "reference.gtf")
  suffix <- c("_genomic.fna.gz", "_genomic.gtf.gz")
  table_write(data.frame(role = c("genome", "annotation"), path = names, sha256 = vapply(file.path(out,
    names), sha256, ""), source = vapply(suffix, function(s) paste0(sources$reference$files[[paste0(ASSEMBLY,
    s)]]$url, "#", ACCESSION), ""), release = ASSEMBLY), file.path(out, "references.tsv"))
  kv_write(list(version = "1", design_terms = "condition", alpha = "0.05", filter = "zero_total", shrinkage = "apeglm",
    reference.condition = "control", type.age = "numeric", type.sex = "categorical", type.pmi = "numeric",
    type.rin = "numeric"), file.path(out, "analysis.tsv"))
  table_write(data.frame(contrast_id = "bipolar_vs_control", factor = "condition", numerator = "bipolar",
    denominator = "control"), file.path(out, "contrasts.tsv"))
  table_write(data.frame(gene_id = sort(names(genes))), file.path(out, "genes.tsv"))
  table_write(do.call(rbind, genes[sort(names(genes))]), file.path(out, "annotation.tsv"))
  provenance$reference_bases <- bases
  provenance$exon_rows <- exon_rows
  provenance$genes <- length(genes)
  provenance$elapsed_seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  provenance$completed_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  write_json(provenance, file.path(out, "provenance.json"))
  writeLines("restricted-reference smoke data prepared; not full-method validation", file.path(out,
    "STAGED"))
  cat("STAGED", out, length(genes), "genes\n")
}
if (sys.nframe() == 0L) main()

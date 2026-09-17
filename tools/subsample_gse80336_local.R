#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "lib_helpers.R"))
args <- parse_cli(list(source = NULL, target = NULL), required = c("source", "target"), positional = c("source",
  "target"))
source <- absolute(args$source)
target <- new_destination(args$target)
dir.create(target)
dir.create(file.path(target, "reads"))
provenance <- read_json(file.path(source, "provenance.json"))
for (i in seq_along(provenance$reads)) {
  read <- provenance$reads[[i]]
  original <- file.path(source, read$path)
  stopifnot(sha256(original) == read$sha256)
  incoming <- file(original, "r")
  outgoing <- file(file.path(target, read$path), "w")
  count <- 0L
  tryCatch({
    for (index in 0:49999) {
      record <- readLines(incoming, n = 4L)
      assert(length(record) == 4L, "short source prefix")
      if (index%%5L == 0L) {
        writeLines(record, outgoing)
        count <- count + 1L
      }
    }
    assert(!length(readLines(incoming, n = 1L)) && count == 10000L, "unexpected source record count")
  }, finally = {
    close(incoming)
    close(outgoing)
  })
  read$source_prefix_records <- 50000
  read$selected_records <- count
  read$prefix_records <- count
  read$selection <- "zero-based records 0,5,...,49995 of 50k archive prefix"
  read$source_prefix_sha256 <- read$sha256
  read$sha256 <- sha256(file.path(target, read$path))
  read$bytes <- file.info(file.path(target, read$path))$size
  provenance$reads[[i]] <- read
}
copy_files(file.path(source, c("samples.tsv", "runs.tsv", "analysis.tsv", "contrasts.tsv", "genes.tsv",
  "annotation.tsv")), target)
refs <- table_read(file.path(source, "references.tsv"))
refs$path <- file.path(source, refs$path)
table_write(refs, file.path(target, "references.tsv"))
provenance$records_per_run <- 10000
provenance$sampling <- "every fifth record of first50k archive prefix; nonrandom"
provenance$parent_stage <- source
provenance$reason <- "50k STAR pilot exceeded budget for six workflows"
write_json(provenance, file.path(target, "provenance.json"))
writeLines("six verified 10k read selections for reduced local check", file.path(target, "STAGED"))
cat("STAGED six samples x 10,000 reads\n")

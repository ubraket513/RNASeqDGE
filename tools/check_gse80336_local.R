#!/usr/bin/env Rscript
# Explicit reduced real-read benchmark: serial jobs and a bounded wall deadline.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "lib_helpers.R"))
main <- function() {
  args <- parse_cli(list(inputs = NULL, out = NULL, seconds = "1800", replicates = "3", `index-cache` = NULL,
    `bin-dir` = file.path(ROOT, ".deps/runtime-tools/bin"), rscript = file.path(ROOT, ".deps/runtime-r/bin/Rscript"),
    `tool-lock` = file.path(ROOT, "config/alignment-linux-64.explicit.txt"), `r-lock` = file.path(ROOT,
      "config/runtime-r-linux-64.explicit.txt")), required = c("inputs", "out"))
  seconds <- as.integer(args$seconds)
  replicates <- as.integer(args$replicates)
  assert(!is.na(seconds) && seconds >= 1L && seconds <= 1800L && !is.na(replicates) && replicates >=
    1L && replicates <= 3L, "seconds must be 1..1800 and replicates 1..3")
  inputs <- absolute(args$inputs)
  out <- new_destination(args$out)
  dir.create(out)
  deadline <- proc.time()[["elapsed"]] + seconds
  native <- file.path(ROOT, "build/rnaseq")
  native_hash <- sha256(native)
  records <- list()
  record_progress <- function() write_json(list(scope = "6 samples, sampled archive-prefix reads, chromosome22 only; not whole-genome inference",
    native_sha256 = native_hash, threads = 2, workers = 1, requested_replicates = replicates, complete = length(records) ==
      2L * replicates && all(vapply(records, function(r) r$status == 0L, TRUE)), runs = records),
    file.path(out, "runs.json"))
  record_progress()
  with_environment(c(clean_environment(), list(PATH = "/usr/bin:/bin")), {
    for (repetition in seq_len(replicates)) for (backend in if (repetition%%2L)
      c("hisat2", "star") else c("star", "hisat2")) {
      remaining <- as.integer(deadline - proc.time()[["elapsed"]])
      assert(remaining > 0L, "local-check deadline reached; partial results retained")
      assert(sha256(native) == native_hash, "native executable changed during check")
      name <- paste0(backend, "-", repetition)
      keys <- c("samples", "runs", "references", "analysis", "contrasts", "genes")
      settings <- setNames(as.list(file.path(inputs, paste0(keys, ".tsv"))), keys)
      settings <- c(settings, list(bin_dir = absolute(args[["bin-dir"]]), rscript = absolute(args$rscript),
        r_script = file.path(ROOT, "run_deg_analysis_offline.R"), tool_lock = absolute(args[["tool-lock"]]),
        r_lock = absolute(args[["r-lock"]]), run_dir = file.path(out, name), index_cache = if (is.null(args[["index-cache"]])) file.path(out,
          "index-cache") else absolute(args[["index-cache"]]), backend = backend, threads = 2,
        workers = 1, star_sa_bases = 11, star_chr_bits = 16))
      config <- file.path(out, paste0(name, ".tsv"))
      kv_write(settings, config)
      command <- c("/usr/bin/timeout", "--signal=TERM", "--kill-after=15s", paste0(remaining, "s"),
        "/usr/bin/time", "-v", "-o", file.path(out, paste0(name, ".time.txt")), native, "workflow-local",
        config)
      cat("START", name, remaining, "seconds remain\n")
      start <- proc.time()[["elapsed"]]
      status <- run(c(file.path(ROOT, "tools/run_guarded.sh"), command), log = file.path(out, paste0(name,
        ".log")), check = FALSE)$status
      records[[length(records) + 1L]] <- list(backend = backend, `repeat` = repetition, status = status,
        wall_seconds = proc.time()[["elapsed"]] - start, command = as.list(command), index_condition = if (!is.null(args[["index-cache"]])) "reuse from pilot" else if (repetition ==
          1L) "build" else "reuse", cache_condition = "uncontrolled OS cache; no cache eviction requested")
      record_progress()
      cat("END", name, "exit=", status, "\n")
      assert(status == 0L, paste(name, "failed; logs retained"))
    }
  })
  cat("All requested reduced workflows completed; scientific comparison remains required.\n")
}
if (sys.nframe() == 0L) main()

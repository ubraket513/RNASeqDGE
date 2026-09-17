# Shared base-R/jsonlite support for repository preparation and independent checks.
script_path <- function() sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
repo_root <- function() normalizePath(file.path(dirname(script_path()), if (basename(dirname(script_path())) ==
  "integration") "../.." else ".."), mustWork = TRUE)
ROOT <- repo_root()
assert <- function(ok, message = "check failed") if (!isTRUE(ok)) stop(message, call. = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
write_json <- function(value, path) jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE,
  null = "null", na = "null", digits = NA)
table_read <- function(path) read.delim(path, colClasses = "character", check.names = FALSE, quote = "\"",
  comment.char = "", na.strings = NULL, stringsAsFactors = FALSE)
table_write <- function(value, path) write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE,
  col.names = TRUE, na = "")
rows_list <- function(value) lapply(seq_len(nrow(value)), function(i) as.list(value[i, , drop = FALSE]))
kv_read <- function(path) {
  x <- table_read(path)
  setNames(as.list(x$value), x$key)
}
kv_write <- function(value, path) table_write(data.frame(key = names(value), value = unlist(value), check.names = FALSE),
  path)
absolute <- function(path) normalizePath(if (startsWith(path, "/")) path else file.path(getwd(), path),
  mustWork = FALSE)
exists_any <- function(path) {
  link <- Sys.readlink(path)
  file.exists(path) || (!is.na(link) && nzchar(link))
}
new_destination <- function(path) {
  path <- absolute(path)
  assert(!exists_any(path), paste("destination already exists:", path))
  assert(dir.exists(dirname(path)), paste("destination parent must exist:", dirname(path)))
  path
}
new_directory <- function(prefix, parent = tempdir()) {
  path <- tempfile(prefix, tmpdir = parent)
  assert(dir.create(path), "cannot create staging directory")
  path
}
sha256 <- function(path) {
  assert(file.exists(path) && !dir.exists(path), paste("missing regular file:", path))
  value <- system2("sha256sum", c("--", shQuote(path)), stdout = TRUE)
  assert(is.null(attr(value, "status")) && length(value) == 1L, "sha256sum failed")
  hash <- sub(" .*$", "", value)
  assert(grepl("^[0-9a-f]{64}$", hash), "invalid sha256 output")
  hash
}
run <- function(argv, log = NULL, timeout = 0, check = TRUE) {
  assert(length(argv) > 0L, "empty command")
  result <- suppressWarnings(system2(argv[1], shQuote(argv[-1]), stdout = if (is.null(log))
    TRUE else log, stderr = if (is.null(log))
    TRUE else log, timeout = timeout))
  status <- if (is.null(log)) {
    x <- attr(result, "status")
    if (is.null(x))
      0L else x
  } else result
  if (check)
    assert(status == 0L, paste("command failed with exit", status, ":", argv[1], if (!is.null(log))
      paste("log:", log)))
  list(status = status, output = if (is.null(log)) result else character())
}
parse_cli <- function(defaults = list(), required = character(), flags = character(), args = commandArgs(TRUE),
  positional = character()) {
  result <- defaults
  used <- character()
  p <- 1L
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--help") {
      cat("Options:", paste(paste0("--", names(defaults)), collapse = " "), "\nRequired:", paste(required,
        collapse = " "), "\n")
      quit(status = 0)
    }
    if (!startsWith(key, "--")) {
      assert(p <= length(positional), "unexpected positional argument")
      name <- positional[p]
      p <- p + 1L
      value <- key
    } else {
      name <- sub("^--", "", key)
      assert(name %in% names(defaults), paste("unknown option", key))
      if (name %in% flags)
        value <- TRUE else {
        i <- i + 1L
        assert(i <= length(args), paste("missing option value", key))
        value <- args[i]
      }
    }
    assert(!name %in% used, paste("duplicate option", name))
    used <- c(used, name)
    result[[name]] <- value
    i <- i + 1L
  }
  for (name in required) assert(!is.null(result[[name]]) && nzchar(as.character(result[[name]])), paste("required:",
    name))
  result
}
clean_environment <- function() list(LD_LIBRARY_PATH = NULL, LD_PRELOAD = NULL, R_HOME = NULL, R_LIBS = NULL,
  R_LIBS_USER = "/dev/null", R_LIBS_SITE = "/dev/null", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
with_environment <- function(values, code) {
  old <- Sys.getenv(names(values), unset = NA_character_)
  on.exit({
    for (name in names(old)) if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv,
      setNames(list(old[[name]]), name))
  })
  for (name in names(values)) if (is.null(values[[name]]))
    Sys.unsetenv(name) else do.call(Sys.setenv, setNames(list(values[[name]]), name))
  force(code)
}
copy_files <- function(files, destination) {
  assert(all(file.copy(files, destination)), "file copy failed")
}
reverse_complement <- function(value) paste(rev(strsplit(chartr("ACGT", "TGCA", value), "")[[1]]), collapse = "")

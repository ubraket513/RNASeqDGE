#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "lib_helpers.R"))
args <- parse_cli(list(prefix = NULL), required = "prefix", positional = "prefix")
prefix <- absolute(args$prefix)
fields <- c("name", "version", "build", "subdir", "url", "sha256", "md5", "license", "depends")
packages <- lapply(sort(list.files(file.path(prefix, "conda-meta"), "[.]json$", full.names = TRUE)),
  function(path) {
    r <- read_json(path)
    setNames(lapply(fields, function(k) r[[k]]), fields)
  })
archives <- lapply(sort(list.files(file.path(prefix, "share"), pattern = "[.]tar[.]gz", recursive = TRUE,
  full.names = TRUE)), function(path) list(path = substring(path, nchar(prefix) + 2L), sha256 = sha256(path)))
sources <- list()
catalog <- file.path(prefix, "share/bioconductor-data-packages/dataURLs.json")
if (file.exists(catalog)) {
  catalog <- read_json(catalog)
  for (record in packages) {
    key <- paste0(sub("^bioconductor-", "", record$name), "-", record$version)
    if (key %in% names(catalog))
      sources[[key]] <- catalog[[key]]
  }
}
cat(jsonlite::toJSON(list(platform = "linux-64", packages = packages, post_link_archives = archives,
  post_link_sources = sources), auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n")

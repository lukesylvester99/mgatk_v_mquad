#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(epiAneufinder)
})

# --------------- CLI ---------------
parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE); kv <- list(); i <- 1
  while (i <= length(args)) {
    tok <- args[[i]]
    if (!startsWith(tok, "--")) { i <- i + 1; next }
    tok <- sub("^--", "", tok)
    if (grepl("=", tok, fixed = TRUE)) {
      parts <- strsplit(tok, "=", fixed = TRUE)[[1]]
      kv[[parts[[1]]]] <- paste(parts[-1], collapse = "="); i <- i + 1
    } else {
      if (i + 1 <= length(args) && !startsWith(args[[i + 1]], "--")) {
        kv[[tok]] <- args[[i + 1]]; i <- i + 2
      } else {
        kv[[tok]] <- TRUE; i <- i + 1
      }
    }
  }
  kv
}
opt <- parse_cli()
must <- function(x) { if (is.null(opt[[x]]) || opt[[x]] == "") stop(sprintf("Missing --%s", x)); opt[[x]] }

sample   <- must("sample")
frags_in <- must("fragments")      # cleaned fragments .tsv(.gz) OK
outdir   <- if (!is.null(opt[["outdir"]])) opt[["outdir"]] else file.path("epi_results", sample)
blacklist<- if (!is.null(opt[["blacklist"]])) opt[["blacklist"]] else "~/.local/hg38.blacklist.bed"
threads  <- if (!is.null(opt[["threads"]])) as.integer(opt[["threads"]]) else 8L
overwrite<- isTRUE(opt[["overwrite"]]) || identical(opt[["overwrite"]], "TRUE") || identical(opt[["overwrite"]], "true")

expand <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
frags  <- expand(frags_in)
bl     <- expand(blacklist)
outdir <- expand(outdir)

# --- Minimal Snakemake alignment ---
# Snakemake expects outputs under: <outdir>/epiAneufinder_results/
outdir <- file.path(outdir, "epiAneufinder_results")
outdir <- expand(outdir)

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("=== epiAneufinder run ===")
message("sample:    ", sample)
message("fragments: ", frags)
message("outdir:    ", outdir)
message("blacklist: ", bl)
message("threads:   ", threads)

# If not overwriting and results exist, honor reuse.existing downstream
reuse <- TRUE
if (overwrite) reuse <- FALSE

# --------------- Run epiAneufinder (no post-processing) ---------------
epiAneufinder(
  input          = frags,
  outdir         = outdir,
  blacklist      = bl,
  windowSize     = 1e5,
  genome         = "BSgenome.Hsapiens.UCSC.hg38",
  exclude        = c("chrX","chrY","chrM"),
  reuse.existing = reuse,
  title_karyo    = paste0("Karyogram - ", sample),
  ncores         = threads,
  minFrags       = 1000,
  minsizeCNV     = 0,
  k              = 4,
  plotKaryo      = TRUE
)

message("epiAneufinder completed.")

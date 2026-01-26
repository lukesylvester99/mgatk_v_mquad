#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(epiAneufinder)
})

# --------------- CLI (same style as your epi runner) ---------------
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

sample  <- if (!is.null(opt[["sample"]])) opt[["sample"]] else "sample"
path_in <- must("results_table")
out_dir <- must("out_dir")

expand <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
path_in <- expand(path_in)
out_dir <- expand(out_dir)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("=== epiAneufinder clone calling ===")
message("sample:       ", sample)
message("results_table:", path_in)
message("out_dir:      ", out_dir)

# --------
# Read results_table robustly (your original logic)
# --------

message("reading results_table...")

first_line <- readLines(path_in, n = 1)
use_tab <- grepl("\t", first_line)

res <- read.table(
  path_in,
  header = TRUE,
  sep = if (use_tab) "\t" else "",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

colnames(res) <- gsub("\\.", "-", colnames(res))

message("results_table is loaded.")

# --------
# Clone calling (unchanged)
# --------

message("splitting subclones...")
subclones <- split_subclones(
  res,
  tree_depth = 4,
  plot_tree  = TRUE,
  plot_path  = file.path(out_dir, "subclones.pdf"),
  plot_width = 4,
  plot_height= 3
)
message("Subclones assigned.")

message("making subclones tsv file...")
# Save clone assignments (portable)
subclone_tsv <- file.path(out_dir, paste0("subclones_", sample, ".tsv"))
write.table(
  subclones,
  file = subclone_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message("subclones tsv written.")

# --------
# Annotated karyotype plot (unchanged)
# --------

message("making annotated karyotype...")
annot_dt <- data.frame(
  cell  = subclones$cell,
  annot = paste0("Clone", subclones$subclone)
)

plot_karyo_annotated(
  res_table = res,
  plot_path = file.path(out_dir, "karyo_annotated.png"),
  annot_dt  = annot_dt
)

message("Annotated karyotype written.")

# --------
# Cancerous barcode list for downstream (epi clone IDs 2, 3, 4)
# --------

message("writing cancerous cell barcodes for mquad...")
cnv_keep <- c(1, 2, 3, 4)
cnv_cells <- subclones$cell[subclones$subclone %in% cnv_keep]

cnv_tsv <- file.path(out_dir, paste0(sample, "_cells_for_mquad.tsv"))
writeLines(c("cell_barcode", cnv_cells), con = cnv_tsv)

message("Done. Outputs written to: ", out_dir)

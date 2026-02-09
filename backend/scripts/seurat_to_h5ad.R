args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript seurat_to_h5ad.R <input.(rds|h5seurat)> <output.h5ad>")
}

input_path <- args[1]
output_path <- args[2]

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

if (grepl("\\.rds$", input_path, ignore.case = TRUE)) {
  obj <- readRDS(input_path)
  tmp_h5seurat <- tempfile(fileext = ".h5seurat")
  SaveH5Seurat(obj, filename = tmp_h5seurat, overwrite = TRUE)
  Convert(tmp_h5seurat, dest = "h5ad", overwrite = TRUE)
  converted <- paste0(tmp_h5seurat, ".h5ad")
  file.copy(converted, output_path, overwrite = TRUE)
} else if (grepl("\\.h5seurat$", input_path, ignore.case = TRUE)) {
  Convert(input_path, dest = "h5ad", overwrite = TRUE)
  converted <- paste0(input_path, ".h5ad")
  file.copy(converted, output_path, overwrite = TRUE)
} else {
  stop("Unsupported file extension. Use .rds or .h5seurat")
}

cat("Converted to:", output_path, "\n")


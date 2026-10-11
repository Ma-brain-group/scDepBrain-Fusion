# Build the ~8,000-cell demo Seurat object shipped with the repository.
# Source: 20k-cell stratified subset of the integrated MDD brain atlas.
# Keeps raw counts only (the pipeline normalises with run_normalize: true),
# genes detected in >= 1 cell, and cell-type annotations + basic QC metrics.
# No donor / sample / batch identifiers are kept.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
set.seed(2026)
SRC <- "path/to/atlas_subset_20k.rds"   # stratified 20k-cell subset of the full atlas
OUT <- "example/data/example_brain_8k.rds"
PER_TYPE <- 400   # cells per subtype (final_anno): ~75% neurons, close to the full atlas (68%)

obj <- readRDS(SRC)
keep <- unlist(lapply(split(colnames(obj), obj$final_anno), function(x)
  if (length(x) > PER_TYPE) sample(x, PER_TYPE) else x))
obj <- obj[, keep]
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- counts[Matrix::rowSums(counts) > 0, ]
stopifnot(all(counts@x %% 1 == 0))

# Anonymise barcodes (cell_0001 ...)
md <- obj@meta.data[, c("anno", "final_anno", "nCount_RNA", "nFeature_RNA"), drop = FALSE]
new_ids <- sprintf("cell_%04d", seq_len(ncol(counts)))
colnames(counts) <- new_ids
rownames(md) <- new_ids
demo <- CreateSeuratObject(counts = counts, meta.data = md, assay = "RNA")
demo$nCount_RNA <- md$nCount_RNA; demo$nFeature_RNA <- md$nFeature_RNA

cat("Cells:", ncol(demo), " Genes:", nrow(demo), "\n")
print(table(demo$anno)); print(table(demo$final_anno))
saveRDS(demo, OUT, compress = "xz")
cat("Saved", OUT, round(file.size(OUT) / 1e6, 1), "MB\n")

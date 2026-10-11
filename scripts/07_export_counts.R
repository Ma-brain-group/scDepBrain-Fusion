# =====================================================================
# scDRS step C1: export RAW counts + metadata from the Seurat object.
#
# Counts are written as the CSC arrays of the genes x cells dgCMatrix
# (x, i, p), which are exactly the CSR arrays of the cells x genes
# matrix that AnnData needs. This avoids SeuratDisk::Convert, which puts
# the log-normalised `data` slot into X (wrong for --flag-raw-count True).
# =====================================================================
source(snakemake@params$common)
setup_logging()
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(hdf5r); library(data.table) })
apply_compat_patches()

p <- snakemake@params
keep_cols <- unlist(p$keep_cols)

ts("Loading ", snakemake@input$rds)
obj <- readRDS(snakemake@input$rds)
ts("Reductions available: ", paste(Reductions(obj), collapse = ", "))

# Batch-corrected embedding for the scDRS neighbour graph (heterogeneity test)
emb <- NULL
if (!is.null(p$reduction) && length(p$reduction) && nzchar(p$reduction)) {
  if (!p$reduction %in% Reductions(obj))
    stop("scdrs.neighbors_reduction '", p$reduction, "' not in Seurat object; available: ",
         paste(Reductions(obj), collapse = ", "))
  emb <- Embeddings(obj, p$reduction)
  ts("Using reduction '", p$reduction, "': ", ncol(emb), " dims")
}
obj <- prep_seurat_rna(obj, p$assay)
obj <- check_meta(obj, p$celltype_col, p$celltype_col)

counts <- GetAssayData(obj, slot = "counts", assay = p$assay)
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
if (length(counts@x) == 0) stop("counts layer is empty — scDRS needs raw counts")
frac_int <- mean(counts@x[seq_len(min(1e6, length(counts@x)))] %% 1 == 0)
ts("Counts: ", nrow(counts), " genes x ", ncol(counts), " cells, nnz = ",
   length(counts@x), ", fraction integer = ", round(frac_int, 4))
if (frac_int < 0.99) stop("counts layer does not look like raw integer counts")

ts("Writing ", snakemake@output$counts)
f <- H5File$new(snakemake@output$counts, mode = "w")
f[["x"]] <- as.numeric(counts@x)
f[["i"]] <- counts@i
f[["p"]] <- counts@p
f[["genes"]] <- rownames(counts)
f[["cells"]] <- colnames(counts)
f$close_all()

md <- obj@meta.data
cols <- intersect(keep_cols, colnames(md))
ts("Metadata columns kept: ", paste(cols, collapse = ", "))
if (length(setdiff(keep_cols, cols)))
  ts("Not found (skipped): ", paste(setdiff(keep_cols, cols), collapse = ", "))
fwrite(data.table(cell_id = rownames(md), md[, cols, drop = FALSE]), snakemake@output$meta)

# Empty file = no embedding given -> build_h5ad computes PCA itself
if (is.null(emb)) {
  file.create(snakemake@output$embedding)
} else {
  emb <- emb[colnames(obj), , drop = FALSE]
  fwrite(data.table(cell_id = rownames(emb), emb), snakemake@output$embedding)
}
ts("DONE")

# =====================================================================
# scPagwas step B1 (GWAS-independent, runs once):
#   Seurat -> Single_data_input -> Pathway_pcascore_run -> pa_meanexpr
#   then split data_mat / pca_scCell_mat into N cell chunks.
# =====================================================================
source(snakemake@params$common)
setup_logging()
suppressPackageStartupMessages({ library(Seurat); library(scPagwas) })
apply_compat_patches()

p <- snakemake@params
group_cols <- unlist(p$group_cols)
pathway_list <- get_resource(p$pathway_list)

ts("Loading ", snakemake@input$rds)
Single_data <- readRDS(snakemake@input$rds)
Single_data <- prep_seurat_rna(Single_data, p$assay)
Single_data <- check_meta(Single_data, p$celltype_col, group_cols)
if (isTRUE(p$run_normalize)) {
  ts("NormalizeData (LogNormalize)")
  Single_data <- NormalizeData(Single_data, assay = p$assay)
}

# Old scPagwas columns would clash with the new ones
old <- grep("^scPagwas\\.|^Random_Correct_BG", colnames(Single_data@meta.data), value = TRUE)
for (col in old) {
  Single_data@meta.data[[paste0(col, "_prev")]] <- Single_data@meta.data[[col]]
  Single_data@meta.data[[col]] <- NULL
}

Idents(Single_data) <- factor(Single_data@meta.data[[p$celltype_col]])
ts("Cells: ", ncol(Single_data), "  Genes: ", nrow(Single_data))
print(table(Idents(Single_data)))

ts("Single_data_input")
Pagwas <- scPagwas::Single_data_input(Pagwas = list(), assay = p$assay,
                                      Single_data = Single_data,
                                      Pathway_list = pathway_list)
Single_data <- Single_data[, colnames(Pagwas$data_mat)]

ts("Pathway_pcascore_run")
Pagwas <- scPagwas::Pathway_pcascore_run(Pagwas = Pagwas, Pathway_list = pathway_list)

ts("pa_meanexpr")
Pagwas <- scPagwas::pa_meanexpr(Pagwas)

ts("data_mat: ", paste(dim(Pagwas$data_mat), collapse = " x "),
   "  pca_scCell_mat: ", paste(dim(Pagwas$pca_scCell_mat), collapse = " x "))
save(Pagwas, Single_data, file = snakemake@output$preprocessed, compress = TRUE)
ts("Saved ", snakemake@output$preprocessed)

# ---- split into chunks along the cell dimension ----
rm(Single_data); gc()
chunk_files <- unlist(snakemake@output$chunks)
n_chunks <- length(chunk_files)
n_cells <- ncol(Pagwas$data_mat)
chunk_idx <- rep(seq_len(n_chunks), each = ceiling(n_cells / n_chunks), length.out = n_cells)
pca_scCell_mat <- Pagwas$pca_scCell_mat
data_mat <- Pagwas$data_mat
for (i in seq_len(n_chunks)) {
  sel <- which(chunk_idx == i)
  chunk <- Pagwas
  chunk$pca_scCell_mat <- pca_scCell_mat[, sel, drop = FALSE]
  chunk$data_mat <- data_mat[, sel, drop = FALSE]
  # chunk files are named chunk_<i>_pagwas.RData; order matches expand()
  save(chunk, file = chunk_files[i], compress = TRUE)
  ts("Chunk ", i, ": n = ", length(sel), " -> ", chunk_files[i])
  rm(chunk); gc()
}
ts("DONE")

# =====================================================================
# scPagwas step B4 (per trait): merge chunks -> gPAS -> PCC top genes
#   -> TRS / upTRS / downTRS -> background-corrected per-cell P
#   -> Merge_celltype_p for every annotation column.
# Uses the validated speed patches (see common.R).
# =====================================================================
source(snakemake@params$common)
setup_logging()
suppressPackageStartupMessages({
  library(Seurat); library(scPagwas); library(Matrix); library(data.table)
})
apply_compat_patches()
patch_scpagwas_fast()

p <- snakemake@params
group_cols <- unlist(p$group_cols)
ntop <- as.integer(p$n_topgenes)

# ---- merge per-chunk regression results (order = chunk 1..N) ----
ts("pcc_method = ", p$pcc_method, "; n_topgenes = ", ntop, "; bg_n_topgenes = ", p$bg_n_topgenes)
ts("Merging ", length(snakemake@input$chunks), " chunks")
pmat_list <- lapply(unlist(snakemake@input$chunks), function(f) { load(f); pmat })
pmat_merge <- do.call(rbind, pmat_list)
rm(pmat_list); gc()
ts("pmat_merge: ", paste(dim(pmat_merge), collapse = " x "))

ts("Loading ", snakemake@input$preprocessed)
load(snakemake@input$preprocessed)   # -> Pagwas, Single_data

# ---- gPAS (same as scPagwas_perform_score(remove_outlier = TRUE)) ----
ts("Merge_gPas")
gpas <- scPagwas::Merge_gPas(Pagwas, pmat_merge)

# ---- gene-gPAS correlation and TRS gene sets ----
if (p$pcc_method == "scGet_PCC") {
  # Lab's original scPagwas workflow: PCC over all cells; top / bottom
  # n_topgenes of the same PCC give TRS and downTRS (no upTRS).
  ts("scGet_PCC (all cells)")
  PCC <- scPagwas::scGet_PCC(scPagwas.gPAS.score = gpas[colnames(Pagwas$data_mat)],
                             data_mat = Pagwas$data_mat)
  pcc <- setNames(PCC[, 1], rownames(PCC))
  topgenes  <- names(sort(pcc, decreasing = TRUE))[seq_len(ntop)]
  downgenes <- names(sort(pcc, decreasing = FALSE))[seq_len(ntop)]
  gene_sets <- list(topgenes, downgenes)
  set_names <- c("scPagwas.TRS.Score", "scPagwas.downTRS.Score")
  fwrite(data.table(rank = seq_len(ntop), top = topgenes, down = downgenes,
                    PCC_top = pcc[topgenes], PCC_down = pcc[downgenes]),
         snakemake@output$topgenes)
} else if (p$pcc_method == "Corr_Random_topbottom") {
  # Corr_Random PCC (random cell subsets); top / bottom n_topgenes of the same
  # PCC give TRS and downTRS (no upTRS).
  ts("Corr_Random (all cells; top / bottom of one PCC)")
  PCC <- scPagwas::Corr_Random(Pagwas$data_mat, gpas, seed = p$seed, random = TRUE,
                               Nrandom = p$corr_nrandom, Nselect = p$corr_nselect)
  topgenes  <- names(sort(PCC, decreasing = TRUE))[seq_len(ntop)]
  downgenes <- names(sort(PCC, decreasing = FALSE))[seq_len(ntop)]
  gene_sets <- list(topgenes, downgenes)
  set_names <- c("scPagwas.TRS.Score", "scPagwas.downTRS.Score")
  fwrite(data.table(rank = seq_len(ntop), top = topgenes, down = downgenes,
                    PCC_top = PCC[topgenes], PCC_down = PCC[downgenes]),
         snakemake@output$topgenes)
} else if (p$pcc_method == "Corr_Random") {
  # Memory-light variant (random cell subsets) with separate up / down PCC
  ts("Corr_Random (all / up / down)")
  PCC <- scPagwas::Corr_Random(Pagwas$data_mat, gpas, seed = p$seed, random = TRUE,
                               Nrandom = p$corr_nrandom, Nselect = p$corr_nselect)
  a1 <- which(gpas >= mean(gpas))
  a2 <- which(gpas < mean(gpas))
  PCC_up   <- scPagwas::Corr_Random(scPagwas.gPAS.score = gpas[a1], data_mat = Pagwas$data_mat[, a1])
  PCC_down <- scPagwas::Corr_Random(scPagwas.gPAS.score = gpas[a2], data_mat = Pagwas$data_mat[, a2])
  topgenes  <- names(sort(PCC, decreasing = TRUE))[seq_len(ntop)]
  upgenes   <- names(sort(PCC_up, decreasing = TRUE))[seq_len(ntop)]
  downgenes <- names(sort(PCC_down, decreasing = FALSE))[seq_len(ntop)]
  gene_sets <- list(topgenes, upgenes, downgenes)
  set_names <- c("scPagwas.TRS.Score", "scPagwas.upTRS.Score", "scPagwas.downTRS.Score")
  fwrite(data.table(rank = seq_len(ntop), top = topgenes, up = upgenes, down = downgenes,
                    PCC_top = PCC[topgenes]),
         snakemake@output$topgenes)
} else {
  stop("scpagwas.pcc_method must be 'scGet_PCC', 'Corr_Random_topbottom' or 'Corr_Random'")
}

# ---- TRS ----
ts("AddModuleScore: ", paste(set_names, collapse = ", "))
Single_data <- Seurat::AddModuleScore(Single_data, assay = p$assay,
                                      features = gene_sets, name = set_names)
for (k in seq_along(set_names)) {
  tmp <- paste0(set_names[k], k)
  Single_data@meta.data[[set_names[k]]] <- Single_data@meta.data[[tmp]]
  Single_data@meta.data[[tmp]] <- NULL
}
Single_data$scPagwas.gPAS.score <- gpas[colnames(Single_data)]

# ---- background-corrected per-cell P ----
ts("Get_CorrectBg_p (iters = ", p$iters_singlecell, ", n_topgenes = ", p$bg_n_topgenes, ")")
correct_pdf <- scPagwas::Get_CorrectBg_p(
  Single_data = Single_data,
  scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score,
  iters_singlecell = p$iters_singlecell,
  n_topgenes = p$bg_n_topgenes,
  scPagwas_topgenes = topgenes)
# Attribute P by cell name (fixes the named-vector bug of the original)
pooled_p <- setNames(correct_pdf$pooled_p, rownames(correct_pdf))
Single_data$Random_Correct_BG_p    <- pooled_p[colnames(Single_data)]
Single_data$Random_Correct_BG_adjp <- setNames(correct_pdf$adj_p, rownames(correct_pdf))[colnames(Single_data)]
Single_data$Random_Correct_BG_z    <- setNames(correct_pdf$pooled_z, rownames(correct_pdf))[colnames(Single_data)]
ts("Cells with valid pooled_p: ", sum(!is.na(Single_data$Random_Correct_BG_p)),
   " / ", ncol(Single_data))

# ---- per-cell table ----
md <- Single_data@meta.data
score_cols <- intersect(c("scPagwas.gPAS.score", set_names, "Random_Correct_BG_p",
                          "Random_Correct_BG_adjp", "Random_Correct_BG_z"), colnames(md))
per_cell <- data.table(cell_id = colnames(Single_data),
                       md[, c(group_cols, score_cols), drop = FALSE])
fwrite(per_cell, snakemake@output$per_cell)

# ---- cell-type level P (Merge_celltype_p) for each annotation ----
out_ct <- unlist(snakemake@output$celltype)
for (k in seq_along(group_cols)) {
  g <- group_cols[k]
  ts("Merge_celltype_p: ", g)
  ct <- scPagwas::Merge_celltype_p(single_p = pooled_p,
                                   celltype = md[names(pooled_p), g])
  ct <- as.data.frame(ct)
  ct$n_cell <- as.integer(table(md[[g]])[as.character(ct$celltype)])
  ct$fdr <- p.adjust(ct$pvalue, method = "BH")
  print(ct)
  fwrite(ct, out_ct[k])
}

if (isTRUE(p$save_seurat)) {
  ts("Saving Seurat with scPagwas scores -> ", p$seurat_out)
  saveRDS(Single_data, p$seurat_out)
}
ts("DONE")

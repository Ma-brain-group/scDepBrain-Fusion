# =====================================================================
# Shared helpers for the scDepBrain-Fusion R steps.
# Sourced via snakemake@params$common.
# =====================================================================

ts <- function(...) cat("[", format(Sys.time()), "] ", ..., "\n", sep = "")

# Redirect stdout + messages to the Snakemake log file
setup_logging <- function() {
  if (exists("snakemake") && length(snakemake@log) > 0) {
    con <- file(snakemake@log[[1]], open = "wt")
    sink(con)
    sink(con, type = "message")
  }
}

# Single-thread BLAS (bigstatsr check) + silence SeuratObject deprecations
apply_compat_patches <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
  }
  options(default.nproc.blas = 1)
  options(bigstatsr.check.parallel.blas = FALSE)
  suppressMessages(utils::assignInNamespace(
    ".Deprecate", function(...) invisible(NULL), ns = "SeuratObject"))
}

# Resolve a scPagwas built-in data object by name, or read an .rds path
get_resource <- function(x) {
  if (file.exists(x)) return(readRDS(x))
  utils::data(list = x, package = "scPagwas", envir = environment())
  get(x, envir = environment())
}

# Keep only the RNA assay, collapsed to a v3 Assay (scPagwas needs slot API)
prep_seurat_rna <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  if (inherits(obj[[assay]], "Assay5")) {
    ts("Assay5 -> JoinLayers + v3 Assay")
    obj <- JoinLayers(obj, assay = assay)
    obj[[assay]] <- as(obj[[assay]], "Assay")
  }
  for (a in setdiff(Assays(obj), assay)) obj[[a]] <- NULL
  obj@reductions <- list()
  gc()
  obj
}

# Validate that metadata columns exist; drop cells with NA in `celltype_col`
check_meta <- function(obj, celltype_col, group_cols) {
  miss <- setdiff(group_cols, colnames(obj@meta.data))
  if (length(miss)) stop("Metadata columns missing from Seurat object: ",
                         paste(miss, collapse = ", "))
  na <- is.na(obj@meta.data[[celltype_col]])
  if (any(na)) {
    ts("Dropping ", sum(na), " cells with NA in ", celltype_col)
    obj <- obj[, !na]
  }
  obj
}

# ---------------------------------------------------------------------
# Speed patches for scPagwas (numerically equivalent to the originals)
#  1. Get_CorrectBg_p: sparse-aware row variance (10+ h -> seconds)
#  2. get_p_from_empi_null: sort + findInterval instead of O(n*m) sapply
# ---------------------------------------------------------------------
patch_scpagwas_fast <- function() {
  patched_Get_CorrectBg_p <- function(Single_data, scPagwas.TRS.Score,
                                      iters_singlecell, n_topgenes,
                                      scPagwas_topgenes, assay = "RNA") {
    gene_matrix <- GetAssayData(Single_data, slot = "data", assay = assay)
    mat_ctrl_raw_score <- matrix(0, nrow = ncol(gene_matrix),
                                 ncol = iters_singlecell)
    dic_ctrl_list <- list()
    for (i in seq_len(iters_singlecell)) {
      set.seed(i)
      dic_ctrl_list[[i]] <- sample(rownames(Single_data), n_topgenes)
      Single_data <- Seurat::AddModuleScore(Single_data, assay = assay,
                                            list(dic_ctrl_list[[i]]),
                                            name = c("contr_genes"))
      mat_ctrl_raw_score[, i] <- Single_data$contr_genes1
      Single_data$contr_genes1 <- NULL
      if (i %% 10 == 0) ts("  control iteration ", i, "/", iters_singlecell)
    }

    genes <- intersect(rownames(Single_data), rownames(gene_matrix))
    scPagwas_topgenes <- intersect(scPagwas_topgenes, genes)
    gene_matrix <- gene_matrix[genes, ]

    n  <- ncol(gene_matrix)
    mu <- Matrix::rowMeans(gene_matrix)
    ss <- Matrix::rowSums(gene_matrix^2)
    df_gene <- data.frame(gene = genes, var = as.numeric((ss - n * mu^2) / (n - 1)))
    rownames(df_gene) <- df_gene$gene

    v_var_ratio_c2t <- vapply(seq_len(iters_singlecell),
                              function(j) sum(df_gene[dic_ctrl_list[[j]], "var"]),
                              numeric(1))
    v_var_ratio_c2t <- v_var_ratio_c2t / sum(df_gene[scPagwas_topgenes, "var"])

    correct_pdf <- scPagwas:::correct_background(scPagwas.TRS.Score,
                                                 mat_ctrl_raw_score,
                                                 v_var_ratio_c2t)
    rownames(correct_pdf) <- colnames(Single_data)
    correct_pdf
  }

  patched_get_p_from_empi_null <- function(v_t, v_t_null) {
    sorted_null <- sort(v_t_null)
    v_pos <- findInterval(v_t, sorted_null)
    (length(v_t_null) - v_pos + 1) / (length(v_t_null) + 1)
  }

  utils::assignInNamespace("Get_CorrectBg_p", patched_Get_CorrectBg_p, ns = "scPagwas")
  utils::assignInNamespace("get_p_from_empi_null", patched_get_p_from_empi_null, ns = "scPagwas")
  ts("scPagwas speed patches applied (Get_CorrectBg_p, get_p_from_empi_null)")
}

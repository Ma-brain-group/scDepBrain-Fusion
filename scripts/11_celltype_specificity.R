# =====================================================================
# Cell-type expression specificity for MAGMA and LDSC (per annotation level).
#   avg  = AverageExpression (mean normalised expression per cell type)
#   spec = avg[g, ct] / sum_ct avg[g, ct]                    (Bryois et al. 2020)
# Outputs
#   magma/covar_<group>.txt           GENE (Entrez) + one continuous column per cell type
#                                     (input for `magma --gene-covar`)
#   ldsc/genesets/<group>__<ct>.geneset.txt   top `top_pct`% genes by specificity (symbols)
#   ldsc/genesets/all_gene.geneset.txt        all expressed genes (LDSC-SEG control set)
#   specificity/celltype_names.csv    group, celltype, safe name used in file / column names
# =====================================================================
source(snakemake@params$common)
setup_logging()
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
apply_compat_patches()

p <- snakemake@params
group_cols <- unlist(p$group_cols)
safe <- function(x) { s <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x)); ifelse(grepl("^[0-9]", s), paste0("ct_", s), s) }

gl <- read.table(p$gene_loc, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)[, c(1, 6)]
colnames(gl) <- c("entrez", "symbol")
sym2ent <- setNames(as.character(gl$entrez), gl$symbol)

ts("Loading ", snakemake@input$rds)
obj <- readRDS(snakemake@input$rds)
obj <- prep_seurat_rna(obj, p$assay)
obj <- check_meta(obj, p$celltype_col, group_cols)
if (isTRUE(p$run_normalize)) obj <- NormalizeData(obj, assay = p$assay)

dir.create(snakemake@output$genesets, recursive = TRUE, showWarnings = FALSE)
covar_files <- unlist(snakemake@output$covar)
name_map <- list(); all_genes <- character(0)

for (k in seq_along(group_cols)) {
  g <- group_cols[k]
  lv <- sort(unique(as.character(obj@meta.data[[g]])))
  ts("AverageExpression by ", g, " (", length(lv), " groups)")
  avg <- as.matrix(AverageExpression(obj, assays = p$assay, group.by = g, return.seurat = FALSE)[[1]])
  # AverageExpression can alter group names (e.g. "_" -> "-"); map back by alphanumeric key
  key <- function(x) gsub("[^A-Za-z0-9]", "", x)
  m <- match(key(colnames(avg)), key(lv))
  # Seurat 5 prefixes numeric group names with "g" (e.g. clusters "0" -> "g0")
  m[is.na(m)] <- match(key(sub("^g", "", colnames(avg)[is.na(m)])), key(lv))
  colnames(avg) <- lv[m]
  if (anyNA(colnames(avg))) stop("Could not match AverageExpression columns to levels of ", g)
  avg <- avg[rowSums(avg) > 0, , drop = FALSE]
  all_genes <- union(all_genes, rownames(avg))
  sn <- safe(colnames(avg))
  if (anyDuplicated(sn)) stop("Cell-type names of ", g, " collide after sanitising: ", paste(colnames(avg), collapse = ", "))
  name_map[[g]] <- data.frame(group = g, celltype = colnames(avg), safe = sn,
                              set = sprintf("%s__%s", safe(g), sn))

  # LDSC: top genes by specificity among all expressed genes (symbols)
  spec <- sweep(avg, 1, rowSums(avg), "/")
  n_top <- floor(nrow(spec) * p$top_pct / 100)
  for (j in seq_len(ncol(spec))) {
    top <- rownames(spec)[order(spec[, j], decreasing = TRUE)[seq_len(n_top)]]
    writeLines(top, file.path(snakemake@output$genesets, sprintf("%s__%s.geneset.txt", safe(g), sn[j])))
  }
  ts("  LDSC gene sets: top ", p$top_pct, "% = ", n_top, " genes per cell type")

  # MAGMA: genes with an Entrez ID; duplicated Entrez -> mean specificity
  a2 <- avg[rownames(avg) %in% names(sym2ent), , drop = FALSE]
  sp2 <- sweep(a2, 1, rowSums(a2), "/")
  ent <- sym2ent[rownames(sp2)]
  sp2 <- rowsum(sp2, ent) / as.vector(table(ent)[unique(sort(ent))])
  colnames(sp2) <- sn
  write.table(data.frame(GENE = rownames(sp2), sp2, check.names = FALSE), covar_files[k],
              quote = FALSE, sep = " ", row.names = FALSE)
  ts("  MAGMA covariates: ", nrow(a2), " genes mapped -> ", nrow(sp2), " Entrez IDs")
}
writeLines(all_genes, file.path(snakemake@output$genesets, "all_gene.geneset.txt"))
write.csv(do.call(rbind, name_map), snakemake@output$names, row.names = FALSE)
ts("all_gene control set: ", length(all_genes), " genes. DONE")

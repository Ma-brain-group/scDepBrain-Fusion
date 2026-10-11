# =====================================================================
# scPagwas step B5 (per trait): cell-type-level analysis on the
# celltype_col used as Idents in sc_prep.
#   Link_pathway_blocks_gwas(celltype = TRUE) -> Pagwas_perform_regression
#   -> Boot_evaluate (bootstrap P per cell type)
# =====================================================================
source(snakemake@params$common)
setup_logging()
suppressPackageStartupMessages({ library(Seurat); library(scPagwas); library(Matrix) })
apply_compat_patches()

p <- snakemake@params
block_annotation <- get_resource(p$block_annotation)
chrom_ld <- get_resource(p$chrom_ld)
dir.create(p$backingpath, showWarnings = FALSE, recursive = TRUE)

ts("Loading preprocessed single-cell Pagwas")
load(snakemake@input$preprocessed)   # -> Pagwas, Single_data
rm(Single_data); gc()
Pagwas_sc <- Pagwas
load(snakemake@input$gwas)           # -> Pagwas (GWAS side)
Pagwas <- c(Pagwas, Pagwas_sc)
rm(Pagwas_sc); gc()

ts("Pathway_annotation_input")
Pagwas <- scPagwas::Pathway_annotation_input(Pagwas = Pagwas,
                                             block_annotation = block_annotation)

ts("Link_pathway_blocks_gwas (singlecell = FALSE, celltype = TRUE)")
Pagwas <- scPagwas::Link_pathway_blocks_gwas(Pagwas = Pagwas, chrom_ld = chrom_ld,
                                             singlecell = FALSE, celltype = TRUE,
                                             backingpath = p$backingpath, n.cores = 1)

ts("Pagwas_perform_regression")
Pagwas$lm_results <- scPagwas::Pagwas_perform_regression(
  Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
print(Pagwas$lm_results)

ts("Boot_evaluate (iters = ", p$bootstrap_iters, ", part = ", p$bootstrap_part, ")")
Pagwas <- scPagwas::Boot_evaluate(Pagwas, bootstrap_iters = p$bootstrap_iters,
                                  part = p$bootstrap_part)
boot <- as.data.frame(Pagwas$bootstrap_results)
print(boot)
write.csv(boot, snakemake@output[[1]], row.names = FALSE)

unlink(p$backingpath, recursive = TRUE)
ts("DONE")

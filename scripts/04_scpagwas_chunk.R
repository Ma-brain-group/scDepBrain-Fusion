# =====================================================================
# scPagwas step B3 (per trait x chunk): Link_pathway_blocks_gwas in
# single-cell mode -> Pathway_sclm_results for this chunk of cells.
# Memory dominated by LD matrices (~200-400 GB for ~60k cells).
# =====================================================================
source(snakemake@params$common)
setup_logging()
suppressPackageStartupMessages({ library(scPagwas) })
apply_compat_patches()

p <- snakemake@params
block_annotation <- get_resource(p$block_annotation)
chrom_ld <- get_resource(p$chrom_ld)
dir.create(p$backingpath, showWarnings = FALSE, recursive = TRUE)

ts("Loading GWAS Pagwas")
load(snakemake@input$gwas)          # -> Pagwas
Pagwas_gwas <- Pagwas; rm(Pagwas)

ts("Loading ", snakemake@input$chunk)
load(snakemake@input$chunk)         # -> chunk
Pagwas <- c(Pagwas_gwas, chunk)
rm(Pagwas_gwas, chunk); gc()
ts("data_mat: ", paste(dim(Pagwas$data_mat), collapse = " x "))

ts("Pathway_annotation_input")
Pagwas <- scPagwas::Pathway_annotation_input(Pagwas = Pagwas,
                                             block_annotation = block_annotation)

ts("Link_pathway_blocks_gwas (singlecell = TRUE, celltype = FALSE)")
Pagwas <- scPagwas::Link_pathway_blocks_gwas(Pagwas = Pagwas, chrom_ld = chrom_ld,
                                             singlecell = TRUE, celltype = FALSE,
                                             backingpath = p$backingpath, n.cores = 1)

pmat <- Pagwas$Pathway_sclm_results
ts("Pathway_sclm_results: ", paste(dim(pmat), collapse = " x "))
save(pmat, file = snakemake@output[[1]], compress = TRUE)
unlink(p$backingpath, recursive = TRUE)
ts("DONE")

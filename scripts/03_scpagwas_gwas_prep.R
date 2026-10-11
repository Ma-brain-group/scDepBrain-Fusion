# =====================================================================
# scPagwas step B2 (per trait): GWAS_summary_input + SnpToGene
# =====================================================================
source(snakemake@params$common)
setup_logging()
suppressPackageStartupMessages({ library(scPagwas); library(bigreadr) })
apply_compat_patches()

p <- snakemake@params
block_annotation <- get_resource(p$block_annotation)

ts("Reading ", snakemake@input[[1]])
gwas_data <- bigreadr::fread2(snakemake@input[[1]])
ts("Rows: ", nrow(gwas_data), "  Cols: ", paste(colnames(gwas_data), collapse = " "))

Pagwas <- scPagwas::GWAS_summary_input(Pagwas = list(), gwas_data = gwas_data,
                                       maf_filter = p$maf_filter)
ts("After maf_filter = ", p$maf_filter, ": ", nrow(Pagwas$gwas_data), " SNPs")

ts("SnpToGene (marg = ", p$marg, ")")
Pagwas$snp_gene_df <- scPagwas::SnpToGene(gwas_data = Pagwas$gwas_data,
                                          block_annotation = block_annotation,
                                          marg = p$marg)
Pagwas$rawPathway_list <- get_resource(p$pathway_list)

save(Pagwas, file = snakemake@output[[1]], compress = TRUE)
ts("DONE")

# Output files

Everything is written under `outdir` from the config. `<trait>` is the trait key, and `<group>` is each entry of `group_cols`.

## `integrated/`: start here

### `integrated/four_method_<group>.csv`

One row per trait × cell type.

| Column | Meaning |
|---|---|
| `trait`, `cell_type` | identifiers |
| `MAGMA_P` | MAGMA gene-property P (continuous specificity, `--gene-covar`) |
| `LDSC_P` | S-LDSC one-sided P of the cell-type coefficient (NA if LDSC is off for the trait) |
| `scPagwas_P` | scPagwas `Merge_celltype_p` |
| `scDRS_P` | scDRS group association MC P |
| `n_sig`, `n_avail` | methods with P < `integration.alpha` / methods with a P value |
| `significant_methods` | combination label used for the colour, e.g. `MAGMA+scDRS`, `All four`, `None` |
| `mean_neglog10P` | mean of −log10 P over the available methods (bar length) |

### `integrated/four_method_<group>.png|pdf`

Bar chart with one panel per trait. Cell types are sorted by their mean over traits; colour shows the significant-method combination and `(k/n)` the number of significant / available methods. `four_method_all.csv` stacks all annotation levels.

## `magma/`, `ldsc/`, `specificity/`

| File | Content |
|---|---|
| `magma/<trait>_<group>.gsa.out` | MAGMA output: `VARIABLE` (cell type, sanitised name), `BETA`, `BETA_STD`, `SE`, `P` |
| `ldsc/<trait>/ldsc_celltype_results.csv` | `group`, `celltype`, `Coefficient`, `Coefficient_se`, `Coefficient_z`, `LDSC_P` (one-sided), `Enrichment`, `Enrichment_p`, `Prop_SNPs` |
| `ldsc/<trait>/results/<set>.results` | raw S-LDSC output per cell type |
| `ldsc/sumstats/<trait>.sumstats.gz` | munged GWAS (HapMap3) |
| `ldsc/ldscores/<set>/` | per-gene-set annotations and LD scores (reusable for other GWAS on the same single-cell data) |
| `specificity/magma_covar_<group>.txt` | specificity matrix (Entrez × cell type) |
| `specificity/ldsc_genesets/` | top-specific gene sets and the `all_gene` control |
| `specificity/celltype_names.csv` | original cell-type names ↔ sanitised names used in files |

## `summary/`: scPagwas × scDRS details

### `summary/<trait>/summary_<group>.csv`

One row per cell type, sorted by scDRS association P.

| Column | Meaning |
|---|---|
| `trait`, `group_col`, `celltype` | identifiers |
| `scPagwas_p`, `scPagwas_fdr` | scPagwas `Merge_celltype_p`: per-cell background-corrected P combined per cell type; BH FDR across cell types |
| `n_cell` | cells in the group (scPagwas) |
| `scPagwas_boot_t`, `scPagwas_boot_p`, `scPagwas_boot_fdr` | scPagwas cell-type regression bootstrap: t value, one-sided P for a positive effect, FDR. Only for `celltype_col`. |
| `scDRS_n_cell` | cells in the group after scDRS filtering |
| `scDRS_assoc_p`, `scDRS_assoc_z`, `scDRS_assoc_fdr` | scDRS cell-type association (MC test of the top-5% score quantile) |
| `scDRS_hetero_p`, `scDRS_hetero_fdr` | scDRS heterogeneity within the cell type (Geary's C on the KNN graph) |
| `scDRS_n_cell_fdr0.1` | cells with per-cell scDRS FDR < 0.1 |

Note on resolution: scDRS MC P values have a floor of 1 / (n_ctrl + 1), e.g. ≈ 0.001 with 1,000 controls.

### `summary/<trait>/summary_<group>.pdf|png`

Heatmap of −log10 P with methods as rows and cell types as columns; `*` marks BH FDR < 0.05.

### `summary/<trait>/scPagwas_vs_scDRS_percell_corr.csv`

Spearman correlation between the scPagwas TRS and the scDRS normalised score across cells, overall (`ALL`) and within each `celltype_col` group. Positive correlations mean the two methods rank the same cells as trait-relevant.

### `summary/all_traits_celltype_summary.csv`

All `summary_<group>.csv` tables stacked.

## `scpagwas/<trait>/`

| File | Content |
|---|---|
| `celltype_p_<group>.csv` | `celltype`, `pvalue` (`Merge_celltype_p`), `n_cell`, `fdr` |
| `celltype_bootstrap_<celltype_col>.csv` | `Boot_evaluate` output: `bootstrap_estimate`, `bootstrap_error`, `bt_value`, `bp_value`, `bias_corrected_estimate`, `CI_lo`, `CI_hi` (row `Intercept` included) |
| `per_cell_scores.csv.gz` | per cell: annotations, `scPagwas.gPAS.score`, `scPagwas.TRS.Score`, `scPagwas.downTRS.Score` (or `upTRS` too with `Corr_Random`), `Random_Correct_BG_p`, `Random_Correct_BG_adjp`, `Random_Correct_BG_z` |
| `scPagwas_topgenes.csv` | genes used for TRS (top) and downTRS (bottom) with their PCC |
| `gwas_pagwas.RData`, `chunks/` | intermediate objects, kept so later steps can resume |
| `scPagwas_seurat.rds` | only if `save_seurat: true` |

## `scpagwas/shared/`

`pagwas_preprocessed.RData` and `chunk_<i>_pagwas.RData` are GWAS-independent intermediates reused by all traits. They can be large and can be deleted after the run.

## `scdrs/`

| File | Content |
|---|---|
| `data/sc_counts.h5ad` | raw-count AnnData used by scDRS |
| `data/cell_meta.csv.gz`, `data/cov.tsv` | metadata; covariates (`const` + `cov_cols`) |
| `gs/<trait>.gene_pval.tsv`, `gs/<trait>.gs` | MAGMA P values per symbol; scDRS gene set (top `n_max`, z-score weights) |
| `scores/<trait>.score.gz` | per cell: `raw_score`, `norm_score`, `mc_pval`, `pval`, `nlog10_pval`, `zscore` |
| `scores/<trait>.full_score.gz` | as above, plus normalised control scores (used by downstream) |
| `downstream/<trait>.scdrs_group.<group>` | `group`, `n_cell`, `n_ctrl`, `assoc_mcp`, `assoc_mcz`, `hetero_mcp`, `hetero_mcz`, `n_fdr_0.05/0.1/0.2` |

## `gwas/`, `logs/`

- `gwas/<trait>.scpagwas_input.txt`: the formatted GWAS (`chrom pos rsid se beta maf`).
- `logs/<step>/…log`: one log per job. Read these first when something fails. SLURM stdout and stderr go to `logs/slurm/` in the repository folder.

# Configuration reference

All settings live in one YAML file passed with `--configfile` (via `run_pipeline.sh`). Start from [`config/config_template.yaml`](../config/config_template.yaml). Relative paths resolve against the repository folder.

## Top level

| Key | Default | Description |
|---|---|---|
| `methods.magma` | `true` | Run the MAGMA gene-property analysis (needs `magma.binary`, `magma.gene_loc`, `traits.*.magma_genes_raw`). |
| `methods.ldsc` | `true` | Run S-LDSC / LDSC-SEG (needs the `ldsc:` block). scPagwas and scDRS always run. |
| `outdir` | — | Output folder. Created if missing. |
| `env_activate` | `""` | Shell snippet run before every job, e.g. `source ~/miniconda3/etc/profile.d/conda.sh && conda activate scdepbrain-fusion`. Leave empty when jobs inherit the active environment (local runs; SLURM with default `--export=ALL`). |

## `single_cell`

| Key | Default | Description |
|---|---|---|
| `rds` | — | Seurat object (`saveRDS`). |
| `assay` | `RNA` | Assay holding `counts` and `data`. |
| `run_normalize` | `false` | Run `NormalizeData()` before scPagwas. Needed if `data` is not log-normalised. |
| `celltype_col` | — | Main annotation: scPagwas `Idents`, cell-type regression, bootstrap. |
| `group_cols` | — | All annotations tested by `Merge_celltype_p` and scDRS group analysis. Must include `celltype_col`. |
| `keep_meta` | `[]` | Extra metadata copied into the scDRS `.h5ad`. Columns that don't exist are skipped. |

## `traits.<name>`

| Key | Default | Description |
|---|---|---|
| `sumstats` | — | GWAS file. |
| `magma_genes` | — | MAGMA `.genes.out` for the same GWAS (used by scDRS). |
| `magma_genes_raw` | — | MAGMA `.genes.raw` from the same run (used by the MAGMA gene-property analysis). |
| `ldsc` | `true` | Run S-LDSC for this trait. Set to `false` for non-European or trans-ancestry GWAS. |
| `n` | `null` | Constant GWAS sample size, used by LDSC when the file has no N column (`gwas_columns.n`). |
| `format` | `table` | `table` or `vcf` (IEU GWAS-VCF). |
| `columns` | `{}` | Overrides of `gwas_columns` for this trait. |

| `scpagwas` | — | Optional per-trait scPagwas settings that override the global `scpagwas:` block. Allowed keys: `maf_filter`, `marg`, `pcc_method`, `n_topgenes`, `bg_n_topgenes`, `iters_singlecell`, `corr_nrandom`, `corr_nselect`, `seed`. Example: `scpagwas: {marg: 20000, pcc_method: Corr_Random}`. |

Trait names may contain only letters, digits and underscores.

## `gwas_columns`, `gwas_filter`

See [inputs.md](inputs.md#2-gwas-summary-statistics).

| Key | Default | Description |
|---|---|---|
| `gwas_filter.autosomes_only` | `true` | Keep chromosomes 1–22. |
| `gwas_filter.min_maf` | `0.01` | Pre-filter during formatting. |
| `gwas_filter.exclude_region` | chr6 25–35 Mb | MHC; `null` to disable. Match it to the build. |

## `scpagwas`

| Key | Default | Description |
|---|---|---|
| `pathway_list` | `Genes_by_pathway_kegg` | scPagwas data object, or `.rds` of a named list of gene vectors. Others shipped with scPagwas: `genes.by.gobp.pathway`, `genes.by.reactome.pathway`, `genes.by.regulatory.pathway`, `genes.by.tft.pathway`. |
| `block_annotation` | `block_annotation_hg37` | Gene coordinates: `block_annotation_hg37` (GRCh37) or `block_annotation` (GRCh38). |
| `chrom_ld` | `chrom_ld` | LD reference shipped with scPagwas. |
| `maf_filter` | `0.1` | `GWAS_summary_input(maf_filter=)`. |
| `marg` | `10000` | `SnpToGene(marg=)`: bp around a gene that assign a SNP to it. |
| `n_chunks` | `5` | Number of cell chunks for `Link_pathway_blocks_gwas`. Results don't depend on it; it only sets memory per job. |
| `pcc_method` | `scGet_PCC` | `scGet_PCC`: PCC between gPAS and each gene over all cells; top / bottom `n_topgenes` → TRS / downTRS. `Corr_Random_topbottom`: same, but PCC from random cell subsets (`corr_nrandom` × `corr_nselect`), which is lighter on memory. `Corr_Random`: random-subset PCC with separate up / down sets (gPAS ≥ / < mean) → TRS / upTRS / downTRS. |
| `n_topgenes` | `500` | Genes per TRS gene set. |
| `bg_n_topgenes` | `1000` | Control-gene set size in `Get_CorrectBg_p`. |
| `iters_singlecell` | `100` | Control sets in `Get_CorrectBg_p`. |
| `seed` | `1234` | Seed for `Corr_Random`. |
| `run_celltype_bootstrap` | `true` | Run the cell-type regression + `Boot_evaluate`. |
| `bootstrap_iters`, `bootstrap_part` | `200`, `0.5` | `Boot_evaluate` settings. |
| `save_seurat` | `false` | Also save `scPagwas_seurat.rds` with all scores in `meta.data`. Large. |

## `magma`, `scdrs`

| Key | Default | Description |
|---|---|---|
| `magma.binary` | `magma` | MAGMA executable (name on `PATH` or full path). |
| `magma.gene_id` | `auto` | `auto`, `entrez` or `symbol`. |
| `magma.gene_loc` | — | MAGMA gene.loc file (build of the GWAS). Used to map Entrez ↔ symbol for scDRS and the MAGMA covariates, and by default for the LDSC gene coordinates. |
| `scdrs.weight` | `zscore` | `scdrs munge-gs --weight`. |
| `scdrs.n_max` | `1000` | Top MAGMA genes in the disease gene set. |
| `scdrs.n_ctrl` | `1000` | Matched control gene sets (`compute-score --n-ctrl`). |
| `scdrs.cov_cols` | `[]` | Covariates regressed out by scDRS. Numeric columns are used as they are; categorical ones are one-hot encoded. An intercept is always included. |
| `scdrs.neighbors_reduction` | `null` | KNN graph for the heterogeneity test. `null` uses the scDRS default (PCA with 20 PCs, k = 15, on the scDRS-normalised data). Give a Seurat reduction name (e.g. `harmony`) to use a batch-corrected embedding. |
| `scdrs.n_neighbors` | `15` | k when `neighbors_reduction` is set. |

## `ldsc`

| Key | Default | Description |
|---|---|---|
| `env_activate` | — | Command that activates the LDSC (Python 2) environment for LDSC steps. |
| `software_dir` | — | Folder with `ldsc.py`, `make_annot.py`, `munge_sumstats.py`. |
| `plink_prefix` | — | 1000G EUR PLINK prefix up to the chromosome number (`.../1000G.EUR.QC.`). |
| `baseline_prefix` | — | baselineLD v2.2 prefix (`.../baselineLD.`). |
| `frq_prefix` | — | allele-frequency prefix (`.../1000G.EUR.QC.`). |
| `weights_prefix` | — | regression-weight prefix (`.../weights.hm3_noMHC.`). |
| `hm3_snplist` | — | HapMap3 SNP list (`w_hm3.snplist`); used by munge and for the LD-score SNP set. |
| `gene_coord_file` | `null` | `GENE CHR START END` file; `null` derives it from `magma.gene_loc`. |
| `window_kb` | `100` | Window around each gene for the annotation. |
| `top_pct` | `10` | Top % most specific genes per cell type in each LDSC gene set. |

GWAS columns used by LDSC: `gwas_columns.a1` (effect allele), `a2` (other allele) and `n` (per-SNP N, or `traits.<name>.n`). For `format: vcf`, ALT/REF and `SS` (or `n`) are used.

## `integration`

| Key | Default | Description |
|---|---|---|
| `alpha` | `0.05` | Raw P threshold for "significant" in the integrated table and figure colours. |

## `resources`

`<rule>: {mem_gb, time_min, cpus}`, plus `default` for rules not listed. With the SLURM profile these become `--mem`, `-t` and `-c`. With the local profile, `mem_mb` is limited by the profile's `resources: mem_mb=`.

### Scaling resources

Measured peaks (the whole 8k demo peaked at 13.6 GB and took 75 min on 4 cores):

| Rule | 8k cells (1 chunk) | 20k cells (2 chunks) | 310k cells (5 chunks) |
|---|---|---|---|
| `sc_prep` | < 14 GB | ~11 GB | ~105 GB at load; set 300–500 GB |
| `scpagwas_chunk` | < 14 GB | ~15 GB, 30–50 min | 200–400 GB, 4–24 h |
| `scpagwas_finalize` | < 14 GB | < 15 GB | 400–800 GB |
| `scpagwas_celltype` | < 14 GB | ~8 GB, 30 min | 150–300 GB |
| `export_counts` | < 14 GB | ~2 GB | ~70 GB |
| `scdrs_score` | < 14 GB | ~2 GB | 300–600 GB |

Extra rules: `celltype_specificity` (loads the Seurat object; memory similar to `export_counts`), `magma_covar` (< 4 GB, minutes), `ldsc_l2` (one job per chromosome, ~8 GB), `ldsc_munge` (~16 GB), `ldsc_h2` (one regression per cell type, run `cpus` at a time, ~2–4 GB each).

Rules of thumb: `scpagwas_chunk` memory depends mostly on the GWAS (LD blocks) and only weakly on cell count, while `finalize` and `scdrs_score` grow with cell count. If a job is killed for memory (SLURM state `OUT_OF_MEMORY`), raise its `mem_gb` and re-run the same command; finished steps are kept.

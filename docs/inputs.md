# Preparing the inputs

The pipeline needs three inputs per analysis: one Seurat object, plus a GWAS summary-statistics file and MAGMA gene results (`.genes.out` and `.genes.raw`) for each trait. MAGMA and LDSC also need reference files ([references.md](references.md)).

## 1. Single-cell data: Seurat object

| Requirement | Why | Checked? |
|---|---|---|
| Saved with `saveRDS()`, Seurat v4 or v5 | v5 `Assay5` objects are converted automatically (`JoinLayers` → v3 `Assay`). | — |
| Human genes as **gene symbols** | They must match the scPagwas pathway genes and the MAGMA gene symbols. | — |
| `counts` layer holds **raw integer UMI counts** | scDRS models raw counts (`--flag-raw-count True`). | yes, the step stops if values are not integers |
| `data` layer holds **log-normalised** values (`NormalizeData`, LogNormalize) | scPagwas uses `data`. Set `run_normalize: true` if `data` is empty or equal to counts. | — |
| Cell-type column(s) in `meta.data` | `celltype_col` and `group_cols` in the config | yes (missing columns stop the run) |

Notes:
- Cells with `NA` in `celltype_col` are dropped.
- Extra assays (SCT, integrated) and reductions are dropped inside the pipeline. Your file is never modified.
- `celltype_col` is the annotation used for the scPagwas regression and bootstrap. Use your main cell-type level (e.g. 8–20 types). Finer levels can go in `group_cols`.
- Very small groups (< ~20 cells) give unstable cell-type P values.

Quick check in R:

```r
obj <- readRDS("my_object.rds")
DefaultAssay(obj) <- "RNA"
x <- LayerData(obj, "counts")            # Seurat v5 (v4: GetAssayData(obj, slot = "counts"))
all(x@x[1:1e5] %% 1 == 0)                # TRUE -> raw counts
table(obj$cell_type, useNA = "ifany")     # your annotation column
```

## 2. GWAS summary statistics

### Genome build (important)

Gene coordinates must match the GWAS build. If they don't, SNPs are assigned to the wrong genes and nothing reports an error.

| GWAS build | `scpagwas.block_annotation` | `gwas_filter.exclude_region` (MHC) | MAGMA gene.loc |
|---|---|---|---|
| GRCh37 / hg19 | `block_annotation_hg37` | chr6: 25,000,000–35,000,000 | `NCBI37.3.gene.loc` |
| GRCh38 / hg38 | `block_annotation` | chr6: 28,500,000–33,500,000 | `NCBI38.gene.loc` |

Most PGC, UK Biobank (Neale) and IEU OpenGWAS releases are GRCh37.

### Format `table`: delimited text

Any tab- or whitespace-delimited file, gzipped or not. Lines starting with `##` are skipped. Map your column names under `gwas_columns` (shared default) or `traits.<name>.columns` (per-trait override). Set a field to `null` if the column is absent.

| Field | Meaning | Required |
|---|---|---|
| `chrom`, `pos`, `rsid` | chromosome (with or without `chr`), base position, rsID | yes |
| `beta` *or* `or` | effect size, or odds ratio (converted with `log`) | one of them |
| `se` *or* `p` | standard error, or P value (SE = \|beta\| / z(P/2)) | one of them |
| `maf` *or* `eaf` *or* `fcas`+`fcon`+`ncas`+`ncon` | allele frequency; the case/control version is size-weighted | one of them |
| `a1`, `a2` | effect allele, other allele | for LDSC |
| `n` | per-SNP sample size (or a constant `traits.<name>.n`) | for LDSC |

Examples:

```yaml
# PGC (default mapping)
columns: {}

# GWAS Catalog harmonised
columns: {chrom: hm_chrom, pos: hm_pos, rsid: hm_rsid, beta: hm_beta, se: standard_error,
          p: p_value, eaf: hm_effect_allele_frequency, fcas: null, fcon: null, ncas: null, ncon: null}

# Typical PLINK / LDSC-style with OR
columns: {chrom: CHR, pos: BP, rsid: SNP, beta: null, or: OR, se: SE, p: P,
          eaf: FRQ, fcas: null, fcon: null, ncas: null, ncon: null}
```

### Format `vcf`: IEU OpenGWAS GWAS-VCF

Set `format: vcf`. The pipeline reads `ES` (beta), `SE` and `AF` from the sample column using each row's FORMAT field (e.g. `ES:SE:LP:AF:ID`), and takes the rsID from the `ID` column.

### Filters applied (`gwas_format`)

Filters run in this order, and each step's count is written to `logs/gwas_format/<trait>.log`:
1. rsID present (`rs…`); scPagwas LD (`chrom_ld`) is indexed by rsID
2. autosomes 1–22 (`autosomes_only`)
3. outside `exclude_region` (MHC)
4. valid beta and SE > 0, frequency present
5. MAF > `min_maf`

scPagwas then applies its own `maf_filter` (default 0.1) inside `GWAS_summary_input`.

## 3. MAGMA gene results

Use the outputs of one MAGMA gene analysis of the **same GWAS**:
- `.genes.out` → scDRS disease gene set (columns `GENE`, `P`);
- `.genes.raw` → MAGMA gene-property analysis against cell-type specificity (`--gene-covar`).

Typical MAGMA run (GRCh37, European LD reference):

```bash
# 1) SNP -> gene annotation (window in kb upstream,downstream)
magma --annotate window=50,50 \
      --snp-loc  my_gwas.snploc          `# columns: SNP CHR BP` \
      --gene-loc NCBI37.3.gene.loc \
      --out      my_gwas

# 2) gene analysis (N = sample size column, or N=<number>)
magma --bfile    g1000_eur/g1000_eur \
      --pval     my_gwas.pval use=SNP,P ncol=N \
      --gene-annot my_gwas.genes.annot \
      --out      my_gwas
# -> my_gwas.genes.out and my_gwas.genes.raw
```

Reference files (`NCBI37.3.gene.loc`, `g1000_eur`) are on the [MAGMA website](https://cncr.nl/research/magma/).

Gene IDs: with `magma.gene_id: auto`, the pipeline detects Entrez IDs (all digits) and maps them to symbols with `magma.gene_loc`, keeping the smallest P for duplicated symbols. Files that already use symbols need no gene.loc.

The demo MAGMA file used `window=50,50`, `g1000_eur` and `ncol=N`, with Entrez IDs converted to symbols.

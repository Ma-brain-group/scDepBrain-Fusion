# Methods

This page describes what the pipeline computes, in a form that can be adapted for a manuscript Methods section. Defaults are shown in brackets.

## GWAS processing

Summary statistics were restricted to autosomal variants with an rsID. Variants in the extended MHC region [chr6:25–35 Mb, GRCh37] were removed, as were variants with MAF ≤ [0.01]. For case–control meta-analyses reporting case and control allele frequencies, MAF was computed from the sample-size-weighted allele frequency. When needed, effect sizes were derived from odds ratios (β = log OR) and standard errors from P values (SE = |β| / Φ⁻¹(1 − P/2)).

## Cell-type specificity

For each annotation level, mean normalised expression per cell type was computed with `Seurat::AverageExpression`. Specificity was defined as each cell type's mean divided by the sum over cell types (Bryois et al. 2020). Genes with zero expression in all cell types were removed.

## MAGMA gene-property analysis

Gene-level association statistics from a MAGMA v1.10 gene analysis of the GWAS (`.genes.raw`) were tested against the continuous specificity of each cell type with `magma --gene-results <.genes.raw> --gene-covar <specificity>` (de Leeuw et al. 2015). Gene symbols were mapped to Entrez IDs with the MAGMA gene-location file, and duplicated IDs were averaged. MAGMA's default gene-property model was used, which adjusts for gene-level confounders such as gene size and gene density; see the MAGMA manual for the exact covariates and test direction of your MAGMA version.

## Stratified LD score regression (LDSC-SEG)

For each cell type, the top [10]% most specific genes formed a gene set (Finucane et al. 2018). Together with a control set of all expressed genes, each gene set was extended by ±[100] kb and mapped to 1000 Genomes Phase 3 European SNPs (`make_annot.py`). LD scores were computed with a 1 cM window and printed for HapMap3 SNPs. All annotations of a chromosome were computed in one `ldsc.py --l2` run and split per gene set; this is numerically identical to separate runs. GWAS summary statistics were processed with `munge_sumstats.py` (HapMap3 SNPs). For each cell type, S-LDSC was run with the baseline-LD v2.2 model, the all-gene control and the cell-type annotation (`--overlap-annot`, HapMap3 weights without the MHC). The reported P value is the one-sided test of the cell-type coefficient (coefficient > 0). LDSC was run only for European-ancestry GWAS.

## scPagwas

scPagwas v2.0.0 (Ma et al. 2023) was run with KEGG pathways (`Genes_by_pathway_kegg`), the LD reference shipped with the package (`chrom_ld`), and gene coordinates matching the GWAS build [`block_annotation_hg37`].

1. **Single-cell input.** Log-normalised expression and the main cell-type annotation (`Idents`) were passed to `Single_data_input`. Pathway activity was summarised with `Pathway_pcascore_run` (first principal component of each pathway's genes), and pathway mean expression with `pa_meanexpr`.
2. **SNP-to-gene mapping.** GWAS variants passing `maf_filter` [0.1] were assigned to genes within ±`marg` [10 kb] (`SnpToGene`), then to pathway blocks (`Pathway_annotation_input`).
3. **Single-cell regression.** `Link_pathway_blocks_gwas` (single-cell mode) related each pathway block's GWAS effects, accounting for LD, to each cell's pathway activity. To keep memory manageable this ran on [5] cell chunks, and the per-cell results were concatenated.
4. **Cell scores.** Per-cell genetically associated pathway activation scores (gPAS) were computed as the sum over pathways of regression coefficient × pathway activity, with outliers removed (`Merge_gPas`; identical to `scPagwas_perform_score(remove_outlier = TRUE)`).
5. **Trait-relevant genes and TRS.** The Pearson correlation between gPAS and each gene's expression was computed across all cells (`scGet_PCC`). The top [500] positively correlated genes defined the trait-relevance score (TRS), and the top [500] negatively correlated genes defined a downTRS, both computed with `Seurat::AddModuleScore`.
6. **Per-cell significance.** Each cell's TRS was compared with [100] random control gene sets of [1,000] genes and corrected for gene-set variance (`Get_CorrectBg_p`), giving a pooled per-cell P value.
7. **Cell-type significance.** Per-cell P values were combined within each cell type (`Merge_celltype_p`). This was done for every annotation level, and Benjamini–Hochberg FDR was computed across cell types.
8. **Cell-type regression.** Separately, `Link_pathway_blocks_gwas` in cell-type mode, followed by `Pagwas_perform_regression`, estimated each cell type's association with the trait. `Boot_evaluate` [200 bootstraps, 50% subsampling] gave one-sided bootstrap P values for positive effects.

**Implementation note.** Two scPagwas internals were replaced with numerically equivalent, faster code. In `Get_CorrectBg_p`, gene variances are computed from sparse row sums: var = (Σx² − n·x̄²)/(n − 1). In `get_p_from_empi_null`, empirical P values are computed by sorting the null distribution once and using binary search (`findInterval`). Per-cell P values are matched to cells by name.

## scDRS

scDRS v1.0.2 (Zhang et al. 2022) was run on raw UMI counts exported from the Seurat `counts` layer.

1. **Disease gene set.** MAGMA gene-level P values (de Leeuw et al. 2015) were converted to gene symbols. The top [1,000] genes, weighted by MAGMA z-score, formed the disease gene set (`scdrs munge-gs --weight zscore --n-max 1000`).
2. **Cell scores.** `scdrs compute-score` (`--flag-filter-data True --flag-raw-count True`, [1,000] matched control gene sets) produced per-cell normalised disease scores and Monte Carlo P values.
3. **Cell-type tests.** `scdrs perform-downstream --group-analysis` tested each cell type for association (top-5% score quantile against controls) and heterogeneity (Geary's C on a KNN graph; default: PCA with 20 PCs, k = 15). This was done for each annotation level. BH FDR was computed across cell types.

## Integration

For each trait, annotation level and cell type, the P values of the available methods (MAGMA, LDSC, scPagwas `Merge_celltype_p`, scDRS association) were summarised as the mean of −log10 P. Cell types were labelled by the combination of methods with nominal P < 0.05. This summary is descriptive, not a formal meta-analysis: the methods share the same GWAS and are not independent.

In addition, scPagwas and scDRS cell-type results were merged by cell-type name. Agreement between methods at the single-cell level was measured as the Spearman correlation between the scPagwas TRS and the scDRS normalised score, across all cells and within each cell type.

## Software

MAGMA v1.10, LDSC v1.0.1 (Python 2.7), Snakemake 7.32.4, R 4.3.3 (Seurat 5.4.0, SeuratObject 5.3.0, scPagwas 2.0.0, hdf5r, bigreadr, data.table), Python 3.10 (scdrs 1.0.2, anndata 0.11.4, scanpy 1.10.1, h5py, statsmodels, matplotlib).

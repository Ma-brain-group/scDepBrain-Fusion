# How the demo data were made

`make_demo_subset.R` builds `example/data/example_brain_8k.rds` from an integrated
human brain snRNA-seq atlas (healthy donors; 10 datasets / brain regions):

- 400 nuclei randomly sampled (seed 2026) per subtype (`final_anno`, 20 subtypes) -> 8,000 nuclei (~75% neurons; full atlas 68%)
- raw UMI counts only (`run_normalize: true` in the demo config), genes detected in >= 1 cell
- metadata kept: `anno`, `final_anno` (20 subtypes), `nCount_RNA`, `nFeature_RNA`
- barcodes replaced by `cell_0001` ... ; no donor, sample or batch identifiers

`example/data/magma/MDD_EUR_PGC2025.genes.symbol.out` is the MAGMA v1.10 gene analysis
of the same PGC MDD 2025 EUR GWAS (`--annotate window=50,50`, NCBI37.3 gene locations,
`--bfile g1000_eur`, `ncol=N`), with Entrez IDs replaced by gene symbols.

The script is kept for transparency; it needs the full atlas and cannot be run from this repository.

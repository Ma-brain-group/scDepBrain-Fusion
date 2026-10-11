# External software and reference files (MAGMA, LDSC)

scPagwas and scDRS run with the conda environment alone. Two optional methods need extra downloads:

| Method | You need | Config keys |
|---|---|---|
| MAGMA | MAGMA executable, gene-location file, and `.genes.raw` / `.genes.out` from your own MAGMA gene analysis | `magma.binary`, `magma.gene_loc`, `traits.<name>.magma_genes(_raw)` |
| LDSC | LDSC software (Python 2 environment) and the 1000 Genomes / baselineLD reference panel | `ldsc.*` |

If you can't install one of them, switch it off with `methods: {magma: false}` or `methods: {ldsc: false}`. The integration then uses the methods that ran.

The demo config expects everything under `references/` in the repository. That folder is ignored by git.

```
references/
  magma/
    magma                       MAGMA executable
    NCBI37.3.gene.loc           gene locations, GRCh37
  ldsc/                         clone of the LDSC repository (ldsc.py, make_annot.py, munge_sumstats.py)
  ldsc_ref/
    1000G_EUR_Phase3_plink/     1000G.EUR.QC.{1..22}.{bed,bim,fam}
    baseline_LD_v2.2/           baselineLD.{1..22}.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}
    1000G_Phase3_frq/           1000G.EUR.QC.{1..22}.frq
    1000G_Phase3_weights_hm3_no_MHC/   weights.hm3_noMHC.{1..22}.l2.ldscore.gz
    w_hm3.snplist               HapMap3 SNP list (SNP A1 A2)
```

Symbolic links work, so existing copies elsewhere on your system can be linked in.

## MAGMA

1. Download the MAGMA executable for your platform and the gene-location file matching your GWAS build (`NCBI37.3.gene.loc` for GRCh37) from the MAGMA website: <https://cncr.nl/research/magma/>. Also download a European LD reference such as `g1000_eur` to run the gene analysis.
2. Run the MAGMA gene analysis on your GWAS (see [inputs.md](inputs.md#3-magma-gene-results)). This produces `<prefix>.genes.out` (used by scDRS) and `<prefix>.genes.raw` (used by the MAGMA gene-property analysis in this pipeline). Both must come from the **same** run.

The pipeline then builds a cell-type specificity matrix from your Seurat object and runs

```
magma --gene-results <prefix>.genes.raw --gene-covar <specificity matrix> --out ...
```

once per trait and annotation level.

## LDSC

1. **Software.** Install it in its own Python 2 environment, as described in the LDSC repository:
   ```bash
   git clone https://github.com/bulik/ldsc.git references/ldsc
   conda env create -f references/ldsc/environment.yml      # creates env "ldsc"
   ```
   Set `ldsc.env_activate` to the command that activates this environment, e.g. `source ~/miniconda3/etc/profile.d/conda.sh && conda activate ldsc`. Pipeline jobs switch to it only for the LDSC steps.
2. **Reference panel (GRCh37, European).** These are the standard S-LDSC / LDSC-SEG resources from the Price lab, linked from the LDSC wiki (<https://github.com/bulik/ldsc/wiki>):
   - 1000 Genomes Phase 3 EUR PLINK files (`1000G_EUR_Phase3_plink`)
   - baseline-LD model v2.2 (`baseline_LD_v2.2`)
   - allele frequencies (`1000G_Phase3_frq`)
   - regression weights without MHC (`1000G_Phase3_weights_hm3_no_MHC`)
   - HapMap3 SNP list (`w_hm3.snplist`)
3. **Gene coordinates.** By default they are taken from `magma.gene_loc`. To use another annotation (e.g. GENCODE), set `ldsc.gene_coord_file` to a tab-separated file with columns `GENE CHR START END`.

LDSC models LD with a European reference panel. Run it only for GWAS of European ancestry, and set `ldsc: false` for other or trans-ancestry traits.

### What the pipeline does with them

1. The top 10% most specific genes per cell type (`ldsc.top_pct`) become a gene set, together with an "all expressed genes" control set.
2. Each gene set is extended by ±100 kb (`ldsc.window_kb`), annotated onto 1000G SNPs and turned into LD scores (1 cM window, HapMap3 SNPs). This runs as one job per chromosome that computes all gene sets at once.
3. S-LDSC runs for each cell type with `baselineLD + all_gene + cell type` and `--overlap-annot`. The reported P is the one-sided test of the cell-type coefficient.

## Checking your setup

```bash
bash run_pipeline.sh config/my_config.yaml local -n
```

A dry run lists every missing file, e.g. `ldsc.baseline_prefix: expected file .../baselineLD.22.l2.ldscore.gz`.

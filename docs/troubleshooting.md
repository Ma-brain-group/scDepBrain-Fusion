# Troubleshooting

When a step fails, Snakemake prints the rule name and its log path. Open `<outdir>/logs/<step>/<trait>.log` first. For SLURM, also check `logs/slurm/<rule>_*.err` in the repository folder. Fix the problem and re-run the same command; finished steps are not repeated.

## Configuration and inputs

| Message | Cause / fix |
|---|---|
| `No config given: use --configfile ...` | Call through `run_pipeline.sh <config> <mode>`, or pass `--configfile`. |
| `Config errors: ... not found` | A path in the config is wrong. Relative paths resolve against the repository folder. |
| `celltype_col ... must be in group_cols` | Add the main annotation to `group_cols`. |
| `Columns not found in GWAS header: [...]` | Adjust `gwas_columns` or `traits.<name>.columns` to your header (printed in the log). |
| `No SNPs survived filtering` | Usually a wrong column mapping, non-rsID variant IDs (e.g. `1:12345:A:G`), or a frequency column holding something else. Check the per-filter counts in the log. |
| `Metadata columns missing from Seurat object` | `celltype_col` / `group_cols` are not in `meta.data`. |
| `counts layer does not look like raw integer counts` | The `counts` layer holds normalised values. Rebuild the object with raw counts, which scDRS needs. |
| `MAGMA file lacks column GENE/P` | Pass the `.genes.out` file, not `.genes.raw` or `.gsa.out`. |
| Many `Unmapped Entrez IDs dropped` | The `gene_loc` build or version doesn't match the MAGMA run. |

## Running

| Symptom | Cause / fix |
|---|---|
| `Directory cannot be locked` | Another Snakemake is running in this folder, or a previous run was killed. If none is running: `snakemake --configfile <cfg> --unlock`. |
| SLURM job state `OUT_OF_MEMORY` / `TIMEOUT` | Raise `resources.<rule>.mem_gb` / `time_min` and re-run. See [parameters.md](parameters.md#scaling-resources). |
| Jobs stay pending with `QOSMaxMemoryPerUser` | Your cluster caps total memory per user, so jobs start as memory frees up. Lower `jobs:` in `profiles/slurm/config.yaml` or raise `n_chunks`. |
| `sbatch: error: ... partition` / `account` | Set `export SBATCH_PARTITION=...` (and `SBATCH_ACCOUNT=...`) before `run_pipeline.sh`. |
| Jobs fail with `Rscript: command not found` on nodes | The environment is not exported to jobs. Set `env_activate` in the config. |
| `MissingRuleException` when giving a target | Targets must be the absolute paths shown in the dry run. |

## MAGMA / LDSC

| Message | Cause / fix |
|---|---|
| `traits.X.magma_genes_raw ... not found` | Give the `.genes.raw` file of your MAGMA gene analysis, or set `methods: {magma: false}`. |
| `magma.binary not found` | Put MAGMA on `PATH` or give the full path; see [references.md](references.md). |
| `ldsc.baseline_prefix: expected file ...` | A prefix must end right before the chromosome number, e.g. `.../baselineLD.`. |
| LDSC step: `ImportError` / `SyntaxError` in `ldsc.py` | LDSC runs under Python 2: check `ldsc.env_activate`. |
| `LDSC needs gwas_columns.a1` / `a sample size` | Map the allele and N columns, or set `traits.<name>.n`. |
| `Could not match AverageExpression columns` | Cell-type names collide after removing non-alphanumeric characters; rename them in `meta.data`. |
| All LDSC P near 0.5 / tiny h² | Check that the GWAS and the references use the same build and that most SNPs survived munging (`logs/ldsc/<trait>_munge.log`). |

## Tool-specific

| Message | Cause / fix |
|---|---|
| `AttributeError: np.float_ was removed in the NumPy 2.0 release` | scdrs was called directly. The pipeline runs scdrs through `scripts/scdrs_cli.py`, which restores the alias. |
| `connectivities not found in adata.obsp; run sc.pp.neighbors first` | Informational: scdrs builds its default KNN graph, and the run continues. |
| scPagwas: `there is no package called 'SOAR'` | Run `Rscript envs/install_R_packages.R`. |
| scPagwas very slow in `Get_CorrectBg_p` | Make sure `scripts/common.R` is sourced (the log should show `scPagwas speed patches applied`). |
| Strong signal in unexpected cell types | Check the genome build (`block_annotation_hg37` vs `block_annotation`) and that the MAGMA run used the same GWAS. |

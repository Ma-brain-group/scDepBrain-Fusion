# =====================================================================
# scDepBrain-Fusion: four-method GWAS x single-cell cell-type association
#   GWAS summary statistics + MAGMA gene results + Seurat object
#     -> MAGMA gene-property, S-LDSC (LDSC-SEG), scPagwas, scDRS
#     -> four-method integration: mean -log10(P) per cell type + figure
#
#   bash run_pipeline.sh config/config_example.yaml local     # run locally
#   bash run_pipeline.sh config/my_config.yaml slurm          # run on SLURM
#   bash run_pipeline.sh config/my_config.yaml local -n       # dry run
# See README.md and docs/ for details.
# =====================================================================
import os
import re
import shutil

# Config is passed with --configfile (run_pipeline.sh defaults to config/config.yaml).
# No `configfile:` directive here, otherwise a user config would be MERGED with
# the example config and inherit its traits.
if not config:
    raise ValueError("No config given: use --configfile config/<your>.yaml")

BASE = workflow.basedir
SCRIPTS = os.path.join(BASE, "scripts")
COMMON_R = os.path.join(SCRIPTS, "common.R")
SCDRS = f"python {os.path.join(SCRIPTS, 'scdrs_cli.py')}"


def abspath(p):
    if p is None:
        return None
    p = os.path.expanduser(str(p))
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(BASE, p))


OUT = abspath(config["outdir"])
SC = config["single_cell"]
PG = config["scpagwas"]
SD = config["scdrs"]

TRAITS = list(config["traits"].keys())
GROUPS = list(SC["group_cols"])
CELLTYPE = SC["celltype_col"]
MT = config.get("methods", {})
RUN_MAGMA = bool(MT.get("magma", False))
RUN_LDSC = bool(MT.get("ldsc", False))
MG = config["magma"]
LD = config.get("ldsc", {}) or {}
LDSC_TRAITS = [t for t in config["traits"] if RUN_LDSC and config["traits"][t].get("ldsc", True)]
CHROMS = list(range(1, 23))
N_CHUNKS = int(PG["n_chunks"])
CHUNKS = list(range(1, N_CHUNKS + 1))


# ---------------------------------------------------------------------
# Input validation (fails at DAG-build time with a readable message)
# ---------------------------------------------------------------------
def _check_config():
    errs = []
    if not os.path.exists(abspath(SC["rds"])):
        errs.append(f"single_cell.rds not found: {abspath(SC['rds'])}")
    if CELLTYPE not in GROUPS:
        errs.append(f"single_cell.celltype_col '{CELLTYPE}' must be in group_cols {GROUPS}")
    if not TRAITS:
        errs.append("no traits defined")
    for t, tc in config["traits"].items():
        if not re.fullmatch(r"[A-Za-z0-9_]+", t):
            errs.append(f"trait name '{t}' must match [A-Za-z0-9_]+")
        for key in ("sumstats", "magma_genes"):
            p = abspath(tc.get(key))
            if p is None or not os.path.exists(p):
                errs.append(f"traits.{t}.{key} not found: {p}")
    if RUN_MAGMA or RUN_LDSC:
        if not MG.get("gene_loc") or not os.path.exists(abspath(MG["gene_loc"])):
            errs.append(f"magma.gene_loc not found (needed for MAGMA / LDSC): {MG.get('gene_loc')}")
    if RUN_MAGMA:
        b = MG.get("binary", "magma")
        if not (os.path.exists(abspath(b)) or shutil.which(b)):
            errs.append(f"magma.binary not found: {b}")
        for t, tc in config["traits"].items():
            p_ = abspath(tc.get("magma_genes_raw"))
            if p_ is None or not os.path.exists(p_):
                errs.append(f"traits.{t}.magma_genes_raw (.genes.raw, needed for MAGMA gene-covar) not found: {p_}")
    if RUN_LDSC:
        if not os.path.exists(os.path.join(abspath(LD.get("software_dir", "")), "ldsc.py")):
            errs.append(f"ldsc.software_dir has no ldsc.py: {LD.get('software_dir')}")
        for k, suffix in [("plink_prefix", "22.bim"), ("baseline_prefix", "22.l2.ldscore.gz"),
                          ("frq_prefix", "22.frq"), ("weights_prefix", "22.l2.ldscore.gz")]:
            f = str(abspath(LD.get(k, ""))) + suffix
            if not os.path.exists(f):
                errs.append(f"ldsc.{k}: expected file {f}")
        if not os.path.exists(abspath(LD.get("hm3_snplist", ""))):
            errs.append(f"ldsc.hm3_snplist not found: {LD.get('hm3_snplist')}")
    if errs:
        raise ValueError("Config errors:\n  - " + "\n  - ".join(errs))


_check_config()

wildcard_constraints:
    trait="|".join(TRAITS),
    group="|".join(map(re.escape, GROUPS)),
    chunk=r"\d+",
    chr=r"\d+",


def res(rule_name):
    d = dict(config["resources"]["default"])
    d.update(config["resources"].get(rule_name, {}))
    return d


# Per-trait scPagwas overrides: traits.<name>.scpagwas.<key> beats scpagwas.<key>
PER_TRAIT_SCPAGWAS = {"maf_filter", "marg", "pcc_method", "n_topgenes", "bg_n_topgenes",
                      "iters_singlecell", "corr_nrandom", "corr_nselect", "seed"}


def pgt(trait, key):
    over = config["traits"][trait].get("scpagwas") or {}
    bad = set(over) - PER_TRAIT_SCPAGWAS
    if bad:
        raise ValueError(f"traits.{trait}.scpagwas: keys {sorted(bad)} cannot be set per trait "
                         f"(allowed: {sorted(PER_TRAIT_SCPAGWAS)})")
    return over.get(key, PG.get(key))


def gwas_columns(trait):
    cols = dict(config["gwas_columns"])
    cols.update(config["traits"][trait].get("columns") or {})
    return cols


ENV = config.get("env_activate", "")
LDSC_ENV = LD.get("env_activate", "")
LDSC_PRE = f"set +u; {LDSC_ENV}; set -u; " if LDSC_ENV else ""
shell.prefix(f"set -eo pipefail; {ENV}; set -u; " if ENV else "set -euo pipefail; ")

# ---------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------
SCPAGWAS_BOOT = (
    expand(f"{OUT}/scpagwas/{{trait}}/celltype_bootstrap_{CELLTYPE}.csv", trait=TRAITS)
    if PG["run_celltype_bootstrap"] else []
)

rule all:
    input:
        expand(f"{OUT}/scpagwas/{{trait}}/celltype_p_{{group}}.csv", trait=TRAITS, group=GROUPS),
        SCPAGWAS_BOOT,
        expand(f"{OUT}/scdrs/downstream/{{trait}}.scdrs_group.{{group}}", trait=TRAITS, group=GROUPS),
        expand(f"{OUT}/summary/{{trait}}/summary_{{group}}.csv", trait=TRAITS, group=GROUPS),
        f"{OUT}/summary/all_traits_celltype_summary.csv",
        expand(f"{OUT}/integrated/four_method_{{group}}.csv", group=GROUPS),
        f"{OUT}/integrated/four_method_all.csv",


# =====================================================================
# Part A — GWAS formatting
# =====================================================================
rule gwas_format:
    input:
        lambda wc: abspath(config["traits"][wc.trait]["sumstats"]),
    output:
        f"{OUT}/gwas/{{trait}}.scpagwas_input.txt",
    params:
        columns=lambda wc: gwas_columns(wc.trait),
        fmt=lambda wc: config["traits"][wc.trait].get("format", "table"),
        filt=config["gwas_filter"],
    log:
        f"{OUT}/logs/gwas_format/{{trait}}.log",
    threads: res("gwas_format")["cpus"]
    resources:
        mem_mb=res("gwas_format")["mem_gb"] * 1024,
        runtime=res("gwas_format")["time_min"],
    script:
        "scripts/01_gwas_to_scpagwas.py"


# =====================================================================
# Part B — scPagwas
# =====================================================================
rule sc_prep:
    """GWAS-independent: Single_data_input + pathway PCA, then split cells into chunks."""
    input:
        rds=abspath(SC["rds"]),
    output:
        preprocessed=f"{OUT}/scpagwas/shared/pagwas_preprocessed.RData",
        chunks=expand(f"{OUT}/scpagwas/shared/chunk_{{chunk}}_pagwas.RData", chunk=CHUNKS),
    params:
        common=COMMON_R,
        assay=SC["assay"],
        celltype_col=CELLTYPE,
        group_cols=GROUPS,
        run_normalize=SC["run_normalize"],
        pathway_list=PG["pathway_list"],
    log:
        f"{OUT}/logs/scpagwas/sc_prep.log",
    threads: res("sc_prep")["cpus"]
    resources:
        mem_mb=res("sc_prep")["mem_gb"] * 1024,
        runtime=res("sc_prep")["time_min"],
    script:
        "scripts/02_scpagwas_sc_prep.R"


rule scpagwas_gwas_prep:
    input:
        f"{OUT}/gwas/{{trait}}.scpagwas_input.txt",
    output:
        f"{OUT}/scpagwas/{{trait}}/gwas_pagwas.RData",
    params:
        common=COMMON_R,
        maf_filter=lambda wc: pgt(wc.trait, "maf_filter"),
        marg=lambda wc: pgt(wc.trait, "marg"),
        pathway_list=PG["pathway_list"],
        block_annotation=PG["block_annotation"],
    log:
        f"{OUT}/logs/scpagwas/{{trait}}_gwas_prep.log",
    threads: res("scpagwas_gwas_prep")["cpus"]
    resources:
        mem_mb=res("scpagwas_gwas_prep")["mem_gb"] * 1024,
        runtime=res("scpagwas_gwas_prep")["time_min"],
    script:
        "scripts/03_scpagwas_gwas_prep.R"


rule scpagwas_chunk:
    """Heavy step: Link_pathway_blocks_gwas (single-cell mode) on one chunk of cells."""
    input:
        gwas=f"{OUT}/scpagwas/{{trait}}/gwas_pagwas.RData",
        chunk=f"{OUT}/scpagwas/shared/chunk_{{chunk}}_pagwas.RData",
    output:
        f"{OUT}/scpagwas/{{trait}}/chunks/chunk_{{chunk}}_results.RData",
    params:
        common=COMMON_R,
        block_annotation=PG["block_annotation"],
        chrom_ld=PG["chrom_ld"],
        backingpath=lambda wc: f"{OUT}/scpagwas/{wc.trait}/chunks/temp_{wc.chunk}",
    log:
        f"{OUT}/logs/scpagwas/{{trait}}_chunk_{{chunk}}.log",
    threads: res("scpagwas_chunk")["cpus"]
    resources:
        mem_mb=res("scpagwas_chunk")["mem_gb"] * 1024,
        runtime=res("scpagwas_chunk")["time_min"],
    script:
        "scripts/04_scpagwas_chunk.R"


rule scpagwas_finalize:
    """Merge chunks -> gPAS, TRS, background-corrected per-cell P, Merge_celltype_p per group."""
    input:
        preprocessed=f"{OUT}/scpagwas/shared/pagwas_preprocessed.RData",
        chunks=expand(f"{OUT}/scpagwas/{{{{trait}}}}/chunks/chunk_{{chunk}}_results.RData", chunk=CHUNKS),
    output:
        per_cell=f"{OUT}/scpagwas/{{trait}}/per_cell_scores.csv.gz",
        topgenes=f"{OUT}/scpagwas/{{trait}}/scPagwas_topgenes.csv",
        celltype=expand(f"{OUT}/scpagwas/{{{{trait}}}}/celltype_p_{{group}}.csv", group=GROUPS),
    params:
        common=COMMON_R,
        assay=SC["assay"],
        group_cols=GROUPS,
        n_topgenes=lambda wc: pgt(wc.trait, "n_topgenes"),
        pcc_method=lambda wc: pgt(wc.trait, "pcc_method") or "scGet_PCC",
        bg_n_topgenes=lambda wc: pgt(wc.trait, "bg_n_topgenes"),
        iters_singlecell=lambda wc: pgt(wc.trait, "iters_singlecell"),
        corr_nrandom=lambda wc: pgt(wc.trait, "corr_nrandom"),
        corr_nselect=lambda wc: pgt(wc.trait, "corr_nselect"),
        seed=lambda wc: pgt(wc.trait, "seed"),
        save_seurat=PG["save_seurat"],
        seurat_out=lambda wc: f"{OUT}/scpagwas/{wc.trait}/scPagwas_seurat.rds",
    log:
        f"{OUT}/logs/scpagwas/{{trait}}_finalize.log",
    threads: res("scpagwas_finalize")["cpus"]
    resources:
        mem_mb=res("scpagwas_finalize")["mem_gb"] * 1024,
        runtime=res("scpagwas_finalize")["time_min"],
    script:
        "scripts/05_scpagwas_finalize.R"


rule scpagwas_celltype:
    """Cell-type-level regression (Pagwas_perform_regression) + Boot_evaluate."""
    input:
        preprocessed=f"{OUT}/scpagwas/shared/pagwas_preprocessed.RData",
        gwas=f"{OUT}/scpagwas/{{trait}}/gwas_pagwas.RData",
    output:
        f"{OUT}/scpagwas/{{trait}}/celltype_bootstrap_{CELLTYPE}.csv",
    params:
        common=COMMON_R,
        block_annotation=PG["block_annotation"],
        chrom_ld=PG["chrom_ld"],
        bootstrap_iters=PG["bootstrap_iters"],
        bootstrap_part=PG["bootstrap_part"],
        backingpath=lambda wc: f"{OUT}/scpagwas/{wc.trait}/celltype_temp",
    log:
        f"{OUT}/logs/scpagwas/{{trait}}_celltype.log",
    threads: res("scpagwas_celltype")["cpus"]
    resources:
        mem_mb=res("scpagwas_celltype")["mem_gb"] * 1024,
        runtime=res("scpagwas_celltype")["time_min"],
    script:
        "scripts/06_scpagwas_celltype.R"


# =====================================================================
# Part C — scDRS
# =====================================================================
rule export_counts:
    """Write RAW counts (CSC arrays, HDF5) + metadata from the Seurat object."""
    input:
        rds=abspath(SC["rds"]),
    output:
        counts=temp(f"{OUT}/scdrs/data/counts.h5"),
        meta=f"{OUT}/scdrs/data/cell_meta.csv.gz",
        embedding=f"{OUT}/scdrs/data/embedding.csv.gz",
    params:
        common=COMMON_R,
        reduction=SD.get("neighbors_reduction"),
        assay=SC["assay"],
        celltype_col=CELLTYPE,
        keep_cols=list(dict.fromkeys(GROUPS + list(SC.get("keep_meta") or []) + list(SD.get("cov_cols") or []))),
    log:
        f"{OUT}/logs/scdrs/export_counts.log",
    threads: res("export_counts")["cpus"]
    resources:
        mem_mb=res("export_counts")["mem_gb"] * 1024,
        runtime=res("export_counts")["time_min"],
    script:
        "scripts/07_export_counts.R"


rule build_h5ad:
    input:
        counts=f"{OUT}/scdrs/data/counts.h5",
        meta=f"{OUT}/scdrs/data/cell_meta.csv.gz",
        embedding=f"{OUT}/scdrs/data/embedding.csv.gz",
    output:
        h5ad=f"{OUT}/scdrs/data/sc_counts.h5ad",
        cov=f"{OUT}/scdrs/data/cov.tsv",
    params:
        group_cols=GROUPS,
        cov_cols=list(SD.get("cov_cols") or []),
        n_neighbors=SD.get("n_neighbors", 15),
    log:
        f"{OUT}/logs/scdrs/build_h5ad.log",
    threads: res("build_h5ad")["cpus"]
    resources:
        mem_mb=res("build_h5ad")["mem_gb"] * 1024,
        runtime=res("build_h5ad")["time_min"],
    script:
        "scripts/08_build_h5ad.py"


rule magma_to_pval:
    input:
        lambda wc: abspath(config["traits"][wc.trait]["magma_genes"]),
    output:
        f"{OUT}/scdrs/gs/{{trait}}.gene_pval.tsv",
    params:
        gene_id=config["magma"]["gene_id"],
        gene_loc=abspath(config["magma"]["gene_loc"]),
    log:
        f"{OUT}/logs/scdrs/{{trait}}_magma_to_pval.log",
    resources:
        mem_mb=res("default")["mem_gb"] * 1024,
        runtime=res("default")["time_min"],
    script:
        "scripts/09_magma_to_pval.py"


rule scdrs_munge_gs:
    input:
        f"{OUT}/scdrs/gs/{{trait}}.gene_pval.tsv",
    output:
        f"{OUT}/scdrs/gs/{{trait}}.gs",
    params:
        weight=SD["weight"],
        n_max=SD["n_max"],
    log:
        f"{OUT}/logs/scdrs/{{trait}}_munge_gs.log",
    resources:
        mem_mb=res("default")["mem_gb"] * 1024,
        runtime=res("default")["time_min"],
    shell:
        SCDRS + " munge-gs --out-file {output} --pval-file {input} "
        "--weight {params.weight} --n-max {params.n_max} > {log} 2>&1"


rule scdrs_score:
    input:
        h5ad=f"{OUT}/scdrs/data/sc_counts.h5ad",
        gs=f"{OUT}/scdrs/gs/{{trait}}.gs",
        cov=f"{OUT}/scdrs/data/cov.tsv",
    output:
        score=f"{OUT}/scdrs/scores/{{trait}}.score.gz",
        full=f"{OUT}/scdrs/scores/{{trait}}.full_score.gz",
    params:
        outdir=f"{OUT}/scdrs/scores",
        n_ctrl=SD["n_ctrl"],
        cov_opt=lambda wc, input: f"--cov-file {input.cov}" if SD.get("cov_cols") else "",
    log:
        f"{OUT}/logs/scdrs/{{trait}}_score.log",
    threads: res("scdrs_score")["cpus"]
    resources:
        mem_mb=res("scdrs_score")["mem_gb"] * 1024,
        runtime=res("scdrs_score")["time_min"],
    shell:
        SCDRS + " compute-score --h5ad-file {input.h5ad} --h5ad-species human "
        "--gs-file {input.gs} --gs-species human --out-folder {params.outdir} "
        "{params.cov_opt} --flag-filter-data True --flag-raw-count True "
        "--n-ctrl {params.n_ctrl} > {log} 2>&1"


rule scdrs_downstream:
    input:
        h5ad=f"{OUT}/scdrs/data/sc_counts.h5ad",
        full=f"{OUT}/scdrs/scores/{{trait}}.full_score.gz",
    output:
        expand(f"{OUT}/scdrs/downstream/{{{{trait}}}}.scdrs_group.{{group}}", group=GROUPS),
    params:
        outdir=f"{OUT}/scdrs/downstream",
        groups=",".join(GROUPS),
    log:
        f"{OUT}/logs/scdrs/{{trait}}_downstream.log",
    threads: res("scdrs_downstream")["cpus"]
    resources:
        mem_mb=res("scdrs_downstream")["mem_gb"] * 1024,
        runtime=res("scdrs_downstream")["time_min"],
    shell:
        SCDRS + " perform-downstream --h5ad-file {input.h5ad} --score-file {input.full} "
        "--out-folder {params.outdir} --group-analysis {params.groups} "
        "--flag-filter-data True --flag-raw-count True > {log} 2>&1"


# =====================================================================
# Part D — MAGMA gene-property analysis and S-LDSC (LDSC-SEG)
# =====================================================================
SPEC = f"{OUT}/specificity"

rule celltype_specificity:
    """Bryois specificity per annotation level -> MAGMA covariates + LDSC gene sets."""
    input:
        rds=abspath(SC["rds"]),
    output:
        covar=expand(f"{SPEC}/magma_covar_{{group}}.txt", group=GROUPS),
        genesets=directory(f"{SPEC}/ldsc_genesets"),
        names=f"{SPEC}/celltype_names.csv",
    params:
        common=COMMON_R,
        assay=SC["assay"],
        celltype_col=CELLTYPE,
        group_cols=GROUPS,
        run_normalize=SC["run_normalize"],
        gene_loc=abspath(MG.get("gene_loc")),
        top_pct=LD.get("top_pct", 10),
    log:
        f"{OUT}/logs/specificity/celltype_specificity.log",
    threads: res("celltype_specificity")["cpus"]
    resources:
        mem_mb=res("celltype_specificity")["mem_gb"] * 1024,
        runtime=res("celltype_specificity")["time_min"],
    script:
        "scripts/11_celltype_specificity.R"


rule magma_covar:
    input:
        raw=lambda wc: abspath(config["traits"][wc.trait].get("magma_genes_raw")),
        covar=f"{SPEC}/magma_covar_{{group}}.txt",
    output:
        f"{OUT}/magma/{{trait}}_{{group}}.gsa.out",
    params:
        magma=MG.get("binary", "magma"),
        prefix=lambda wc: f"{OUT}/magma/{wc.trait}_{wc.group}",
    log:
        f"{OUT}/logs/magma/{{trait}}_{{group}}.log",
    resources:
        mem_mb=res("magma_covar")["mem_gb"] * 1024,
        runtime=res("magma_covar")["time_min"],
    shell:
        "{params.magma} --gene-results {input.raw} --gene-covar {input.covar} "
        "--out {params.prefix} > {log} 2>&1"


rule ldsc_coords:
    """Gene coordinates for make_annot.py (GENE CHR START END) from the MAGMA gene.loc."""
    input:
        abspath(LD.get("gene_coord_file") or MG.get("gene_loc") or "MISSING"),
    output:
        f"{OUT}/ldsc/gene_coords.tsv",
    params:
        from_magma=not bool(LD.get("gene_coord_file")),
    run:
        import pandas as pd
        if params.from_magma:
            d = pd.read_csv(input[0], sep=r"\s+", header=None, usecols=[1, 2, 3, 5],
                            names=["CHR", "START", "END", "GENE"], dtype=str)
            d = d[["GENE", "CHR", "START", "END"]].drop_duplicates("GENE")
            d = d[d.CHR.isin([str(i) for i in range(1, 23)])]
        else:
            d = pd.read_csv(input[0], sep=r"\s+", dtype=str)[["GENE", "CHR", "START", "END"]]
        d.to_csv(output[0], sep="\t", index=False)


rule ldsc_print_snps:
    input:
        abspath(LD.get("hm3_snplist") or "MISSING"),
    output:
        f"{OUT}/ldsc/hm3.print_snps.txt",
    shell:
        "awk 'NR > 1 {{print $1}}' {input} > {output}"


rule ldsc_l2_chr:
    """All gene sets on one chromosome: make_annot + one multi-column ldsc --l2, then split."""
    input:
        genesets=f"{SPEC}/ldsc_genesets",
        coords=f"{OUT}/ldsc/gene_coords.tsv",
        snps=f"{OUT}/ldsc/hm3.print_snps.txt",
    output:
        touch(f"{OUT}/ldsc/ldscores/.chr{{chr}}.done"),
    params:
        plink=abspath(LD.get("plink_prefix", "")),
        window=int(LD.get("window_kb", 100)) * 1000,
        ldsc=abspath(LD.get("software_dir", "")),
        outdir=f"{OUT}/ldsc/ldscores",
    log:
        f"{OUT}/logs/ldsc/l2_chr{{chr}}.log",
    resources:
        mem_mb=res("ldsc_l2")["mem_gb"] * 1024,
        runtime=res("ldsc_l2")["time_min"],
    shell:
        LDSC_PRE + "bash " + os.path.join(SCRIPTS, "12_ldsc_l2_chr.sh") +
        " {wildcards.chr} {input.genesets} {input.coords} {params.plink} {input.snps} "
        "{params.window} {params.ldsc} {params.outdir} " + SCRIPTS + " > {log} 2>&1"


rule ldsc_prep:
    """GWAS -> SNP A1 A2 N BETA P for munge_sumstats.py."""
    input:
        lambda wc: abspath(config["traits"][wc.trait]["sumstats"]),
    output:
        temp(f"{OUT}/ldsc/sumstats/{{trait}}.pre.tsv.gz"),
    params:
        columns=lambda wc: gwas_columns(wc.trait),
        fmt=lambda wc: config["traits"][wc.trait].get("format", "table"),
        n=lambda wc: config["traits"][wc.trait].get("n"),
    log:
        f"{OUT}/logs/ldsc/{{trait}}_prep.log",
    resources:
        mem_mb=res("ldsc_munge")["mem_gb"] * 1024,
        runtime=res("ldsc_munge")["time_min"],
    script:
        "scripts/13_ldsc_prep_sumstats.py"


rule ldsc_munge:
    input:
        f"{OUT}/ldsc/sumstats/{{trait}}.pre.tsv.gz",
    output:
        f"{OUT}/ldsc/sumstats/{{trait}}.sumstats.gz",
    params:
        ldsc=abspath(LD.get("software_dir", "")),
        hm3=abspath(LD.get("hm3_snplist") or "MISSING"),
        prefix=lambda wc: f"{OUT}/ldsc/sumstats/{wc.trait}",
    log:
        f"{OUT}/logs/ldsc/{{trait}}_munge.log",
    resources:
        mem_mb=res("ldsc_munge")["mem_gb"] * 1024,
        runtime=res("ldsc_munge")["time_min"],
    shell:
        LDSC_PRE + "python {params.ldsc}/munge_sumstats.py --sumstats {input} "
        "--snp SNP --a1 A1 --a2 A2 --N-col N --signed-sumstats BETA,0 --p P "
        "--merge-alleles {params.hm3} --chunksize 500000 --out {params.prefix} > {log} 2>&1"


rule ldsc_h2:
    input:
        done=expand(f"{OUT}/ldsc/ldscores/.chr{{chr}}.done", chr=CHROMS),
        sumstats=f"{OUT}/ldsc/sumstats/{{trait}}.sumstats.gz",
    output:
        touch(f"{OUT}/ldsc/{{trait}}/.h2.done"),
    params:
        ldscores=f"{OUT}/ldsc/ldscores",
        baseline=abspath(LD.get("baseline_prefix", "")),
        frq=abspath(LD.get("frq_prefix", "")),
        weights=abspath(LD.get("weights_prefix", "")),
        ldsc=abspath(LD.get("software_dir", "")),
        outdir=lambda wc: f"{OUT}/ldsc/{wc.trait}/results",
    log:
        f"{OUT}/logs/ldsc/{{trait}}_h2.log",
    threads: res("ldsc_h2")["cpus"]
    resources:
        mem_mb=res("ldsc_h2")["mem_gb"] * 1024,
        runtime=res("ldsc_h2")["time_min"],
    shell:
        LDSC_PRE + "bash " + os.path.join(SCRIPTS, "14_ldsc_h2.sh") +
        " {input.sumstats} {params.ldscores} {params.baseline} {params.frq} {params.weights} "
        "{params.ldsc} {params.outdir} {threads} > {log} 2>&1"


rule ldsc_collect:
    input:
        done=f"{OUT}/ldsc/{{trait}}/.h2.done",
        names=f"{SPEC}/celltype_names.csv",
    output:
        f"{OUT}/ldsc/{{trait}}/ldsc_celltype_results.csv",
    params:
        res_dir=lambda wc: f"{OUT}/ldsc/{wc.trait}/results",
    log:
        f"{OUT}/logs/ldsc/{{trait}}_collect.log",
    script:
        "scripts/15_ldsc_collect.py"


# =====================================================================
# Part E — Four-method integration
# =====================================================================
def _magma_out(t, g):
    return f"{OUT}/magma/{t}_{g}.gsa.out" if RUN_MAGMA else ""


def _ldsc_out(t):
    return f"{OUT}/ldsc/{t}/ldsc_celltype_results.csv" if t in LDSC_TRAITS else ""


rule integrate:
    input:
        names=f"{SPEC}/celltype_names.csv",
        magma=lambda wc: [_magma_out(t, wc.group) for t in TRAITS] if RUN_MAGMA else [],
        ldsc=[_ldsc_out(t) for t in LDSC_TRAITS],
        scpagwas=expand(f"{OUT}/scpagwas/{{trait}}/celltype_p_{{{{group}}}}.csv", trait=TRAITS),
        scdrs=expand(f"{OUT}/scdrs/downstream/{{trait}}.scdrs_group.{{{{group}}}}", trait=TRAITS),
    output:
        table=f"{OUT}/integrated/four_method_{{group}}.csv",
        png=f"{OUT}/integrated/four_method_{{group}}.png",
        pdf=f"{OUT}/integrated/four_method_{{group}}.pdf",
    params:
        group=lambda wc: wc.group,
        traits=TRAITS,
        alpha=config.get("integration", {}).get("alpha", 0.05),
        magma=lambda wc: {t: _magma_out(t, wc.group) for t in TRAITS},
        ldsc={t: _ldsc_out(t) for t in TRAITS},
        scpagwas=lambda wc: {t: f"{OUT}/scpagwas/{t}/celltype_p_{wc.group}.csv" for t in TRAITS},
        scdrs=lambda wc: {t: f"{OUT}/scdrs/downstream/{t}.scdrs_group.{wc.group}" for t in TRAITS},
    log:
        f"{OUT}/logs/integrated/{{group}}.log",
    resources:
        mem_mb=res("summarize")["mem_gb"] * 1024,
        runtime=res("summarize")["time_min"],
    script:
        "scripts/16_integrate_four_methods.R"


rule integrate_all:
    input:
        expand(f"{OUT}/integrated/four_method_{{group}}.csv", group=GROUPS),
    output:
        f"{OUT}/integrated/four_method_all.csv",
    run:
        import pandas as pd
        frames = []
        for g, f in zip(GROUPS, input):
            d = pd.read_csv(f)
            d.insert(0, "annotation", g)
            frames.append(d)
        pd.concat(frames, ignore_index=True).to_csv(output[0], index=False)


# =====================================================================
# Part F — Per-trait scPagwas / scDRS summary (heatmap, per-cell agreement)
# =====================================================================
rule summarize_trait:
    input:
        scp_ct=expand(f"{OUT}/scpagwas/{{{{trait}}}}/celltype_p_{{group}}.csv", group=GROUPS),
        scp_boot=lambda wc: (
            [f"{OUT}/scpagwas/{wc.trait}/celltype_bootstrap_{CELLTYPE}.csv"]
            if PG["run_celltype_bootstrap"] else []
        ),
        scp_cell=f"{OUT}/scpagwas/{{trait}}/per_cell_scores.csv.gz",
        scdrs_grp=expand(f"{OUT}/scdrs/downstream/{{{{trait}}}}.scdrs_group.{{group}}", group=GROUPS),
        scdrs_full=f"{OUT}/scdrs/scores/{{trait}}.full_score.gz",
    output:
        tables=expand(f"{OUT}/summary/{{{{trait}}}}/summary_{{group}}.csv", group=GROUPS),
        corr=f"{OUT}/summary/{{trait}}/scPagwas_vs_scDRS_percell_corr.csv",
        figs=expand(f"{OUT}/summary/{{{{trait}}}}/summary_{{group}}.pdf", group=GROUPS),
    params:
        groups=GROUPS,
        celltype_col=CELLTYPE,
    log:
        f"{OUT}/logs/summary/{{trait}}.log",
    resources:
        mem_mb=res("summarize")["mem_gb"] * 1024,
        runtime=res("summarize")["time_min"],
    script:
        "scripts/10_summarize_trait.py"


localrules: all, summarize_all, integrate_all, ldsc_coords, ldsc_print_snps


rule summarize_all:
    input:
        expand(f"{OUT}/summary/{{trait}}/summary_{{group}}.csv", trait=TRAITS, group=GROUPS),
    output:
        f"{OUT}/summary/all_traits_celltype_summary.csv",
    run:
        import pandas as pd
        pd.concat([pd.read_csv(f) for f in input], ignore_index=True).to_csv(output[0], index=False)

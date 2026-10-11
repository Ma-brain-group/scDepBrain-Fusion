"""
scDRS step C3: MAGMA .genes.out -> scDRS pval file (GENE <tab> <trait>).
Entrez IDs are mapped to symbols with the MAGMA gene.loc file;
duplicated symbols keep the smallest P.
"""
import sys

import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

trait = snakemake.wildcards.trait
d = pd.read_csv(snakemake.input[0], sep=r"\s+", dtype={"GENE": str})
for c in ("GENE", "P"):
    if c not in d:
        raise ValueError(f"MAGMA file lacks column {c}; columns = {list(d.columns)}")
print(f"MAGMA genes: {len(d):,}")

mode = snakemake.params.gene_id
if mode == "auto":
    mode = "entrez" if d["GENE"].str.fullmatch(r"\d+").mean() > 0.9 else "symbol"
print(f"Gene ID type: {mode}")

if mode == "entrez":
    gloc = pd.read_csv(snakemake.params.gene_loc, sep=r"\s+", header=None, dtype=str,
                       names=["entrez", "chr", "start", "end", "strand", "symbol"])
    d["GENE"] = d["GENE"].map(dict(zip(gloc["entrez"], gloc["symbol"])))
    n_na = d["GENE"].isna().sum()
    print(f"Unmapped Entrez IDs dropped: {n_na:,}")
    d = d.dropna(subset=["GENE"])

out = (d[["GENE", "P"]].sort_values("P").drop_duplicates("GENE")
       .rename(columns={"P": trait}).set_index("GENE"))
out.to_csv(snakemake.output[0], sep="\t")
print(f"Wrote {len(out):,} genes; P < 2.5e-6: {(out[trait] < 2.5e-6).sum()}")

"""
Prepare GWAS summary statistics for LDSC munge_sumstats.py:
    SNP  A1  A2  N  BETA  P      (A1 = effect allele; tab-separated, gzipped)

format table : columns from gwas_columns / traits.<name>.columns
               (rsid, a1, a2, beta or or, se and/or p, n) ; N may be a constant (traits.<name>.n)
format vcf   : IEU OpenGWAS GWAS-VCF; A1 = ALT, A2 = REF, BETA = ES, SE = SE,
               P from LP (-log10 P) if present, else from BETA / SE; N from SS or traits.<name>.n
Only autosomal rsIDs are kept; munge_sumstats.py then restricts to HapMap3 SNPs.
"""
import gzip
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.stderr = sys.stdout = open(snakemake.log[0], "w")
src, out = snakemake.input[0], snakemake.output[0]
fmt = snakemake.params.fmt
USED = {"chrom", "rsid", "a1", "a2", "n", "beta", "or", "se", "p"}
cols = {k: v for k, v in snakemake.params.columns.items() if v and k in USED}
n_const = snakemake.params.n


def num(x):
    return pd.to_numeric(x, errors="coerce")


opener = gzip.open if src.endswith(".gz") else open
n_skip = 0
with opener(src, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            n_skip += 1
            continue
        header = line.rstrip("\n")
        break
sep = "\t" if "\t" in header else r"\s+"
hcols = header.split("\t") if sep == "\t" else header.split()

if fmt == "vcf":
    sample = hcols[-1]
    use = ["#CHROM", "ID", "REF", "ALT", "FORMAT", sample]
else:
    for req in ("chrom", "rsid", "a1", "a2"):
        if req not in cols:
            raise ValueError(f"LDSC needs gwas_columns.{req}")
    if "n" not in cols and n_const is None:
        raise ValueError("LDSC needs a sample size: gwas_columns.n or traits.<name>.n")
    use = list(cols.values())
    missing = [c for c in use if c not in hcols]
    if missing:
        raise ValueError(f"Columns not in GWAS header: {missing}")

parts = []
for chunk in pd.read_csv(src, sep=sep, skiprows=n_skip, usecols=use, dtype=str, chunksize=1_000_000, engine="c"):
    if fmt == "vcf":
        rows = []
        for fstr, sub in chunk.groupby("FORMAT", sort=False):
            keys = fstr.split(":")
            v = sub[sample].str.split(":", expand=True).iloc[:, :len(keys)]
            v.columns = keys[:v.shape[1]]
            beta, se = num(v["ES"]), num(v["SE"])
            p = 10 ** (-num(v["LP"])) if "LP" in v else 2 * norm.sf(np.abs(beta / se))
            n = num(v["SS"]) if "SS" in v else pd.Series(n_const, index=sub.index, dtype=float)
            rows.append(pd.DataFrame({"chrom": sub["#CHROM"], "SNP": sub["ID"], "A1": sub["ALT"], "A2": sub["REF"],
                                      "N": n, "BETA": beta, "P": p}, index=sub.index))
        d = pd.concat(rows)
    else:
        inv = {v: k for k, v in cols.items()}
        c = chunk.rename(columns=inv)
        beta = num(c["beta"]) if "beta" in c else np.log(num(c["or"]))
        if "p" in c:
            p = num(c["p"])
        else:
            p = 2 * norm.sf(np.abs(beta / num(c["se"])))
        n = num(c["n"]) if "n" in c else pd.Series(n_const, index=c.index, dtype=float)
        d = pd.DataFrame({"chrom": c["chrom"], "SNP": c["rsid"], "A1": c["a1"].str.upper(), "A2": c["a2"].str.upper(),
                          "N": n, "BETA": beta, "P": p})
    chrom = num(d["chrom"].astype(str).str.replace(r"^chr", "", case=False, regex=True))
    keep = chrom.between(1, 22) & d["SNP"].astype(str).str.startswith("rs") & d[["N", "BETA", "P"]].notna().all(axis=1)
    parts.append(d.loc[keep, ["SNP", "A1", "A2", "N", "BETA", "P"]])

res = pd.concat(parts, ignore_index=True)
if res.empty:
    raise ValueError("No SNPs left for LDSC")
res.to_csv(out, sep="\t", index=False, compression="gzip")
print(f"Wrote {len(res):,} SNPs; median N = {res['N'].median():,.0f}")

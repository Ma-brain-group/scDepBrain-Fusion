"""
Convert GWAS summary statistics into the 6-column scPagwas input format:
    chrom pos rsid se beta maf   (space-delimited)

Two input formats (snakemake.params.fmt):
  table : delimited text; column names given by snakemake.params.columns
            beta : BETA, or log(OR) if only OR is given
            se   : SE, or |beta| / z(p) if only P is given
            maf  : MAF -> EAF -> case/control-size-weighted FCAS/FCON
  vcf   : IEU OpenGWAS GWAS-VCF (FORMAT ES:SE:LP:AF:ID)
            beta = ES, se = SE, maf = min(AF, 1 - AF), rsid = ID column
Filters: autosomes, rsID present (rs*), optional excluded region, MAF > min_maf.
"""
import gzip
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

src = snakemake.input[0]
out = snakemake.output[0]
fmt = snakemake.params.fmt
filt = snakemake.params.filt
# only the columns this step uses (a1 / a2 / n are for LDSC)
USED = {"chrom", "pos", "rsid", "beta", "or", "se", "p", "maf", "eaf", "fcas", "fcon", "ncas", "ncon"}
cols = {k: v for k, v in snakemake.params.columns.items() if v and k in USED}

# ---- locate header (skip '##' meta lines) ----
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
header_cols = header.split("\t") if sep == "\t" else header.split()
print(f"Input: {src}\nFormat: {fmt}\nHeader ({len(header_cols)} cols): {header_cols}")


def num(x):
    return pd.to_numeric(x, errors="coerce")


# ---- per-format: produce a standard frame with chrom,pos,rsid,beta,se,f ----
if fmt == "vcf":
    sample_col = header_cols[-1]
    usecols = ["#CHROM", "POS", "ID", "FORMAT", sample_col]

    def standardise(chunk):
        # FORMAT can differ between rows (e.g. ES:SE:LP:AF:ID vs ES:SE:LP:AF),
        # so parse each FORMAT group with its own field order.
        parts = []
        for fmt_str, sub in chunk.groupby("FORMAT", sort=False):
            keys = fmt_str.split(":")
            for req in ("ES", "SE", "AF"):
                if req not in keys:
                    raise ValueError(f"VCF FORMAT lacks {req}: {keys}")
            vals = sub[sample_col].str.split(":", expand=True)
            vals = vals.iloc[:, :len(keys)]
            vals.columns = keys[:vals.shape[1]]
            parts.append(pd.DataFrame({
                "chrom": sub["#CHROM"], "pos": num(sub["POS"]), "rsid": sub["ID"],
                "beta": num(vals["ES"]), "se": num(vals["SE"]), "f": num(vals["AF"]),
            }, index=sub.index))
        return pd.concat(parts).loc[chunk.index]

elif fmt == "table":
    missing = [v for v in cols.values() if v not in header_cols]
    if missing:
        raise ValueError(f"Columns not found in GWAS header: {missing}. "
                         f"Fix gwas_columns / traits.<name>.columns in the config.")
    for req in ("chrom", "pos", "rsid"):
        if req not in cols:
            raise ValueError(f"gwas_columns.{req} is required")
    if "beta" not in cols and "or" not in cols:
        raise ValueError("Need gwas_columns.beta or gwas_columns.or")
    if "se" not in cols and "p" not in cols:
        raise ValueError("Need gwas_columns.se or gwas_columns.p")
    if not ("maf" in cols or "eaf" in cols
            or all(k in cols for k in ("fcas", "fcon", "ncas", "ncon"))):
        raise ValueError("Need an allele-frequency source: maf, eaf, or fcas+fcon+ncas+ncon")
    usecols = list(cols.values())
    inv = {v: k for k, v in cols.items()}

    def standardise(chunk):
        d = chunk.rename(columns=inv)
        beta = num(d["beta"]) if "beta" in d else np.log(num(d["or"]))
        if "se" in d:
            se = num(d["se"])
        else:
            p = num(d["p"]).clip(lower=1e-300)
            se = beta.abs() / norm.isf(p / 2)
        if "maf" in d:
            f = num(d["maf"])
        elif "eaf" in d:
            f = num(d["eaf"])
        else:
            fa, fc, na, nc = (num(d[k]) for k in ("fcas", "fcon", "ncas", "ncon"))
            f = (fa * na + fc * nc) / (na + nc)
        return pd.DataFrame({"chrom": d["chrom"], "pos": num(d["pos"]), "rsid": d["rsid"],
                             "beta": beta, "se": se, "f": f})
else:
    raise ValueError(f"Unknown GWAS format '{fmt}' (use 'table' or 'vcf')")

reader = pd.read_csv(src, sep=sep, skiprows=n_skip, usecols=usecols, dtype=str,
                     chunksize=1_000_000, engine="c")

n_in = 0
drop = {"no_rsid": 0, "non_autosome": 0, "excluded_region": 0, "bad_values": 0, "low_maf": 0}
stats = {}
first = True
reg = filt.get("exclude_region")
min_maf = float(filt.get("min_maf", 0))
with open(out, "w") as fo:
    for chunk in reader:
        n_in += len(chunk)
        d = standardise(chunk)
        chrom_n = num(d["chrom"].str.replace(r"^chr", "", case=False, regex=True))
        maf = np.minimum(d["f"], 1 - d["f"])

        ok_rsid = d["rsid"].notna() & d["rsid"].str.startswith("rs")
        ok_chr = chrom_n.between(1, 22) if filt.get("autosomes_only", True) else chrom_n.notna()
        in_reg = ((chrom_n == int(reg["chrom"])) & d["pos"].between(reg["start"], reg["end"])
                  if reg else pd.Series(False, index=d.index))
        ok_val = d["pos"].notna() & d["beta"].notna() & d["se"].notna() & (d["se"] > 0) & maf.notna()
        ok_maf = maf > min_maf

        # sequential filter accounting
        m = ok_rsid.copy();          drop["no_rsid"] += int((~ok_rsid).sum())
        drop["non_autosome"] += int((m & ~ok_chr).sum());   m &= ok_chr
        drop["excluded_region"] += int((m & in_reg).sum()); m &= ~in_reg
        drop["bad_values"] += int((m & ~ok_val).sum());     m &= ok_val
        drop["low_maf"] += int((m & ~ok_maf).sum());        m &= ok_maf

        res = pd.DataFrame({
            "chrom": chrom_n[m].astype(int),
            "pos": d["pos"][m].astype(np.int64),
            "rsid": d["rsid"][m],
            "se": d["se"][m],
            "beta": d["beta"][m],
            "maf": maf[m],
        })
        res.to_csv(fo, sep=" ", index=False, header=first, float_format="%.6g")
        first = False
        for c, n in res["chrom"].value_counts().items():
            stats[c] = stats.get(c, 0) + n

n_out = sum(stats.values())
print(f"Rows in: {n_in:,}   rows out: {n_out:,}")
print("Dropped (sequential):", {k: f"{v:,}" for k, v in drop.items()})
print("Per-chromosome:", dict(sorted(stats.items())))
if n_out == 0:
    raise ValueError("No SNPs survived filtering — check format / column mapping / filters")

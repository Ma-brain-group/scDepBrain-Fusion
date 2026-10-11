"""
Helper for 12_ldsc_l2_chr.sh (runs in the LDSC Python 2 environment; also py3-compatible).

merge <tmp> <combined.annot.gz> <names...>
    Combine single-column thin annotations (<tmp>/<name>.annot.gz, column ANNOT)
    into one thin annotation with one column per gene set.
split <tmp> <out> <chr> <names...>
    Split <tmp>/combined.<chr>.l2.ldscore.gz / .l2.M / .l2.M_5_50 into
    <out>/<name>/<name>.<chr>.{annot.gz, l2.ldscore.gz, l2.M, l2.M_5_50}.
"""
from __future__ import print_function
import os
import shutil
import sys

import pandas as pd


def merge(tmp, out, names):
    cols = {}
    for n in names:
        d = pd.read_csv(os.path.join(tmp, n + ".annot.gz"), sep="\t", compression="gzip")
        cols[n] = d["ANNOT"].values
    df = pd.DataFrame(cols, columns=names)
    df.to_csv(out, sep="\t", index=False, compression="gzip")
    print("merged", len(names), "annotations,", len(df), "SNPs")


def split(tmp, out, chrom, names):
    pre = os.path.join(tmp, "combined.%s" % chrom)
    ld = pd.read_csv(pre + ".l2.ldscore.gz", sep="\t", compression="gzip")
    M = open(pre + ".l2.M").read().split()
    M550 = open(pre + ".l2.M_5_50").read().split()
    if len(M) != len(names) or len(M550) != len(names):
        raise ValueError("M files do not match the number of annotations")
    for i, n in enumerate(names):
        d = os.path.join(out, n)
        if not os.path.isdir(d):
            os.makedirs(d)
        p = os.path.join(d, "%s.%s" % (n, chrom))
        col = n + "L2"
        if col not in ld.columns:
            raise ValueError("column %s missing in combined LD scores" % col)
        ld[["CHR", "SNP", "BP", col]].to_csv(p + ".l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        with open(p + ".l2.M", "w") as f:
            f.write(M[i] + "\n")
        with open(p + ".l2.M_5_50", "w") as f:
            f.write(M550[i] + "\n")
        # single-column annotation, renamed so names stay unique across --ref-ld-chr prefixes
        a = pd.read_csv(os.path.join(tmp, n + ".annot.gz"), sep="\t", compression="gzip")
        a.columns = [n]
        a.to_csv(p + ".annot.gz", sep="\t", index=False, compression="gzip")
    print("split", len(names), "annotations for chr", chrom)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "merge":
        merge(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif mode == "split":
        split(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
    else:
        sys.exit("mode must be merge or split")

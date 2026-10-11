"""Collect S-LDSC results for one trait: the last annotation row of each
<set>.results is the cell-type annotation (after baseline and all_gene).
Coefficient P is one-sided (H1: coefficient > 0), as in LDSC-SEG."""
import os
import sys

import pandas as pd
from scipy.stats import norm

sys.stderr = sys.stdout = open(snakemake.log[0], "w")
names = pd.read_csv(snakemake.input.names)
rows = []
for _, r in names.iterrows():
    f = os.path.join(snakemake.params.res_dir, r["set"] + ".results")
    if not os.path.exists(f):
        raise FileNotFoundError(f)
    d = pd.read_csv(f, sep="\t").iloc[-1]
    z = float(d["Coefficient_z-score"])
    rows.append({"group": r["group"], "celltype": r["celltype"], "Coefficient": d["Coefficient"],
                 "Coefficient_se": d["Coefficient_std_error"], "Coefficient_z": z, "LDSC_P": norm.sf(z),
                 "Enrichment": d["Enrichment"], "Enrichment_p": d["Enrichment_p"], "Prop_SNPs": d["Prop._SNPs"]})
out = pd.DataFrame(rows)
out.to_csv(snakemake.output[0], index=False)
print(out.to_string(index=False))

"""
Step D: integrate scPagwas and scDRS cell-type results for one trait.

Per annotation column (group):
  summary_<group>.csv : one row per cell type with
      scPagwas Merge_celltype_p (+FDR), scPagwas bootstrap (celltype_col only),
      scDRS assoc_mcp / hetero_mcp (+FDR), n_fdr_0.1 cells
  summary_<group>.pdf : -log10 P heatmap (methods x cell types), FDR<0.05 marked
Also: per-cell Spearman correlation of scPagwas TRS vs scDRS norm_score
(overall and within each celltype_col group).
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

trait = snakemake.wildcards.trait
groups = list(snakemake.params.groups)
ct_col = snakemake.params.celltype_col


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


boot = None
if snakemake.input.scp_boot:
    boot = pd.read_csv(snakemake.input.scp_boot[0])
    boot = boot[boot["annotation"] != "Intercept"].rename(columns={"annotation": "celltype"})
    boot = boot[["celltype", "bt_value", "bp_value"]].rename(
        columns={"bt_value": "scPagwas_boot_t", "bp_value": "scPagwas_boot_p"})
    boot["scPagwas_boot_fdr"] = bh(boot["scPagwas_boot_p"])

for g, f_scp, f_drs, f_tab, f_fig in zip(groups, snakemake.input.scp_ct, snakemake.input.scdrs_grp,
                                          snakemake.output.tables, snakemake.output.figs):
    scp = pd.read_csv(f_scp).rename(columns={"pvalue": "scPagwas_p", "fdr": "scPagwas_fdr"})
    drs = pd.read_csv(f_drs, sep="\t").rename(columns={"group": "celltype"})
    drs = drs[["celltype", "n_cell", "assoc_mcp", "assoc_mcz", "hetero_mcp", "n_fdr_0.1"]].rename(
        columns={"n_cell": "scDRS_n_cell", "assoc_mcp": "scDRS_assoc_p", "assoc_mcz": "scDRS_assoc_z",
                 "hetero_mcp": "scDRS_hetero_p", "n_fdr_0.1": "scDRS_n_cell_fdr0.1"})
    drs["scDRS_assoc_fdr"] = bh(drs["scDRS_assoc_p"])
    drs["scDRS_hetero_fdr"] = bh(drs["scDRS_hetero_p"])

    tab = scp.merge(drs, on="celltype", how="outer")
    if boot is not None and g == ct_col:
        tab = tab.merge(boot, on="celltype", how="left")
    tab.insert(0, "group_col", g)
    tab.insert(0, "trait", trait)
    tab = tab.sort_values("scDRS_assoc_p", na_position="last")
    tab.to_csv(f_tab, index=False)
    print(f"\n=== {trait} / {g} ===\n{tab.to_string(index=False)}")

    # ---- heatmap: methods (rows) x cell types (cols), sequential single hue ----
    methods = [("scPagwas_p", "scPagwas_fdr", "scPagwas\nMerge_celltype_p"),
               ("scDRS_assoc_p", "scDRS_assoc_fdr", "scDRS\nassociation"),
               ("scDRS_hetero_p", "scDRS_hetero_fdr", "scDRS\nheterogeneity")]
    if "scPagwas_boot_p" in tab:
        methods.insert(1, ("scPagwas_boot_p", "scPagwas_boot_fdr", "scPagwas\nbootstrap"))
    cts = tab["celltype"].astype(str).tolist()
    M = np.vstack([-np.log10(tab[p].astype(float).clip(lower=1e-300)) for p, _, _ in methods])
    F = np.vstack([tab[q].astype(float) < 0.05 for _, q, _ in methods])

    fig, ax = plt.subplots(figsize=(max(5, 0.45 * len(cts) + 2.5), 0.6 * len(methods) + 2.6))
    im = ax.imshow(M, cmap="Blues", aspect="auto", vmin=0, vmax=max(2.0, np.nanmax(M)))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if F[i, j]:
                ax.text(j, i, "*", ha="center", va="center", fontsize=13,
                        color="white" if M[i, j] > 0.6 * im.get_clim()[1] else "#1a1a1a")
    ax.set_xticks(range(len(cts)), cts, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(methods)), [m[2] for m in methods], fontsize=8)
    ax.set_xticks(np.arange(-.5, len(cts), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(methods), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("-log10 P", fontsize=8)
    cb.outline.set_visible(False)
    ax.set_title(f"{trait} — {g}   (* BH FDR < 0.05)", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(f_fig)
    fig.savefig(f_fig.replace(".pdf", ".png"), dpi=200)
    plt.close(fig)

# ---- per-cell scPagwas TRS vs scDRS norm_score ----
cell = pd.read_csv(snakemake.input.scp_cell, usecols=["cell_id", ct_col, "scPagwas.TRS.Score"])
# first column = cell barcode (named "index" or "cell_id" depending on the h5ad)
drs = pd.read_csv(snakemake.input.scdrs_full, sep="\t", index_col=0, usecols=lambda c: True)
drs = drs[["norm_score"]].rename_axis("cell_id").reset_index()
m = cell.merge(drs, on="cell_id")
print(f"\nCells matched for correlation: {len(m):,}")
rows = []
for name, sub in [("ALL", m)] + list(m.groupby(ct_col, observed=True)):
    if len(sub) >= 10:
        r, p = spearmanr(sub["scPagwas.TRS.Score"], sub["norm_score"])
        rows.append({"trait": trait, "group": ct_col, "celltype": name,
                     "n_cell": len(sub), "spearman_rho": r, "p": p})
pd.DataFrame(rows).to_csv(snakemake.output.corr, index=False)
print(pd.DataFrame(rows).to_string(index=False))

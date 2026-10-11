"""
scDRS step C2: assemble an AnnData (X = raw counts, cells x genes) from
the CSC arrays written by 07_export_counts.R, plus the covariate file.
"""
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import os

sys.stderr = sys.stdout = open(snakemake.log[0], "w")


def read_str(ds):
    try:
        return ds.asstr()[:]
    except (AttributeError, TypeError):
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in ds[:]])


with h5py.File(snakemake.input.counts, "r") as f:
    genes = read_str(f["genes"])
    cells = read_str(f["cells"])
    X = sp.csr_matrix(
        (f["x"][:].astype(np.float32), f["i"][:].astype(np.int64), f["p"][:].astype(np.int64)),
        shape=(len(cells), len(genes)),
    )
print(f"X: {X.shape[0]:,} cells x {X.shape[1]:,} genes, nnz = {X.nnz:,}")

obs = pd.read_csv(snakemake.input.meta, index_col="cell_id", low_memory=False)
obs.index = obs.index.astype(str)
obs = obs.loc[cells]
for g in snakemake.params.group_cols:
    obs[g] = obs[g].astype(str).astype("category")
    print(f"{g}: {obs[g].nunique()} groups")

adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
adata.var_names_make_unique()

# ---- neighbour graph (scDRS heterogeneity / Geary's C) ----
# Only when a Seurat reduction is given. Otherwise no graph is stored and
# `scdrs perform-downstream` builds its default one (PCA 20 PCs, k = 15 on the
# scDRS-normalised data), which is what the original scDRS analysis used.
if os.path.getsize(snakemake.input.embedding) > 0:
    k = int(snakemake.params.n_neighbors)
    emb = pd.read_csv(snakemake.input.embedding, index_col="cell_id")
    emb.index = emb.index.astype(str)
    adata.obsm["X_emb"] = emb.loc[cells].to_numpy(dtype=np.float32)
    print(f"Neighbours on provided embedding: {adata.obsm['X_emb'].shape[1]} dims, k = {k}")
    sc.pp.neighbors(adata, use_rep="X_emb", n_neighbors=k)
else:
    print("No embedding given: scdrs perform-downstream will build its default KNN graph")
print("obsp keys:", list(adata.obsp.keys()))
adata.write_h5ad(snakemake.output.h5ad)
print(f"Wrote {snakemake.output.h5ad}")

# Covariates for scDRS: intercept + numeric columns + one-hot categoricals
cov = pd.DataFrame({"const": 1.0}, index=obs.index)
for c in snakemake.params.cov_cols:
    if c not in obs:
        raise ValueError(f"cov column '{c}' not in metadata")
    v = obs[c]
    if pd.api.types.is_numeric_dtype(v):
        cov[c] = v.fillna(v.mean()).astype(float)
    else:
        cov = cov.join(pd.get_dummies(v.astype(str), prefix=c, drop_first=True, dtype=float))
cov.index.name = "index"
cov.to_csv(snakemake.output.cov, sep="\t")
print(f"Covariates: {list(cov.columns)}")

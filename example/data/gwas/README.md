# Demo GWAS (download required)

Place the PGC MDD 2025 European-ancestry summary statistics (excluding 23andMe) here:

    pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz

Source: Psychiatric Genomics Consortium results page,
https://pgc.unc.edu/for-researchers/download-results/  (MDD; Adams et al., Cell 2025).
The file is not redistributed in this repository because of the PGC data-use terms.

Format: tab-separated, GRCh37, columns
`#CHROM POS ID EA NEA BETA SE PVAL FCAS FCON IMPINFO NEFF NCAS NCON HETI HETDF HETPVAL`
(the default `gwas_columns` mapping in the configs).

A symbolic link works too:

    ln -s /path/to/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz example/data/gwas/

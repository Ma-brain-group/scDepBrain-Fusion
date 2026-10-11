# Install R packages that are not on conda: SOAR (archived on CRAN, required
# by scPagwas) and scPagwas itself (GitHub). Run inside the activated env:
#   Rscript envs/install_R_packages.R
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

if (!requireNamespace("SOAR", quietly = TRUE)) {
  remotes::install_url("https://cran.r-project.org/src/contrib/Archive/SOAR/SOAR_0.99-11.tar.gz",
                       upgrade = "never")
}
if (!requireNamespace("scPagwas", quietly = TRUE)) {
  # Developed / tested with scPagwas 2.0.0
  remotes::install_github("sulab-wmu/scPagwas", upgrade = "never")
}

# Sanity check
for (p in c("scPagwas", "Seurat", "SeuratObject", "hdf5r", "bigreadr", "data.table", "Matrix")) {
  cat(sprintf("%-14s %s\n", p, as.character(packageVersion(p))))
}
suppressMessages(library(scPagwas))
utils::data("block_annotation_hg37", package = "scPagwas")
cat("scPagwas data OK:", exists("block_annotation_hg37"), "\n")

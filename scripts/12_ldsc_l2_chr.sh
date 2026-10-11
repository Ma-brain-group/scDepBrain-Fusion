#!/bin/bash
# =====================================================================
# LDSC annotations + LD scores for ALL gene sets on one chromosome.
#   1. make_annot.py per gene set (gene +/- window -> SNPs of the 1000G bim)
#   2. merge into one multi-column thin annotation and run `ldsc.py --l2` ONCE
#      (LD scores are computed independently per annotation column)
#   3. split back into one annot / ldscore / M / M_5_50 file set per gene set,
#      the layout S-LDSC (--ref-ld-chr, --overlap-annot) expects
# Usage: 12_ldsc_l2_chr.sh <chr> <genesets_dir> <coords> <plink_prefix> <print_snps>
#                          <window_bp> <ldsc_dir> <out_dir> <scripts_dir>
# =====================================================================
set -euo pipefail
CHR=$1; GS_DIR=$2; COORD=$3; PLINK=$4; PRINT=$5; WIN=$6; LDSC=$7; OUT=$8; SCRIPTS=$9
TMP=$OUT/tmp_chr$CHR; mkdir -p $TMP
names=()
for f in $(ls $GS_DIR/*.geneset.txt | sort); do
  n=$(basename $f .geneset.txt); names+=($n)
  python $LDSC/make_annot.py --gene-set-file $f --gene-coord-file $COORD --windowsize $WIN \
         --bimfile ${PLINK}${CHR}.bim --annot-file $TMP/$n.annot.gz
done
echo "[$(date)] chr $CHR: ${#names[@]} annotations made"
python $SCRIPTS/ldsc_split.py merge $TMP $TMP/combined.$CHR.annot.gz "${names[@]}"
python $LDSC/ldsc.py --l2 --bfile ${PLINK}${CHR} --ld-wind-cm 1 --annot $TMP/combined.$CHR.annot.gz \
       --thin-annot --print-snps $PRINT --out $TMP/combined.$CHR
python $SCRIPTS/ldsc_split.py split $TMP $OUT $CHR "${names[@]}"
rm -rf $TMP
echo "[$(date)] chr $CHR done"

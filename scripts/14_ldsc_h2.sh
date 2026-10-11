#!/bin/bash
# =====================================================================
# S-LDSC (LDSC-SEG): one regression per cell-type gene set,
#   --ref-ld-chr baselineLD, all_gene (control), cell type   --overlap-annot
# Usage: 14_ldsc_h2.sh <sumstats.gz> <ldscore_dir> <baseline_prefix> <frq_prefix>
#                      <weights_prefix> <ldsc_dir> <out_dir> <threads>
# =====================================================================
set -euo pipefail
SS=$1; LD=$2; BASE=$3; FRQ=$4; W=$5; LDSC=$6; OUT=$7; NT=$8
mkdir -p $OUT
sets=$(ls -d $LD/*/ | xargs -n1 basename | grep -v -x -e all_gene -e tmp_chr.* | sort)
echo "$sets" | xargs -P $NT -I{} bash -c '
  s={}; o='"$OUT"'/$s
  [[ -s $o.results ]] && exit 0
  python '"$LDSC"'/ldsc.py --h2 '"$SS"' \
    --ref-ld-chr '"$BASE"','"$LD"'/all_gene/all_gene.,'"$LD"'/$s/$s. \
    --frqfile-chr '"$FRQ"' --w-ld-chr '"$W"' \
    --overlap-annot --print-coefficients --out $o > /dev/null
  echo "done $s"'
echo "[$(date)] $(ls $OUT/*.results | wc -l) regressions finished"

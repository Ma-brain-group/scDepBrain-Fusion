#!/bin/bash
# =====================================================================
# Launch the scDepBrain-Fusion Snakemake pipeline.
#
#   bash run_pipeline.sh <config.yaml> <local|slurm> [extra snakemake args]
#
# Examples
#   bash run_pipeline.sh config/config_example.yaml local -n    # dry run
#   bash run_pipeline.sh config/config_example.yaml local       # run here
#   bash run_pipeline.sh config/my_config.yaml slurm            # SLURM
#
# local : runs Snakemake in the foreground (profiles/local).
# slurm : submits the Snakemake controller as a small, long SLURM job; the
#         controller submits every rule as its own job (profiles/slurm).
#         Set SBATCH_PARTITION / SBATCH_ACCOUNT first if your cluster needs
#         them. Controller options can be changed with CONTROLLER_SBATCH_OPTS.
# The software environment must be active (conda activate scdepbrain-fusion), or
# set `env_activate` in the config.
# =====================================================================
set -euo pipefail
if [[ $# -lt 2 ]]; then
  sed -n '3,20p' "$0"; exit 1
fi
CFG="$1"; MODE="$2"; shift 2
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPE"
[[ -f "$CFG" ]] || { echo "Config not found: $CFG"; exit 1; }
command -v snakemake >/dev/null || { echo "snakemake not found: activate the conda env first"; exit 1; }

DRY=false
for a in "$@"; do [[ "$a" == "-n" || "$a" == "--dry-run" || "$a" == "--dryrun" ]] && DRY=true; done

case "$MODE" in
  local)
    exec snakemake --profile profiles/local --configfile "$CFG" "$@"
    ;;
  slurm)
    if $DRY; then
      exec snakemake --profile profiles/slurm --configfile "$CFG" "$@"
    fi
    mkdir -p logs
    OPTS="${CONTROLLER_SBATCH_OPTS:--c 1 --mem=8G -t 14-00:00:00}"
    # shellcheck disable=SC2086
    sbatch -J scDepBrain_smk $OPTS -o logs/controller_%j.out -e logs/controller_%j.err \
      --wrap "cd '$PIPE' && snakemake --profile profiles/slurm --configfile '$CFG' $*"
    echo "Controller submitted. Follow progress with: tail -f $PIPE/logs/controller_*.err"
    ;;
  *)
    echo "Mode must be 'local' or 'slurm'"; exit 1 ;;
esac

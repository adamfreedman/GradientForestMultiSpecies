#!/bin/bash
#SBATCH --job-name=gfms_snakemake
#SBATCH --partition=sapphire
#SBATCH --time=3-00:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/snakemake_head_%j.log
#
# Snakemake controller job: submits every rule as its own SLURM job via
# profiles/cannon. Run from the repository root:
#   sbatch run_slurm.sh                 # whole workflow
#   sbatch run_slurm.sh -n              # extra arguments are passed to snakemake
#   sbatch run_slurm.sh --use-all-predictors
#       no predictor selection: every candidate variable goes into the model;
#       results go to <outdir>_all_predictors (e.g. results_all_predictors)
#   sbatch run_slurm.sh --importance-floor 0.01
#       predictor selection as usual, then drop kept variables whose share of the
#       summed screening importance is below 1%; results go to
#       <outdir>_importance_floor0.01. Default: no floor.
#
# Rules that need internet (install_gf, get_boundaries) are local rules; run
# them once on a login node first:
#   snakemake --use-conda --cores 1 results/.gf_installed results/env/boundaries.geojson
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate my_snakemake
args=()
extra_config=()
want_floor=0
for a in "$@"; do
  if [ "$want_floor" = 1 ]; then
    extra_config+=("importance_floor=$a"); want_floor=0; continue
  fi
  case "$a" in
    --use-all-predictors) extra_config+=(use_all_predictors=true) ;;
    --importance-floor=*) extra_config+=("importance_floor=${a#*=}") ;;
    --importance-floor) want_floor=1 ;;
    *) args+=("$a") ;;
  esac
done
if [ "$want_floor" = 1 ]; then
  echo "--importance-floor needs a value, e.g. --importance-floor 0.01" >&2; exit 1
fi
if [ ${#extra_config[@]} -gt 0 ]; then extra_config=(--config "${extra_config[@]}"); fi
# --config goes last so it cannot swallow targets (it takes a list of values)
snakemake --profile profiles/cannon "${args[@]}" "${extra_config[@]}"

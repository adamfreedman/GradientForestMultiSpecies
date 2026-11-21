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
#
# Rules that need internet (install_gf, get_boundaries) are local rules; run
# them once on a login node first:
#   snakemake --use-conda --cores 1 results/.gf_installed results/env/boundaries.geojson
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate my_snakemake
snakemake --profile profiles/cannon "$@"

#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=youremail@university.edu
#SBATCH -N 1
#SBATCH --account=yourlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G # modify if you need less or more
#SBATCH -J neutrality_search

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake \
  --profile slurm_general \
  --use-conda \
  --conda-prefix /hpc/group/wraylab/ab620/conda \
  --keep-going

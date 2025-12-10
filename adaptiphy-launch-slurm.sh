#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=youremail@university.edu
#SBATCH -N 1
#SBATCH --account=youraccount
#SBATCH --partition=common
#SBATCH --mem=12G
#SBATCH -J adaptiphy
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake \
  --profile slurm_general \
  --use-conda \
  --conda-prefix /hpc/path/your/conda \
  --keep-going

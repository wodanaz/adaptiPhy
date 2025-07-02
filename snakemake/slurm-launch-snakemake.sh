#!/bin/sh
#SBATCH --mail-user=apm58@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm.test.%j.out
#SBATCH -n 23
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=scavenger

module load Anaconda3
#source /hpc/home/apm58/.zshrc
source activate /hpc/group/wraylab/apm58/miniconda3/envs/snakemake

snakemake --workflow-profile slurm_general

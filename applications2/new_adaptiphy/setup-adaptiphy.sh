#!/usr/bin/env bash
# Performs setup and staging required by run-adaptiphy.sh
#
set -e

while getopts "g:" OPTION; do
    case $OPTION in
    g)
        export GENOME=$OPTARG
        ;;
    esac
done

if [ -z "$GENOME" ]
then
   echo "ERROR: You must specify a genome via the '-g <GENOME.maf> argument."
   echo "Example: $0 -g chromosome.maf"
   exit 1
fi

# create logs directory if necessary
mkdir -p logs


echo "Setup for adaptiphy-pipeline - Starting"

echo "Activating 'adaptiphy' conda environment"

if [ "$USE_MODULES" == "Y" ]
then
  # set default value for ANACONDAMODULE to "Anaconda3/2019.10-gcb02"
   ANACONDAMODULE="${ANACONDAMODULE-Anaconda3/2019.10-gcb02}"
   module load $ANACONDAMODULE
else
  # make sure conda is setup within this script
  source $CONDA_PREFIX/etc/profile.d/conda.sh
fi

conda activate adaptiphy


echo ""

echo "Setup part 1 - Extract query alignments"
sbatch --wait scripts/extract-queries-genome.sh
echo "Setup part 1 - Done"
echo ""


echo "Setup part 2 - Create a dictionary file for using picard tools"
sbatch --wait scripts/create-references.sh
echo "Setup part 2 - Done"
echo ""


echo "Setup for adaptiphy-pipeline - Done"

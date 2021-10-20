#!/usr/bin/env bash
set -e

ShowHelp()
{
   # Display Help
   echo "Runs a Slurm pipeline determining branch-specific positive selection from a BED file of peak regions, optionally staging data in and out."
   echo
   echo "usage: $0 -g genome -d datadir -i inputproject -t tree -b branch [-w workdir] [-j numjobs] [-k] [-x]"
   echo "options:"
   echo "-g genome        *.maf genome to use - required"
   echo "-d datadir       directory used to hold input and output files - required"
   echo "-i inputproject  project name  - required"
   echo "-f bedfile       bed file  - required"
   echo "-t tree          Phylogenetic tree in Newick format, example: '(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19)))' - required"
   echo "-b branch        branch to test, example: hg19  - required"
   echo "-w workdir       directory that will hold a tempdir - defaults to current directory"
   echo "-j numjobs       number of array jobs to run in parallel - defaults to 10"
   echo "-e email         email address to notify on pipeline completion - defaults to empty(no email sent)"
   echo "-k               somatic chromosome number"
   exho "-x               sex chromosomes?"
   echo ""
   echo "NOTE: The input genome must first be indexed by running ./setup-variants-pipeline.sh."
   echo "NOTE: The genome and datadir must be shared across the slurm cluster."
   echo ""
}


###################################
### Read command line arguments ###
###################################


# set default argument values
export WORKDIR=$(pwd)
export OUTDIR=$(pwd)
export LOGSUFFIX=$$
export DELETE_EVTMPDIR=Y
export PROJECTNAME=
export EVMODE=""
export DATETAB=""
export NEUTRALBED="neutral.bed"
export MAX_ARRAY_JOBS=10
export DELETE_EVTMPDIR=Y
export DOWNLOAD_INPUT_DATA=Y
export UPLOAD_OUTPUT_DATA=Y



while getopts "g:d:i:m:w:j:D:sSke:" OPTION; do
    case $OPTION in
    g)
        export GENOME=$OPTARG
        ;;
    d)
        export DATADIR=$(readlink -e $OPTARG)
        ;;
    m)
        export EVMODE=$OPTARG
        ;;
    i)
        export PROJECTNAME=$OPTARG
        ;;
    j)
        export MAX_ARRAY_JOBS=$OPTARG
        ;;
    w)
        export WORKDIR=$(readlink -e $OPTARG)
        ;;
    D)
        export DATETAB=$OPTARG
        ;;
    e)
        EMAIL=$OPTARG
        ;;
    k)
        export DELETE_EVTMPDIR=N
        ;;
    s)
        export DOWNLOAD_INPUT_DATA=N
        ;;
    S)
        export UPLOAD_OUTPUT_DATA=N
        ;;
    esac
done




####################################
### Check command line arguments ###
####################################

# declare directory names
INPUT_DATADIR=$DATADIR/input
export INPUTDIR=$INPUT_DATADIR/$PROJECTNAME
OUTPUT_DATADIR=$DATADIR/output
export OUTDIR=$OUTPUT_DATADIR/$PROJECTNAME
export LOGDIR="$OUTDIR/logs"


# check required arguments
if [ -z "$GENOME" ]
then
   echo "ERROR: Missing required '-g genome' argument."
   echo ""
   ShowHelp
   exit 1
fi
if [ -z "$DATADIR" ]
then
   echo "ERROR: Missing required '-d datadir' argument."
   echo ""
   ShowHelp
   exit 1
fi
if [ -z "$PROJECTNAME" ]
then
   echo "ERROR: Missing required '-i inputproject' argument."
   echo ""
   ShowHelp
   exit 1
fi

if [[ "$EVMODE" != "duke" && "$EVMODE" != "nc_state" ]]
then
   echo "ERROR: Required '-m mode' argument must be 'duke' or 'nc_state'."
   echo ""
   ShowHelp
   exit 1
fi

if [ -z "$DATETAB" ]
then
   echo "ERROR: Missing required '-D datetab' argument."
   echo ""
   ShowHelp
   exit 1
fi

# check that the genome dictionary has been created by ./setup-escape-variants.sh
GENOME_BASE_NAME=$(basename $GENOME .fasta)
GENOME_DICTIONARY="$GENOME_BASE_NAME.dict"
if [[ ! -f "$GENOME_DICTIONARY" ]]
then
    echo "ERROR: Genome dictionary file $GENOME_DICTIONARY not found."
    echo "To fix run ./setup-escape-variants.sh -g $GENOME"
    echo ""
    exit 1
fi

# check that there are no duplicate sample names in the input fastq files
BADSAMPLES=$(basename $INPUTDIR/*.fastq.gz | sed -e 's/_.*//' | uniq -d)
if [ "$BADSAMPLES" ]
then
   echo "Warning: Found duplicate sample names!"
   DUPDIR="${INPUTDIR}_dups"
   mkdir -p $DUPDIR
   for BADSAMPLE in $BADSAMPLES
   do
       echo "Moving $BADSAMPLE to $DUPDIR."
       mv $INPUTDIR/${BADSAMPLE}*.fastq.gz $DUPDIR/.
   done
fi

# check that the config.sh has been setup
if [ -f config.sh ]
then
   source config.sh
else
   echo "ERROR: Missing config.sh config file."
   echo "To fix run:"
   echo "  cp example-config.sh config.sh"
   echo ""
   exit 1
fi

# enable USE_MODULES if not explicitly turned off by config.sh
if [ "$USE_MODULES" != "N" ]
then
   export USE_MODULES=Y
fi

##########################
### Create directories ###
##########################

# Create directories
# make base input directory if necessary
mkdir -p $INPUT_DATADIR
# make base output directory if necessary
mkdir -p $OUTPUT_DATADIR
# make output results directory if necessary
mkdir -p $OUTDIR
# make output logs directory if necessary
mkdir -p $LOGDIR

###########################
### Run Escape Variants ###
###########################

SBATCH_FLAGS="$SBATCH_FLAGS --output=$LOGDIR/run-staging-pipeline-%j.out"
if [ ! -z "$EMAIL" ]
then
   echo "Emailing $EMAIL on pipeline completion."
   SBATCH_FLAGS="$SBATCH_FLAGS --mail-type=END --mail-user=$EMAIL"
fi

JOBID=$(sbatch --parsable ${SBATCH_FLAGS} scripts/run-staging-pipeline.sh)
echo "Submitted batch job $JOBID"
echo "To monitor main log run:"
echo "tail -f $LOGDIR/run-staging-pipeline-$JOBID.out"

# AdaptiPhy: Implementation with snakemake

This code updates the existing AdaptiPhy 2.0 pipeline to run using snakemake, allowing the user to plug in data at the beginning of the pipeline & performing the intervening steps automatically. We have deprecated the original version and has been moved into a new directory named 'deprecated'.

## Installing the AdaptiPhy pipeline to run with snakemake ##

### 1a. Recommended: Download the tarball release from the "Releases" page ###

add information here when this option is ready

### 1b. Not recommended: Clone the git repo ###

To clone this repo from the command line into your working directory, use:

```bash
git clone https://github.com/wodanaz/adaptiPhy
```

You will need to add your files to the ```data/``` directory before running the snakemake pipeline for the first time. Read on for more info about the necessary file structure in this folder!

### 2. Installing snakemake as a conda environment ###
The majority of the conda packages required in the pipeline will be loaded automatically, and you will not need to do any sort of manual install. However, you will only need to create an environment that contains snakemake and python in order to run this pipeline.

Please install it in your directory specified for conda environments:

```
nano snakemake.yml 
```

```
name: snakemake
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - snakemake=9.1
  - python=3.11
  - snakemake-executor-plugin-slurm
  - bedtools
  - bedops
  - bzip2
```


If you have conda in your system, please install this package:

```
conda env create --file snakemake.yml 
```

To invoke or activate do:

```
conda activate snakemake
```

Advanced: If you don't set a `--prefix`, make sure that your `.condarc` file has a specified location to save environments to. If you are using a shared computing space, not specifying a stable install location can result in a loss of environment files.

### 3. Change the necessary file parameters ###

To run the snakemake pipeline either interactively or through a job manager like SLURM, you will need to update some file paths and other information in the cloned repository.

1. ```./config.yaml```: This is the file where you will need to update the most information. You will need:
   * ```windows```: path to a list of ATAC peaks or similar targets in BED file format
   * ```tree_topology``` and ```foreground_branches```: provide a phylogenetic tree in Newick format and specify which branches are focal. For the first parameter, use standard Newick format. For the second, provide a vector of the focal branch names from the Newick tree.
   * ```maf_pattern``` and ```fa_pattern```: provide paths to files (wildcards permitted) for one or many chromosomes. Note that both the .maf (multi alignment) and .fa (nucleotide) files are required.
   * ```neutral_set```: provide a path to a .txt file that contains paths to neutral proxy files, or set this parameter to "goodalignments.txt" if running AdaptiPhy in local mode. More on this later!
   * ```chromosomes```: provide a vector of chromosomes to examine.
   
     Example:
   
```bash
# INPUT SPLITS ##############################################################################################################
windows: "data/thurman.bed"
num_replicates: 10
min_frac: 0.9

# TREE TOPOLOGY #############################################################################################################
tree_topology: "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))"
foreground_branches: ["hg19"]

# GENOME TARGET FILES #######################################################################################################
# provide the input file to be split by phast's msa_split here. this file can be in a .fasta, phylip, mpm, maf, or ss file
# format. msa_split will try to guess the contents.
# if this fails, the snakefile may need to be modified to have an --in-format parameter specifying the file type. We typically
# provide a MAF file.
maf_pattern: "data/{chrom}.primate.maf"
#if providing a MAF file, provide the reference sequence location here.
fa_pattern: "data/{chrom}.fa"

# LOCAL VS GLOBAL RUN SPECIFICATION ##########################################################################################
# If running a local version of adaptiphy, no neutral sequence is required. Set the parameter below to "goodalignments.txt".
# If running a global version of adaptiphy, provide a neutral set file. Keep in mind that if you perform a local run of
# adaptiphy (meaning that you set neutral_set to "goodalignments.txt") in a global (whole-genome) run, your neutral set
# is a random sampling of the genome, which may not have a significant effect (see Berrio et. al. BMC) but caveat emptor.
neutral_set: "neutral_smk/neutralset.txt"
#options are: local = "goodalignments.txt", global = path to neutral set
chromosomes: ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", "chrX"]
#"chr" if one sequence (i.e. viral genome, one chromosome only in the file provided") or specific chromosomes to target if
# using a multi-chromosome genome (i.e. "chr19", etc)
```
    
 3. ```data/```: your input data lives in this folder. To run the AdaptiPhy pipeline, this folder must contain:
    * a folder containing MAF and .fa files, matching the specified 'pattern' paths in your ```config.yaml``` file from the previous step
    * a file of target windows/peak calls, matching the 'windows' path in  your ```config.yaml``` file from the previous step
    * if running AdaptiPhy in global mode (more on this later), a .txt file containing a list of paths to neutral proxy .fa files and a directory containing those neutral proxy .fa files
   
       Example:
```bash
ls  data/
chr10.fa           chr13.primate.maf  chr17.fa           chr1.primate.maf   chr2.fa           chr5.primate.maf  chr9.fa           ncHAE.v2.bed
chr10.primate.maf  chr14.fa           chr17.primate.maf  chr20.fa           chr2.primate.maf  chr6.fa           chr9.primate.maf
chr11.fa           chr14.primate.maf  chr18.fa           chr20.primate.maf  chr3.fa           chr6.primate.maf  chrX.fa           thurman.bed
chr11.primate.maf  chr15.fa           chr18.primate.maf  chr21.fa           chr3.primate.maf  chr7.fa           chrX.primate.maf
chr12.fa           chr15.primate.maf  chr19.fa           chr21.primate.maf  chr4.fa           chr7.primate.maf  chrY.fa
chr12.primate.maf  chr16.fa           chr19.primate.maf  chr22.fa           chr4.primate.maf  chr8.fa           chrY.primate.maf
chr13.fa           chr16.primate.maf  chr1.fa            chr22.primate.maf  chr5.fa           chr8.primate.maf  hg19
```

 4. ```./adaptiphy-launch-slurm.py ``` (optional): update this script if you are planning on using SLURM as a job manager to run the AdaptiPhy snakemake (preferred).
    * modify the header of this file to point to your snakemake conda env and email.
       Example:

```bash
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=email@university.edu
#SBATCH -N 1
#SBATCH --account=sciencelab
#SBATCH --partition=common
#SBATCH --mem=10G
#SBATCH -J adaptiphy
#SBATCH --time=3-00:00:00

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake \
  --profile slurm_general \
  --use-conda \
  --conda-prefix /path/to/your/conda/directories \
  --keep-going
```
     
 4. ```slurm_general/config.yaml``` (optional): update this file if you are planning on using SLURM as a job manager to run the AdaptiPhy snakemake (preferred). Do not modify this file's name or relative directory location.
    * update your ```slurm_partition``` and ```slurm_account``` to point to the correct partition and account.
    * set the ```tmpdir``` path to point to a scratch or work directory if available/desired.
    * the ```latency_wait``` and ```use-conda``` variables should not be altered. The other parameters can be optimized for your job scheduler system and memory needs.


## Generating data ##

### Running a local vs. global AdaptiPhy test ###

AdaptiPhy can be run in two modes.

**Local mode**: This is appropriate for performing a local scan for selection, e.g. examining one chromosome or a section of one chromosome. To run this mode, your chromosome target in the `config.yaml` file should point to only one chromosome, and the `neutral_set` parameter should additionally be set to `goodalignments.txt`. The neutral proxy in local mode is a series of concatenated random regions from your target window. This test is good for determining whether a large region that is putatively under selection is actually homogenously experiencing a higher mutation rate, etc. or whether there is truly a specific small region driving the signal selection. It is also good for detecting selection if your genome is very small (e.g. viral) and/or no neutral set is available.

**Global mode**: This is appropriate for performing a genome-wide scan for selection. To run this mode, your chromosome target in the `config.yaml` file should point to a subset of chromosomes or the entire list of chromosomes in your genome. The `neutral_set` parameter should additionally be set to a file such as `data/neutralset.txt`, which contains paths to a series of files in a `neutral_set` directory. To build this directory and file, proceed to the "Building an appropriate neutral proxy set" section. This test is good for scanning an entire large genome for regions under selection relative to a background genome-wide neutral evolution rate.

### Building an appropriate neutral proxy set ###

**running separate snakemake to generate  neutral set file (second snakemake)**

### Running the pipeline ###

**Option A:** running the snakemake pipeline as a batch job on SLURM
* It is typically a good idea to do a dry run of the snakemake pipeline interactively (`snakemake -n --latency-wait=300 --use-conda`) to check for syntax errors before submitting the final job to a job scheduler. 
* When you're ready, submit your batch job from the top-level directory as:

  ```bash
   sbatch slurm-launch-snakemake.sh
   ```
  
* This will source parameters for daughter jobs from `slurm_general/config.yaml`. The resulting SLURM log file will be written to `slurm.test.%j.out`. Individual rule logs will primarily be written to `logs/`.

**Option B:** running the snakemake pipeline interactively
* If you are not using a job scheduler and/or wish to run snakemake interactively, first log in to a virtual machine or other processor with sufficient memory to run the pipeline efficiently. Then run the pipeline from the top-level directory as:
   ```bash
   snakemake --latency-wait=300 --use-conda
   ```

**Option C:** running the snakemake pipeline on other job schedulers
 * we currently don't have support for this option. If you'd like to explore this option, check out the documentation for snakemake on other cluster systems [here] (https://snakemake.readthedocs.io/en/v5.6.0/executable.html) with more examples [here] (https://github.com/snakemake-profiles/doc).

### AdaptiPhy output ###

If the snakemake pipeline completes with no errors, your file structure should look something like this:
```bash
ls
config.yaml data/ DONE_SUMMARY.txt HYPHY/ intermediate_files/ logs/ OUTPUT_FINAL/ PhyloFit/ scripts/ slurm_general/ slurm-launch-snakemake.sh slurm.test.1234567.out Snakefile
```
The pipeline has generated the files `DONE_SUMMARY.txt HYPHY/ intermediate_files/ logs/ OUTPUT_FINAL/ PhyloFit/ slurm.test.1234567.out`. 
* _The not important stuff_: The `DONE_SUMMARY.txt` file simply contains the list of files in this directory after the final cleanup step. The `slurm.test.1234567.out` file will only exist if you ran the snakemake as a batch job on SLURM, and contains the breakdown of submitted SLURM jobs with information about step success, step order, and log file locations. The `logs/` directory contains logs from each individual rule run in the pipeline. the `intermediate_files` directory contains all intermediate files generated in the pipeline, which may be helpful for troubleshooting.
* _The semi-important stuff_: this snakemake runs the HyPhy as well as PhyloFit programs. All of the typical output files from these two tools are deposited in `HYPHY/` and `PhyloFit/`, respectively.
* _The important stuff_: the directory `OUTPUT_FINAL/` contains the formatted output tables from AdaptiPhy. The major results table is stored in `merged_summary_table.txt`.

__provide details here on how to interpret this table!__

### Resetting AdaptiPhy to run again ###

If you'd like to rerun AdaptiPhy again, make sure to remove the following files: `DONE_SUMMARY.txt HYPHY/ intermediate_files/ logs/ PhyloFit/`. You can also remove the `slurm.test.1234567.out` if you wish.

DO NOT remove the following directories and files unless you know what you're doing: `data/ config.yaml scripts/ slurm_general/ slurm-launch-snakemake.sh Snakefile`. Remove the `OUTPUT_FINAL` directory only if you're sure you don't need its contents anymore! 

### Citation
If you use this pipeline, please cite:
Berrio, A., Haygood, R. & Wray, G.A. Identifying branch-specific positive selection throughout the regulatory genome using an appropriate proxy neutral. BMC Genomics 21, 359 (2020). https://doi.org/10.1186/s12864-020-6752-4

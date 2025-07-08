# AdaptiPhy: Implementation with snakemake

This code updates the existing AdaptiPhy 2.0 pipeline to run using snakemake, allowing the user to plug in data at the beginning of the pipeline & performing the intervening steps automatically.

## Installing the AdaptiPhy pipeline to run with snakemake ##

### 1. Clone the git repo ###

To clone this repo from the command line into your working directory, use:

```git clone https://github.com/wodanaz/adaptiPhy```

You will need to add your files to the ```data/``` directory before running the snakemake pipeline for the first time.

### 2. Loading a conda env ###
The majority of the conda packages required in the pipeline will be loaded automatically, and you will not need to do any sort of manual install. You will only need to create an environment that contains snakemake and python in order to run this pipeline.

* Dependencies:

  * Anaconda/conda 24.x.x
  * python 3.x.x
  * snakemake 9.x.x
  * Unix/linux environment
  * Slurm (optional)

To create a fresh conda env for this purpose from the command line with Anaconda, use something like:

```conda create --name <my-env> python=3.6 snakemake```

Although this is a minimal environment, it's a good idea to store a .yml copy of your environment somewhere for replicability. To do this, activate your new environment and run:

``` conda env export --from-history > <yourname>.yml```

### 3. Change the necessary file parameters ###

To run the snakemake pipeline either interactively or through a job manager like SLURM, you will need to update some file paths and other information in the cloned repository.

1. ```./config.yaml```: This is the file where you will need to update the most information. You will need:
  * ```windows```: a list of ATAC peaks or similar targets in BED file format
  * ```tree_topology``` and ```foreground_branches```: provide a phylogenetic tree in Newick format and specify which branches are focal. For the first parameter, use standard Newick format. For the second, provide a vector of the focal branch names from the Newick tree.
  * ```maf_pattern``` and ```fa_pattern```: provide paths to files (wildcards permitted) for one or many chromosomes. Note that both the .maf (multi alignment) and .fa (nucleotide) files are required.
  * ```neutral_set```: provide a path to a .txt file that contains paths to neutral proxy files, or set this parameter to "goodalignments.txt" if running AdaptiPhy in local mode. More on this later!
  * ```chromosomes```: provide a vector of chromosomes to examine.

## Generating data ##

### Running a local vs. global AdaptiPhy test ###

### Building an appropriate neutral proxy set ###
running separate snakemake to generate  neutral set file (second snakemake)

### Running the pipeline ###
sbatch the pipeline or run interactively if on a non-slurm cluster (provide link to instructions for snakemake on other clusters)

### AdaptiPhy output ###
understanding the file structure and output table

### Citation
If you use this pipeline, please cite:
Berrio, A., Haygood, R. & Wray, G.A. Identifying branch-specific positive selection throughout the regulatory genome using an appropriate proxy neutral. BMC Genomics 21, 359 (2020). https://doi.org/10.1186/s12864-020-6752-4

# AdaptiPhy: Implementation with snakemake

This code updates the existing AdaptiPhy 2.0 pipeline to run using snakemake, allowing the user to plug in data at the beginning of the pipeline & performing the intervening steps automatically.

## Installing the AdaptiPhy pipeline to run with snakemake ##

### 1. Clone the git repo ###

To clone this repo from the command line into your working directory, use:

```bash
git clone https://github.com/wodanaz/adaptiPhy
```

You will need to add your files to the ```data/``` directory before running the snakemake pipeline for the first time. Read on for more info about the necessary file structure in this folder!

### 2. Loading a conda env ###
The majority of the conda packages required in the pipeline will be loaded automatically, and you will not need to do any sort of manual install. You will only need to create an environment that contains snakemake and python in order to run this pipeline.

* Dependencies:

  * Anaconda/conda 24.x.x
  * python 3.x.x
  * snakemake 9.x.x
  * Unix/linux environment
  * Slurm (optional)

To create a fresh conda env for this purpose from the command line with Anaconda, use something like:

```bash
conda create --name <my-env> python=3.6 snakemake
```

Although this is a minimal environment, it's a good idea to store a .yml copy of your environment somewhere for replicability. To do this, activate your new environment and run:

```bash
conda env export --from-history > <yourname>.yml
```

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
      windows: "data/allpeaks2.bed"
      tree_topology: "(Lv, (Ht, He))"
      foreground_branches: ["Ht", "He"]
      maf_pattern: "data/maf_files/{chrom}.He.maf"
      fa_pattern: "data/maf_files/{chrom}.He.fa"
      neutral_set: "data/neutralset.txt"
      chromosomes: ["chr1"]
   ```
    
 3. ```data/```: your input data lives in this folder. To run the AdaptiPhy pipeline, this folder must contain:
    * a folder containing MAF and .fa files, matching the specified 'pattern' paths in your ```config.yaml``` file from the previous step
    * a file of target windows/peak calls, matching the 'windows' path in  your ```config.yaml``` file from the previous step
    * if running AdaptiPhy in global mode (more on this later), a .txt file containing a list of paths to neutral proxy .fa files and a directory containing those neutral proxy .fa files
   
       Example:
       ```bash
        ls data/
        allpeaks.bed maf_fa/ neutral_set/ neutralset.txt
       ```

 4. ```./slurm-launch-snakemake.sh``` (optional): update this script if you are planning on using SLURM as a job manager to run the AdaptiPhy snakemake (preferred).
    * modify the header of this file to point to your snakemake conda env and email. Example:

      ```bash
       #SBATCH --mail-user=apm58@duke.edu
       ...
       source activate /path/to/conda/envs/env_name
       ```
     
 4. ```slurm_general/config.yaml``` (optional): update this file if you are planning on using SLURM as a job manager to run the AdaptiPhy snakemake (preferred). Do not modify this file's name or relative directory location.
    * update your ```slurm_partition``` and ```slurm_account``` to point to the correct partition and account.
    * set the ```tmpdir``` path to point to a scratch or work directory if available/desired.
    * the ```latency_wait``` and ```use-conda``` variables should not be altered. The other parameters can be optimized for your job scheduler system and memory needs.


## Generating data ##

### Running a local vs. global AdaptiPhy test ###
goodalignments vs neutral_set

### Building an appropriate neutral proxy set ###
running separate snakemake to generate  neutral set file (second snakemake)

### Running the pipeline ###
sbatch the pipeline or run interactively if on a non-slurm cluster (provide link to instructions for snakemake on other clusters)

### AdaptiPhy output ###
understanding the file structure and output table

### Citation
If you use this pipeline, please cite:
Berrio, A., Haygood, R. & Wray, G.A. Identifying branch-specific positive selection throughout the regulatory genome using an appropriate proxy neutral. BMC Genomics 21, 359 (2020). https://doi.org/10.1186/s12864-020-6752-4

# AdaptiPhy: Implementation with snakemake

This code updates the existing AdaptiPhy 2.0 pipeline to run using snakemake, allowing the user to plug in data at the beginning of the pipeline & performing the intervening steps automatically.

## Installing the AdaptiPhy pipeline to run with snakemake ##

### 1a. Recommended: Download the tarball release from the "Releases" page ###

add information here when this option is ready

### 1b. Not recommended: Clone the git repo ###

To clone this repo from the command line into your working directory, use:

```bash
git clone https://github.com/wodanaz/adaptiPhy
```

You will need to add your files to the ```data/``` directory before running the snakemake pipeline for the first time. Read on for more info about the necessary file structure in this folder!

### 2. Loading a conda env ###
The majority of the conda packages required in the pipeline will be loaded automatically, and you will not need to do any sort of manual install. You will only need to create an environment that contains snakemake and python in order to run this pipeline.

* Dependencies:

  * Miniconda/conda 24.9.2
  * python 3.11.7
  * snakemake 9.1.6
  * Unix/linux environment
  * Slurm (optional)

To create a fresh conda env for this purpose from the command line with conda, use something like:

```bash
touch snakemake.yml
```

Your file should look like this:

```
name: snakemake
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - snakemake=9.1
  - python=3.11
```

Now run:

```bash
conda env create -f snakemake.yml --prefix path/to/your/envs/folder/snakemake/
conda activate snakemake
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
    * modify the header of this file to point to your snakemake conda env and email.
       Example:

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
* The pipeline has generated the files `DONE_SUMMARY.txt HYPHY/ intermediate_files/ logs/ OUTPUT_FINAL/ PhyloFit/ slurm.test.1234567.out`. 
* _ _The not important stuff_ _: The `DONE_SUMMARY.txt` file simply contains the list of files in this directory after the final cleanup step. The `slurm.test.1234567.out` file will only exist if you ran the snakemake as a batch job on SLURM, and contains the breakdown of submitted SLURM jobs with information about step success, step order, and log file locations. The `logs/` directory contains logs from each individual rule run in the pipeline. the `intermediate_files` directory contains all intermediate files generated in the pipeline, which may be helpful for troubleshooting.
* _ _The semi-important stuff_ _: this snakemake runs the HyPhy as well as PhyloFit programs. All of the typical output files from these two tools are deposited in `HYPHY/` and `PhyloFit/`, respectively.
* _ _The important stuff_ _: the directory `OUTPUT_FINAL/` contains the formatted output tables from AdaptiPhy. The major results table is stored in `merged_summary_table.txt`.

__provide details here on how to interpret this table!__

### Citation
If you use this pipeline, please cite:
Berrio, A., Haygood, R. & Wray, G.A. Identifying branch-specific positive selection throughout the regulatory genome using an appropriate proxy neutral. BMC Genomics 21, 359 (2020). https://doi.org/10.1186/s12864-020-6752-4

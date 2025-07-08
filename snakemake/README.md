# AdaptiPhy: Implementation with snakemake

This code updates the existing AdaptiPhy 2.0 pipeline to run using snakemake, allowing the user to plug in data at the beginning of the pipeline & performing the intervening steps automatically.

## Installing the AdaptiPhy pipeline to run with snakemake ##
To clone this repo from the command line into your working directory, use:

```git clone https://github.com/wodanaz/adaptiPhy```

You will need to add your files to the ```data/``` directory before running the snakemake pipeline for the first time.

### Loading a conda env ###
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

### Change the necessary file parameters ###
list of files that require changes

### Generating an appropriate neutral proxy set ###
running separate snakemake to generate  neutral set file (second snakemake)

### Running the pipeline ###
sbatch the pipeline or run interactively if on a non-slurm cluster (provide link to instructions for snakemake on other clusters)

### AdaptiPhy output ###
understanding the file structure and output table

### Citation
If you use this pipeline, please cite:
Berrio, A., Haygood, R. & Wray, G.A. Identifying branch-specific positive selection throughout the regulatory genome using an appropriate proxy neutral. BMC Genomics 21, 359 (2020). https://doi.org/10.1186/s12864-020-6752-4

# Adaptiphy Slurm Pipeline

The Adaptiphy Slurm pipeline is based on [legacy.md](https://github.com/wodanaz/adaptiPhy/legacy.md). 

## Requirements
- [bash](https://www.gnu.org/software/bash/)
- [Slurm cluster](https://slurm.schedmd.com/) with shared storage
- [conda](https://conda.io/projects/conda/en/latest/index.html)  

_By default __conda__ is provided via the [Anaconda3/2019.10-gcb02](https://github.com/Duke-GCB/helmod/blob/master/rpmbuild/SPECS/Anaconda3-2019.10-gcb02.spec) environment module. This module name can be overridden via the `ANACONDAMODULE` environment variable._

### Staging Data Requirements
It is required to have a set of putatively neutral alignments and masked genome-wide alignment of the species of interest

### Conda Environment
The pipeline requires a conda environment named `escapevariants`.
See [environment.yml](environment.yml) for details.

## Installation
Clone this repository onto a shared location in your Slurm cluster. 
Create the config.sh file containing global configuration settings such as the slurm account.
Then create the `adaptiphy` conda environment.

### Create config.sh
The config.sh file can be created by copying the example file in place.
```
cp example-config.sh config.sh
```

### HARDAC Installation Instructions
On HARDAC the Anaconda3 module provides the conda command.
From an interactive session you can use the Anaconda3 module to create the `adaptiphy` conda environment by running:
```
module load Anaconda3/2019.10-gcb02
conda env create -f environment.yml
```

### General Installation Instructions
Using an Anaconda or Miniconda installallation create the `adaptiphy` conda environment by running:
```
conda env create -f environment.yml
```

## Setup
Before the pipeline can be run the input genome-alignment MAF files must be processed by [extract_alignments.sh](https://github.com/wodanaz/x/x/x/x/setup-xxxxx.sh).
The setup script will extract the query alignments and create its respective concatenated genome references used by the adaptiphy pipeline. The output files will be created in different directories named query and ref.
For example, to process a genome-wide aligment file named `genomewide_alignment.maf` run:
```
./setup-escape-variants.sh -g genomewide_alignment.maf -b file.bed
```

## Running
The `run-adaptiphys.sh` script is used to run the pipeline.
By default this script will stage input data and upload results, but this functionality can be disabled.

The script requires the following arguments:
1. query alignment files
2. reference files
3. a data directory that will contain input and output files
4. a project name - controls file naming and where to upload and download data from 
6. a tree in wick format 

The data directory will have two subdirectories created to hold various files.
- `/input` - Contains a subdirectory for each input project containing input files.
- `/output` - Contains a subdirectories for each output project containing output files created by the pipeline.


### Run including staging data
Given 
- Use "NC_045512.fasta" as your genome
- Use directory "data" to as your data directory
- DDS project named "CVExample" containing *.fastq.gz files in the top level directory
- Run in "duke" mode
- A datetab(TSV) file named date.tab
- Receive an email at username@example.com when the entire pipeline completes

The pipeline can be run using the NC_045512.fasta genome like so:
```
./run-escape-variants.sh -g NC_045512.fasta -d data -i CVExample -m duke -D date.tab -e username@example.com
```
The above command will do the following
1. Create the following directories
   - `data/input/CVExample` - the files/folders from the input CVExample project will be downloaded here   
   - `data/output/CVExample` - output files created by the pipeline will be created here
   - `data/output/CVExample/logs` - log files created running the pipeline will be created here          
2. Download the data in project named "cv-example" into the `data/input/CVExample` directory.
3. Run the pipeline creating output files in `data/output/CVExample`.
3. Upload the results from `data/output/CVExample` to an output project named CVExample_results.
4. Send an email to username@example.com when the pipeline completes.


### Run without staging data
To run without staging data you must ensure the input fastq.gz files are in a project specific directory.  
So if you want to use "sarscv" as your project name and use "data" as your data directory place your input *.fastq.gz files in a directory named `data/input/sarscv".
Run the pipeline skipping the download and upload steps like so: 
```
./run-escape-variants.sh -g NC_045512.fasta -d data -i sarscv -m duke -D date.tab -s -S
```
The `-s` argument causes run-escape-variants.sh to skip the download input data step.
The `-S` argument causes run-escape-variants.sh to skip the upload output data step.


#### Debugging
The pipeline creates a temp directory that contains all intermediate files. By default this directory is deleted when the pipeline completes successfully.
If you want to preserve the temp directory pass the `-k` argument.

Example:
```
./run-escape-variants.sh -g NC_045512.fasta -d data -i jpb-cv-example2 -m duke -D date.tab -k
```

#### Help
To see command line help by running `./run-escape-variants.sh` without arguments.
There are additional command line arguments for controlling where the logs and temp directory are created.

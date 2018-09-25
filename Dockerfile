FROM continuumio/anaconda

RUN apt-get update && apt-get install -y curl

## Install phast
WORKDIR /tmp

RUN curl -SLO http://compgen.cshl.edu/phast/downloads/phast.v1_5.x86_64.deb
RUN dpkg -i phast.v1_5.x86_64.deb

# Configure bioconda channel per https://bioconda.github.io/index.html

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# phast does two things. generate initial fasta files and substitution rate.

# Now install hyphy from bioconda
RUN conda install hyphy

# Install nano
RUN apt-get install -y nano

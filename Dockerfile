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

# Install nano, awk, sed
RUN apt-get install -y nano
RUN apt-get install -y gawk
RUN apt-get install -y sed 
RUN apt-get install -y util-linux 


#Install R and python

#From node:4
RUN apt-get update && apt-get remove -y python && apt-get install -y python2.7 r-base


# Add data
ADD /query/ /home/query/
ADD /ref/ /home/ref/
ADD /test/ /home/test/


# Define directory

WORKDIR /home/test


# Define maintainer

MAINTAINER Alejandro Berrio alebesc@gmail.com

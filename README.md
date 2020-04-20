# README

This README describes the minimal steps that are necessary to get an analysis of selection up and running.

### What is this repository for? ###

* The scripts associated with this walkthrough can be used to pulldown multiple alignments and analyse positive selection in specific branches.
* Version 2.0

### How do I get set up? ###

* Dependencies:

1. PHAST
1. HYPHY
1. BEDTOOLS
1. BIOPYTHON
1. PYTHON 2.7
1. Unix/linux environment
1. Slurm (optional)


# Selection test walkthrough
											
This walkthough can be used to test for positive selection in regulatory regions 	
from data obtained from functional genomic experiments such as ATAC-seq or  	
ChIP-seq.		

** note: For simplicity, you can download the non-functional regions 'refmasked2.tar.gz', the unmasked maf alignemtns as 'unmasked1.tar.gz', and the masked genomewide alignments for extracting local NFR extraction 'masked1.tar.gz'. After downloading,  just decompress with tar -xzvf and these directories are ready to use for making either local, global concatenated reference sequences, and to extract query regions

** if you download these files, start in step 4 ###





### 1 Download genomewide alignments  ### 
Download a multiZ alignment from UCSC (100way multiple alignment)

```bash
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; 
do wget --timestamping 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz100way/maf/'$chr.maf.gz; done
```

### 2 Download Referennce genome ###
Download the human assembly reference for each MAF from UCSC

```bash
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; 
do wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/'$chr.fa.gz; done ```
```

#### note: for the originl identification of NFRs, we masked the maf and human reference using known annotations. The following link shows how to mask the genome with different features. https://github.com/wodanaz/adaptiPhy/blob/master/Masking_MAF.md

### 3 Subset the 100-way genome-wide alignment for queries ###
To extract a smaller set of MAF files to use as query. This should correspond to the species of interest. In this case: human, chimp, gorilla, orangutan, macaque

```bash
for file in chr*.maf; do root= basename $file .maf; maf_parse $file --seqs hg19,ponAbe2,gorGor3,panTro4,rheMac3 --out-root primate/$root; done
```
Note: Save in a directory called unmasked


### 4 Extract query alignments ###
Now, we can pull down query alignments into multiple fasta files according to your features directory. This is the list of locations of to scan for selection


```bash
mkdir features
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY ; 
	do grep -w $chr yourfeatures.bed  | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V >  features/$chr.feat.bed; 
done
```

To generate alignments for the query regions

```bash
mkdir query
cd /your/directory/with/unmasked/maf/files


for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;
do msa_split $chr.primate.maf --refseq $chr.fa --gap-strip ANY -q --in-format MAF --features /your/test/directory/features/$chr.feat.bed --for-features 
--out-root/your/test/directory/query/$chr; 
done
```

HyPhy doesn't perform well with missing dashes, ambiguities or missing sequence. To remove these sequence features please run the included biopython script

### Note: Modify loading biopython according to your cluster needs 

### Note: if you decide to run local tests of selection, the procedure is a bit different. Please follow the pipeline in:

https://github.com/wodanaz/adaptiPhy/blob/master/Local_Reference.md



```bash

for file in *fa; do echo $file >> all.list;done


module load Anaconda/1.9.2-fasrc01
python prunning.py 

```

Once this is done, please generate a new list of fasta files

```bash
for file in *prunned; do echo $file >> all.prunned.list;done
```

Next, remove anything alignment smaller than 250 bp. Modify if necessary

```bash
module load Anaconda/1.9.2-fasrc01
python filtering.py 
```

### 6 Concatenate NFRs for neutral reference ###
Now, we need to make a putatively neutral reference


```bash
cat query/goodalignments.txt > queries.list

cp queries.list /your/path/to/neutral/alignments/from/masked/fasta

cd  /your/path/to/neutral/alignments/from/masked/fasta
```

Run these commands to generate neutral reference sequences. But make sure you have a set of putatively neutral sequences.
Please refer to findingneutrality.md in order to do this.


```bash
module load Anaconda/1.9.2-fasrc01
python DictGen.py
```



```bash
for file in `cat queries.list`; do echo $file.ref >> ref.tab; done
```

Move alignments to your reference directory

```bash
for file in `cat ref.tab`; do mv $file /your/path/to/reference/alignments ; done
```


### note: python adds a description that should be cleaned because hyphy will experience troubles at running. To do that, please run

```bash
for file in `cat ref.tab`; do root=`basename $file .fa.prunned.ref`; 
awk '{if($1 ~ /^>/){split($1,a,"\t"); print a[1]}else{print}}' $file > $root.ref.fa; done
```

```bash
for file in `cat  ref.tab`; do echo $file >> ref2.list; done
awk -F"." '{print $1 "." $2 }' ref2.list > ref.list


for file in `cat ref.tab`; do rm $file; done
```

Make a list of all batch files that have been generated for the alternative and the null models. This can be split for runs in parallel


```bash
for file in *hg19.null.bf; do echo $file >> null.hg19.list; done
for file in *hg19.alt.bf; do echo $file >> alt.hg19.list; done
```


If you have thousands of alignments to run. Please use this biopython script to generate a shell script to create a set of launchers to hyphy

```bash
module load Anaconda/1.9.2-fasrc01
python shgenerator.py
```

Alternatively:
```bash
for file in `cat null.hg19.list`; do
	root=`basename $file .hg19.null.bf`;
	HYPHYMP $file > HYPHY/$root.hg19.null.out;
done #exit nano ctrl+O ENTER ctrl+x
```


### 7 Run adaptiPhy using HyPhy ###
Now, we are ready to run hyphy


```bash
mkdir HYPHY
mkdir res

module load hyphy
for file in hyphy*sh ; do sbatch $file ; done
```




To extract P-values for each LRT, we need to pull down each likelihood from the null and alternative tests and construct a table. 




```bash
for file in *hg19.null.res; do echo $file >>null.hg19.log; done;
for file in *hg19.alt.res; do echo $file >>alt.hg19.log; done;

for filename in `cat alt.hg19.log`; do grep -H "BEST LOG-L:" $filename >> alt.hg19.tab; done;
for filename in `cat null.hg19.log`; do grep -H "BEST LOG-L:" $filename >> null.hg19.tab; done;

for file in *panTro4.null.res; do echo $file >>null.panTro4.log; done;
for file in *panTro4.alt.res; do echo $file >> alt.panTro4.log; done;

for filename in `cat alt.panTro4.log`; do grep -H "BEST LOG-L:" $filename >> alt.panTro4.tab; done;
for filename in `cat null.panTro4.log`; do grep -H "BEST LOG-L:" $filename >> null.panTro4.tab; done;
```


Consolidating tables for null and alternative models

```bash
awk '{print $1 "\t" $3}' null.hg19.tab | sort -k1,1 -V  > nulls.hg19.tab
awk '{print $1 "\t" $3}' alt.hg19.tab | sort -k1,1 -V > alts.hg19.tab

awk '{print $1 "\t" $3}' null.panTro4.tab | sort -k1,1 -V  > nulls.panTro4.tab
awk '{print $1 "\t" $3}' alt.panTro4.tab | sort -k1,1 -V > alts.panTro4.tab
```


```bash
awk -F"." '{print $1 ":" $2 "\t" $3 }'  nulls.hg19.tab > col1.hg19.tab
paste col1.hg19.tab nulls.hg19.tab  | awk '{print $1 "\t" $2 "\t" $4   }' > lnulls.hg19.tab
paste lnulls.hg19.tab alts.hg19.tab | awk '{print $1 "\t" $2  "\t" $3 "\t" $5  }' |  sort -k1,1 -V > likelihoods.hg19.tab
 
awk -F"." '{print $1 ":" $2 "\t" $3 }'  nulls.panTro4.tab > col1.panTro4.tab
paste col1.panTro4.tab nulls.panTro4.tab  | awk '{print $1 "\t" $2 "\t" $4   }'> lnulls.panTro4.tab
paste lnulls.panTro4.tab alts.panTro4.tab | awk '{print $1 "\t" $2  "\t" $3 "\t" $5  }' |  sort -k1,1 -V > likelihoods.panTro4.tab
 


cat likelihoods.hg19.tab likelihoods.panTro4.tab  |  sort -k1,1 -V | sed 1i"chromosome\tbranch\tlnull\tlalt" > likelihoods.tab

```
Here we should have a table with four columns of data (Genomic region, branch, null and alternative likelihoods). Next, we will add zeta and other sequence features

To get the p-value using a chi^2 function in R. Please execute LRT.R

```bash
Rscript LRT.R
```

This will create a new table containing the P values. This file is called "likelihoods.pvals.tab"

```bash
awk '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' likelihoods.pvals.tab | sort -k1,1 -V | sed 1i"chromosome\tbranch\tlnull\tlalt\tpval"  > P_value.data
```

To calculate zeta we compute substitution rates at each branch of the tree in both query and its corresponding reference.

For queries (modify according to the species in your tree):

```bash
cd query
mkdir -p MODELS_HKY85
for file in `cat goodalignments.txt` ; 
do root=`basename $file .ref.fa`; 
phyloFit $file --tree "(taxon0,(taxon1,(taxon2,(taxon3, taxon4))))" --subst-mod HKY85 --msa-format FASTA --out-root MODELS_HKY85/$root; 
done 
```

For the reference:

```bash

ls *.ref.fa > ref.fa.list
cd ref
mkdir -p MODELS_HKY85


for file in `cat ref.fa.list` ; 
do root=`basename $file .fa`; 
phyloFit $file --tree "(taxon0,(taxon1,(taxon2,(taxon3, taxon4))))" --subst-mod HKY85 --msa-format FASTA --out-root MODELS_HKY85/$root; 
done 

```

When done, we need to transform tables in a format we can use:


In the query directory:

```bash
for filename in *.mod; do grep -H "TREE:" $filename; done > output.hky85.txt
cat output.hky85.txt | awk -F":" '{print $1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11  }' | \
awk -F"," '{print $1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10  }'  | \
awk -F")" '{print $1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10  }'  | \
awk '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $13 }'  > BranchLenghts.tab

# Let's stop here for a while.... we need to check the table and make sure it's looking good

# we need to further modify column 1 alone -> I don't like the dots, it would be nicer to have chromosome and location in separate columns

awk -F"." '{print $1 ":" $2  }' BranchLenghts.tab > chr_pos.tab

paste chr_pos.tab BranchLenghts.tab  | column -s '\t' -t | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }'  > Branches.tab


sed 1i"chromosome\trheMac3\tponAbe2\tgorGor3\tpanTro4\thg19\tPanHomo\tPanHomoGor\tPanHomoGorPon" Branches.tab > Q.hky85.Branches.tab
```

In the Reference directory

```bash
for filename in *.mod; do grep -H "TREE:" $filename; done > output.hky85.txt
cat output.hky85.txt | awk -F":" '{print $1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11  }' | \
awk -F"," '{print $1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10  }'  | \
awk -F")" '{print $1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10  }'  | \
awk '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $13 }'  > BranchLenghts.tab
# Let's stop here for a while.... we need to check the table and make sure it's looking good

# we need to further modify column 1 alone -> I don't like the dots, it would be nicer to have chromosome and location in separate columns

awk -F"." '{print $1 ":" $2  }' BranchLenghts.tab > chr_pos.tab

paste chr_pos.tab BranchLenghts.tab  | column -s '\t' -t | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }'  > Branches.tab


sed 1i"chromosome\trheMac3\tponAbe2\tgorGor3\tpanTro4\thg19\tPanHomo\tPanHomoGor\tPanHomoGorPon" Branches.tab > R.hky85.Branches.tab

```

Next, copy these two results in the main directory to merge t


```bash

cp query/MODELS_HKY85/Q.hky85.Branches.tab .
cp ref/MODELS_HKY85/R.hky85.Branches.tab .


wc -l Q.hky85.Branches.tab
wc -l R.hky85.Branches.tab



sed -i -e "1d" Q.hky85.Branches.tab
sed -i -e "1d" R.hky85.Branches.tab



sort Q.hky85.Branches.tab -k1,1 -V > Q.hky85.Branches.sorted.tab
sort R.hky85.Branches.tab -k1,1 -V > R.hky85.Branches.sorted.tab

awk '{print $1}' Q.hky85.Branches.sorted.tab | awk -F":" '{print $1 "\t" $2 }' | awk -F"-" '{print $1 "\t" $2 "\t" $3}' > Q.hky85.bed
awk '{print $1}' R.hky85.Branches.sorted.tab | awk -F":" '{print $1 "\t" $2 }'  | awk -F"-" '{print $1 "\t" $2 "\t" $3}' > R.hky85.bed


paste Q.hky85.bed Q.hky85.Branches.sorted.tab | column -s '\t' -t | awk '{print $1 "\t" $2 "\t" $3 "\t" $9 }' > Q.hky85.hg19.bed # human
paste R.hky85.bed R.hky85.Branches.sorted.tab | column -s '\t' -t | awk '{print $1 "\t" $2 "\t" $3 "\t" $9 }' > R.hky85.hg19.bed # human

paste Q.hky85.bed Q.hky85.Branches.sorted.tab | column -s '\t' -t | awk '{print $1 "\t" $2 "\t" $3 "\t" $8 }' > Q.hky85.panTro4.bed # chimp
paste R.hky85.bed R.hky85.Branches.sorted.tab | column -s '\t' -t | awk '{print $1 "\t" $2 "\t" $3 "\t" $8 }' > R.hky85.panTro4.bed # chimp





join -t $'\t' -j 1 -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9  Q.hky85.Branches.sorted.tab R.hky85.Branches.sorted.tab  > join_trans.tab


awk '{ rate1 = $6 / $14 ; print $1 "\t" $6 "\t"  $14  "\t" rate1 "\t" "hg19" }' join_trans.tab  > 	PhyloFit.hg19.tab
awk '{  rate2 = $5 / $13 ; print $1 "\t" $5 "\t"  $13  "\t" rate2 "\t" "panTro4"  }' join_trans.tab  > 	PhyloFit.panTro4.tab


cat PhyloFit.hg19.tab PhyloFit.panTro4.tab  |  sort -k1,1 -k5 -V | awk '{print $1 "\t" $4 "\t" $5 }'  | sed 1i"chromosome\tRatioPhyloFit\tbranch" > PhyloFit.data

```


To build consolidated table I prefer to use R to merge tables by column

```bash
module load R
R

data_pval = as.data.frame(read.table("P_value.data", header = F)) # read tab file 
head(data_pval)
colnames(data_pval) <- c("chromosome", "branch", "lnull", "lalt", "pval")
rates = as.data.frame(read.table("PhyloFit.data", header = T)) # read tab file   
head(rates)

selection.data <- merge(data_pval, rates, by= c("chromosome", "branch"), suffixes = c(".human",".chimp"))
selection.data$lnull <- NULL
selection.data$lalt <- NULL
write.table(selection.data, file ="selection.tab.data", row.names=F, col.names=T, quote=F) 



q()
``` 


A shortcut to transform the previous table and eliminate the branch column.


```bash
grep 'hg19' selection.tab.data | awk  '{ print $1 "\t" $2 "\t" $3 "\t" $4 }' > selection.hg19.data
grep 'panTro4' selection.tab.data | awk  '{ print $1 "\t" $2 "\t" $3 "\t" $4 }' > selection.panTro4.data
```


```R
selection_hg19 = as.data.frame(read.table("selection.hg19.data", header = F)) # read tab file 
selection_pantro4 = as.data.frame(read.table("selection.panTro4.data", header = F)) # read tab file 

selection.data <- merge(selection_hg19, selection_pantro4, by= "V1", suffixes = c(".human",".chimp"))

selection.data$V2.human <- NULL     # to remove branch column
selection.data$$V2.chimp <- NULL

colnames(selection.data) <- c('chromosome', 'pval.human', 'zeta.human', 'pval.chimp', 'zeta.chimp')


write.table(selection.data, file ="selection.data", row.names=F, col.names=T, quote=F) 


q()

```


### You can use a similar method to add any other sequence feature of interest. Such as: gene names, distance from closest TSS, CG content, etc. ###


# Congratulations!! Now you are ready to download this table to visualize and analyse your selection data.  


### Who do I talk to? ###

* Please contact me at alebesc@gmail.com

* Other community or team contact:
	Greg Wray gwray at duke.edu 
	Ralph Haygood ralph at ralphhaygood.org




#This repository also contains scripts that can be used to analyse sequence data


GC_content.py is a script that calculates GC content in the query sequences in the human and chimpanzee branches, as well, as the mean accross all branches. It requires biopython.

#module load Anaconda/1.9.2-fasrc01
#python GC_content.py

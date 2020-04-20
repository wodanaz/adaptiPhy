# This walkthrough can be used to mask the human genome and the primate MAF

Copy and transfer all the bed files with boundary information to this new directory


* knownRefSeq.bed: RefSeq Human Known genes 
* knownGene.bed: UCSC Human Known genes   
* lincRNA.bed: long non-coding RNAs
* tRNA.bed: Human tRNA
* wgRNA.bed: Human microRNA and sncRNA
* humanmRNA.bed: Total human mRNA

* OpenChromSynth.master.bed: Synthesis of all human open chromatin, DNAseI-seq, FAIRE-seq, ChIP-seq
* thurman.master.bed: Thurman et al, 2012. 125 different cell types and tissues
* honeybadger.master.bed: Curated list of promoters and enhancers from HoneyBadger https://personal.broadinstitute.org/meuleman/reg2map/HoneyBadger_release/
* vistaEnhancer.bed: Human Vista Enhancers

* CpGislands.bed: known for UCSC genome browser

* microsats.bed: short repeats
* InterruptedRpts.bed: Interrupted Rpts - Fragments of Interrupted Repeats Joined by RepeatMasker ID
* NumtS.bed: Human NumtS mitochondrial sequence
* simplerepeats.bed: Simple Tandem Repeats by TRF




 
```bash
cat 5_UTR_Exons.bed 3_UTR_Exons.bed lincRNA.bed humanmRNA.bed tRNA.bed wgRNA.bed  OpenChromSynth.master.bed thurman.master.bed honeybadger.master.bed vistaEnhancer.bed CpGislands.bed microsats.bed InterruptedRpts.bed NumtS.bed simplerepeats.bed | sort -k1,1 -k2,2n -V | awk '{ print $1 "\t" $2 "\t" $3 }' > functional.feat.bed
```



After concatenating many files, consider using refeature from PHAST, this tool avoids overlapping features.
We don't want to obtain repeated alignments.


```bash
bedtools merge -i functional.feat.bed > functional.feat.merge.bed

mkdir featBEDv3

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY ; 
	do grep -w $chr functional.feat.merge.bed  | awk '{print $1 "\t" $2 "\t" $3 }' > featBEDv3/$chr.feat.bed; 
done

cd ..
```


To sample random genomic regions that can be considered putatively neutral. 
We need to mask genomic regions that are 

To sample random genomic regions that can be considered putatively neutral. 
We need to mask genomic regions that are 


```bash
mkdir maskedv3
for primates in chr*.maf ; do
root=`basename $primates .primate.maf`;
maf_parse $primates --features functional_feat/featBEDv3/$root.feat.bed --mask-features hg19,ponAbe2,gorGor3,panTro4,rheMac3 > maskedv3/$root.masked.maf;   
done
```

#### Don't forget to mask the same genomic regions in the reference chromosomes as well


```bash
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
bedtools maskfasta -fi $chr.fa -bed functional_feat/featBEDv3/$chr.feat.bed -fo maskedv3/$chr.masked.fa;
done 
```


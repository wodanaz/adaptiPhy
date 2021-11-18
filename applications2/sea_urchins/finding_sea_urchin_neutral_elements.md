# Analysis of Positive Selection in Sea Urchins


To find neutral elements, we need to mask the genome-wide alignments (maf) using all known functional sequences and repetitive DNA

```bash

module load Anaconda3

# load your environment
conda activate aleconda

awk '{ print $1 "\t" $2 "\t" $3 }'  transcript_coords_He.txt  > genes.He.bed

cat genes.He.bed Hery_repeats.bed  | awk '{ print $1 "\t" $2 "\t" $3 }'  | sort -k1,1V -k2,2n >  functional.feat.bed

bedtools merge -i functional.feat.bed > functional.feat.merge.bed

mkdir featBED

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ; 
	do grep -w $chr functional.feat.merge.bed  | awk '{print $1 "\t" $2 "\t" $3 }' > featBED/$chr.feat.bed; 
done

```


```bash
mkdir masked_He_genome
for seaurchins in chr*.maf ; do
root=`basename $seaurchins .He.maf`;
maf_parse $seaurchins --features featBED/$root.feat.bed --mask-features He,Ht,Lv > masked_He_genome/$root.masked.maf;   
done


for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21; do
bedtools maskfasta -fi $chr.He.fa -bed featBED/$chr.feat.bed -fo masked_He_genome/$chr.masked.fa;
done 


```


To extract a set of non-functional elements that can be used as  putatively neutral loci


```bash

cd masked_He_genome

bedtools random -n 10000000 -l 300 -g sizes.genome > random.bed

sort -k1,1 -k2,2n random.bed > random.sorted.bed

awk '{ print $1 "\t" $2 - 1 "\t" $3 }' random.sorted.bed > random.300.bed 

wc -l random.300.bed # in masked version
5000000 


mkdir features

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ; 
	do grep -w $chr random.300.bed | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V > features/$chr.feat.bed; 
done

```

Note: to run the previous lines, it is necessary to have a genome file with the lenghts of each chromosome. To create it, do:

```bash

cat sizes.genome
chr1    78554062
chr2    71715399
chr3    56650573
chr4    56843197
chr5    52154090
chr6    54795018
chr7    49373989
chr8    49071437
chr9    53544062
chr10   55791632
chr11   60414877
chr12   51325886
chr13   39571473
chr14   46128125
chr15   42174613
chr16   43151226
chr17   37624446
chr18   41902230
chr19   37017750
chr20   37187759
chr21   48917100
unplaced_scaffold22     172537
unplaced_scaffold23     2000
unplaced_scaffold24     112705
unplaced_scaffold25     93463
unplaced_scaffold26     2000
unplaced_scaffold27     12905
unplaced_scaffold28     10217
unplaced_scaffold29     14520

```


# To create a subset list of putative neutral elements:

1. Extract Alignments from features

```bash

mkdir neutral_alignments

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21; do
echo '#!/usr/bin/env bash' > $chr.do_data.sh;
echo '#SBATCH --mail-type=END' >> $chr.do_data.sh;
echo '#SBATCH --mail-user=alebesc@gmail.com' >> $chr.do_data.sh;
echo '#SBATCH -N 1' >> $chr.do_data.sh;
echo "msa_split $chr.masked.maf  --refseq $chr.masked.fa --gap-strip ANY -q --in-format MAF  --features features/$chr.feat.bed --for-features --out-root  neutral_alignments2/$chr;"    >> $chr.do_data.sh;
done

for file in *do_data.sh ; do sbatch $file ; done


```

First, we need to filter and select the sequences with fully homologous sequences and > 200 bp long:


```bash


cd neutral_alignments2

nano makelists.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21; do
        echo $chr.*fa >> $chr.list;
	sed -ri 's/ /\n/g' $chr.list
done


sbatch makelists.sh

```

Number of regions extracted

```bash
wc -l chr*list


```

```bash
nano prunning.py
import re # for Regular expressions
import sys, sets

infile=sys.argv[1]
#infile = 'all.list'

with open(infile) as f:
    reflist = f.read().splitlines() 



f.close()


for i in reflist:
	myfile= open(i, 'r')
	data=myfile.read().splitlines()
	myfile.close()
	output=open('%s.prunned' % i, 'w')
	results={}
	seq=""
	compName=""
	for line in data:
		if line[0]==">":
			if compName!="":
				results[compName]=seq
			compName=line[1:]
			seq=""
		else:
			seq=seq+line
	results[compName]=seq
	nsize=len(results[results.keys()[0]])
	ns=len(results.keys())
	indels=sets.Set([])
	prunned={}
	for spec in results.keys():
		for i in range(nsize):
			base=results[spec][i]
			if base.upper() not in ["A","C","G","T"]:
				indels.add(i)
	nsize2=nsize-len(indels)
	for spec in results.keys():
		seq=""
		for i in range(nsize):
			if i not in list(indels):
				seq=seq+results[spec][i]
		prunned[spec]=seq
	for f in range(len(prunned.keys())):
		spec=prunned.keys()[f]
		if i>0: output.write("\n")
		output.write(">"+spec+"\n")
		output.write(prunned[spec])
	output.close()
	

		



######### Run this python script
### IMPORTANT  ------>>>> 
#To run this you need Python 2.7.6
```





```bash

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21; do
echo '#!/usr/bin/env bash' > $chr.prunning.sh;
echo '#SBATCH --mail-type=END' >> $chr.prunning.sh;
echo '#SBATCH --mail-user=alebesc@gmail.com' >> $chr.prunning.sh;
echo '#SBATCH -N 1' >> $chr.prunning.sh;
echo "python2 prunning.py ${chr}.list"    >> $chr.prunning.sh;
done

for file in *prunning.sh ; do sbatch $file ; done


```


```bash


nano makelists2.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ; do
        echo $chr.*prunned >> $chr.prun.list;
	sed -ri 's/ /\n/g' $chr.prun.list
done
```





```bash

wc -l chr*prun.list


cat  chr*prun.list > all.prunned.list

nano filtering_p3.py
import re # for Regular expressions
import sys
from Bio import AlignIO

infile=sys.argv[1]

#infile = 'all.prunned.list'

with open(infile) as f:
    reflist = f.read().splitlines() 



f.close()
missing = '[(+*)]'
ambiguous = '[N]'


badseqssNs = []
badseqsAsk = []
goodseqs = []

for i in reflist:
	myfile= open(i, 'r')
	mydata = myfile.read()
	fasta = AlignIO.read(i, "fasta")
	myfile.close()
	if fasta.get_alignment_length() > 200 : 
		maxNs = len(re.findall(ambiguous, str(mydata), re.I))
		maxAsk = len(re.findall(missing, str(mydata), re.I))
		if maxNs > 2:
			print(')-B')
			badseqsNs.append(i)
		elif maxAsk > 1:
			print(')-;')
			badseqsAsk.append(i)
		else:
			print('lol')
			goodseqs.append(i)
		




thefile3 = open('goodalignments.txt', 'w')
for item in goodseqs:
  thefile3.write("%s\n" % item)





thefile3.close()





thefile1 = open('ambiguous.txt', 'w')
for item in badseqsNs:
  thefile1.write("%s\n" % item)




thefile1.close()

thefile2 = open('asterisks.txt', 'w')
for item in badseqsAsk:
  thefile2.write("%s\n" % item)





thefile2.close()


######### Run this python script
# python3

```





```bash

nano filtering.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
python filtering_p3.py  all.prunned.list


sbatch filtering.sh



```

Next, we need to compute the substitution rate among each node and tip of the three. Here, we used PhyloFit to do this.

```bash

mkdir -p MODELS_HKY85
nano domodel.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
for file in *.fa.prunned ; 
do root=`basename $file .fa.prunned`; 
phyloFit $file --tree "(Lv,(Ht,He))" -i FASTA  --subst-mod HKY85 --out-root MODELS_HKY85/$root; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x



cd MODELS_HKY85

for filename in *.mod; do grep -H "TREE:" $filename; done > output.hky85.txt


sed -r 's/\.mod:TREE: \(/\t/g' output.hky85.txt | sed -r 's/:/\t/g' | sed -r 's/,\(/\t/g' | sed -r 's/,/\t/g' | sed -r 's/\);//g' |  sed -r 's/\)//g' |  awk '{ print $1  "\t" $3  "\t" $5  "\t" $7 "\t" $8  }' |  sed -r 's/\./:/' | sed 1i"genome_location\tLv\tHt\tHe\tHt_He" > urchin_branches.mask.data



```

### Next, we need to compute the relative substitution rate of each branch relative to the entire tree.
 
Here we use R to do this

```R
# set working directory to the path where you saved Branches.mask.data
setwd("~/path/to/file")

# Sum up the entire tree branch lenght

fulltree_br <- Genome_branches_mask$Lv  + Genome_branches_mask$Ht_He + Genome_branches_mask$Ht + Genome_branches_mask$He

# Compute relative branch length for each branch:

Genome_branches_mask$Lv_rb <- Genome_branches_mask$Lv / fulltree_br
Genome_branches_mask$He_rb <- Genome_branches_mask$He / fulltree_br
Genome_branches_mask$Ht_rb <- Genome_branches_mask$Ht / fulltree_br

# Plot density distribution of relative branch length
library(ggplot2)


ggplot() +
  geom_density(data =Genome_branches_mask, aes(x= Lv_rb),  colour="black", fill="darkblue") +
  geom_density(data = Genome_branches_mask, aes(x=He_rb),  colour="black", fill="green", alpha= 0.3) +
  geom_density(data =Genome_branches_mask, aes(x=Ht_rb),  colour="black", fill="red", alpha= 0.5) +
  theme_bw() + labs(x = "Branch length", y="Density", title = "Reference Branch Length Distribution (Masked)") +
  scale_x_continuous(limits = c(0, 0.5)) 

  
  
```

![alt text](https://github.com/wodanaz/adaptiPhy/blob/master/applications2/sea_urchins/distribution_plot.png)


```R

# remove trees with very small substitution rate. Change the value of the parameter (0.001) according to the distribution figure. The idea is to remove the peak of sites with a substitution rate near to zero (these are highly conserved sites).

# Remove trees with very small rates
NoMissData_mask <- subset(subset(subset(Genome_branches_mask, He > 0.001 ), Ht > 0.001 ), Lv > 0.001)

dim(NoMissData_mask)
dim(Genome_branches_mask)


ggplot() +
  geom_density(data =NoMissData_mask, aes(x= Lv_rb),  colour="black", fill="darkblue") +
  geom_density(data =NoMissData_mask, aes(x=He_rb),  colour="black", fill="green", alpha= 0.3) +
  geom_density(data =NoMissData_mask, aes(x=Ht_rb),  colour="black", fill="red", alpha= 0.5) +
  theme_bw() + labs(x = "Branch length", y="Density", title = "Reference Branch Length Distribution (Masked)") +
  scale_x_continuous(limits = c(0, 0.5)) 
```

![alt text](https://github.com/wodanaz/adaptiPhy/blob/master/applications2/sea_urchins/No_small_rates.png)

```R
# compute mean, median and first and third quantile for the human branch:


mean_He <- mean(NoMissData_mask$He_rb)
median_He <- median(NoMissData_mask$He_rb)

quantile(NoMissData_mask$He_rb)
q1_he <- matrix(summary(NoMissData_mask$He_rb))[2,]
q3_he <- matrix(summary(NoMissData_mask$He_rb))[5,]

q1_ht <- matrix(summary(NoMissData_mask$Ht_rb))[2,]
q3_ht <- matrix(summary(NoMissData_mask$Ht_rb))[5,]

neutralset <- subset(subset(NoMissData_mask, He_rb>q1_he & He_rb<q3_he), Ht_rb>q1_ht & Ht_rb<q3_ht)
dim(neutralset)


ggplot() +
  geom_density(data =NoMissData_mask, aes(x= Lv_rb),  colour="black", fill="darkblue") +
  geom_density(data =NoMissData_mask, aes(x=He_rb),  colour="black", fill="green", alpha= 0.3) +
  geom_density(data =NoMissData_mask, aes(x=Ht_rb),  colour="black", fill="red", alpha= 0.5) +
  theme_bw() + labs(x = "Branch length", y="Density", title = "Reference Branch Length Distribution (Masked)") +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_vline(xintercept = q1, size = 1, colour = "green",
                                                    linetype = "dashed") +
  geom_vline(xintercept = q3, size = 1, colour = "green",
             linetype = "dashed")


ggplot() +
  geom_density(data =neutralset, aes(x= Lv_rb),  colour="black", fill="darkblue") +
  geom_density(data =neutralset, aes(x=He_rb),  colour="black", fill="green", alpha= 0.3) +
  geom_density(data =neutralset, aes(x=Ht_rb),  colour="black", fill="red", alpha= 0.5) +
  theme_bw() + labs(x = "Branch length", y="Density", title = "Reference Branch Length Distribution (Masked)") +
  scale_x_continuous(limits = c(0, 0.5)) 


![alt text](https://github.com/wodanaz/adaptiPhy/blob/master/applications2/sea_urchins/neutral_rates.png)


# Finally, save this set of sequences. We assume these are non-functional and putatively neutral
neutralset.txt <- paste(neutralset$genome_location, "fa", sep=".")


write.table(neutralset.txt, file ="neutralset.txt", row.names=F, col.names=F, quote=F) 


```

With the previous pipeline, I identified  21318 putatively neutral elements that can be used for testing selection


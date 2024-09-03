# Extract Alignments
To generate REFERENCE alignments for using blind features located 
```bash
module load bedtools2

bedtools random -n 5000000 -l 300 -g hg19 > random.bed

sort -k1,1 -k2,2n random.bed > random.sorted.bed

awk '{ print $1 "\t" $2 - 1 "\t" $3 }' random.sorted.bed > random.300.bed 

wc -l random.300.bed # in masked version
5000000 


mkdir features

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY ; 
	do grep -w $chr random.300.bed | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V > features/$chr.feat.bed; 
done


```

# To create a subset list of putative neutral elements:

1. Extract Alignments from features

```bash

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
echo '#!/usr/bin/env bash' > $chr.do_data.sh;
echo '#SBATCH --mail-type=END' >> $chr.do_data.sh;
echo '#SBATCH --mail-user=your.email.com' >> $chr.do_data.sh;
echo '#SBATCH -N 1' >> $chr.do_data.sh;
echo "msa_split $chr.masked.maf  --refseq $chr.masked.fa --gap-strip ANY -q --in-format MAF  --features features/$chr.feat.bed --for-features --out-root  neutral_alignments/$chr;"    >> $chr.do_data.sh;
done

for file in *do_data.sh ; do sbatch $file ; done



```

First, we need to filter and select the sequences with fully homologous sequences and 300bp long:


```bash

nano makelists.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email.com
#SBATCH -N 1
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
        echo $chr.*fa >> $chr.list;
	sed -ri 's/ /\n/g' $chr.list
done
```

Number of regions extracted

```bash
wc -l chr*list
   216410 chr10.list
   218074 chr11.list
   215652 chr12.list
   185533 chr13.list
   173391 chr14.list
   165023 chr15.list
   146168 chr16.list
   134329 chr17.list
   130789 chr18.list
    94301 chr19.list
   403181 chr1.list
   104160 chr20.list
    75389 chr21.list
    82617 chr22.list
   390814 chr2.list
   320722 chr3.list
   307558 chr4.list
   293864 chr5.list
   276362 chr6.list
   257769 chr7.list
   234542 chr8.list
   223761 chr9.list
   252409 chrX.list
    93081 chrY.list
  4995899 total
```




```bash

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
echo '#!/usr/bin/env bash' > $chr.prunning.sh;
echo '#SBATCH --mail-type=END' >> $chr.prunning.sh;
echo '#SBATCH --mail-user=your.email.com' >> $chr.prunning.sh;
echo '#SBATCH -N 1' >> $chr.prunning.sh;
echo "python2 prunning.py ${chr}.list"    >> $chr.prunning.sh;
done

for file in *prunning.sh ; do sbatch $file ; done


```


```bash


nano makelists2.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email.com
#SBATCH -N 1
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
        echo $chr.*prunned >> $chr.prun.list;
	sed -ri 's/ /\n/g' $chr.prun.list
done


for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
echo '#!/usr/bin/env bash' > $chr.filtering.sh;
echo '#SBATCH --mail-type=END' >> $chr.filtering.sh;
echo '#SBATCH --mail-user=your.email.com' >> $chr.filtering.sh;
echo '#SBATCH -N 1' >> $chr.filtering.sh;
echo "python filtering_p3.py  ${chr}.prun.list" >> $chr.filtering.sh;
done

for file in *filtering.sh ; do sbatch $file ; done



```

Next, we need to compute the substitution rate among each node and tip of the three. Here, we used PhyloFit to do this.

```bash

mkdir -p MODELS_HKY85
nano domodel.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email.com
for file in *.fa ; 
do root=`basename $file .fa`; 
phyloFit $file --tree "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))" --subst-mod HKY85 --out-root MODELS_HKY85/$root; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x



cd MODELS_HKY85

#For Reference
for filename in *.mod; do grep -H "TREE:" $filename; done > output.hky85.txt


sed -r 's/\.mod:TREE: \(/\t/g' output.hky85.txt | sed -r 's/:/\t/g' | sed -r 's/,\(/\t/g' | sed -r 's/,/\t/g' | sed -r 's/\);//g' |  sed -r 's/\)//g' |  awk '{ print $1  "\t" $3  "\t" $5  "\t" $7 "\t" $9  "\t" $11 "\t" $12 "\t" $13 "\t" $14 }' |  sed -r 's/\./:/' | sed 1i"chromosome\trheMac3\tponAbe2\tgorGor3\tpanTro4\thg38\tPanHomo\tPanHomoGor\tPanHomoGorPon" > Branches.mask.data



```

### Next, we need to compute the relative substitution rate of each branch relative to the entire tree.
 
Here we use R to do this

```R
# set working directory to the path where you saved Branches.mask.data
setwd("~/path/to/file")

Genome_branches_mask = read.table("Branches.mask.data", header = TRUE)  # read tab file 

# Compute relative branch lenght for each branch:

rheMac3_rb <- Genome_branches_mask$rheMac3 / (Genome_branches_mask$rheMac3 + Genome_branches_mask$ponAbe2 + Genome_branches_mask$gorGor3 + Genome_branches_mask$panTro4 + Genome_branches_mask$hg19 + Genome_branches_mask$PanHomo+ Genome_branches_mask$PanHomoGor + Genome_branches_mask$PanHomoGorPon)
ponAbe2_rb <- Genome_branches_mask$ponAbe2 / (Genome_branches_mask$rheMac3 + Genome_branches_mask$ponAbe2 + Genome_branches_mask$gorGor3 + Genome_branches_mask$panTro4 + Genome_branches_mask$hg19 + Genome_branches_mask$PanHomo+ Genome_branches_mask$PanHomoGor + Genome_branches_mask$PanHomoGorPon)
gorGor3_rb <- Genome_branches_mask$gorGor3 / (Genome_branches_mask$rheMac3 + Genome_branches_mask$ponAbe2 + Genome_branches_mask$gorGor3 + Genome_branches_mask$panTro4 + Genome_branches_mask$hg19 + Genome_branches_mask$PanHomo+ Genome_branches_mask$PanHomoGor + Genome_branches_mask$PanHomoGorPon)
panTro4_rb <- Genome_branches_mask$panTro4 / (Genome_branches_mask$rheMac3 + Genome_branches_mask$ponAbe2 + Genome_branches_mask$gorGor3 + Genome_branches_mask$panTro4 + Genome_branches_mask$hg19 + Genome_branches_mask$PanHomo+ Genome_branches_mask$PanHomoGor + Genome_branches_mask$PanHomoGorPon)
hg19_rb    <- Genome_branches_mask$hg19    / (Genome_branches_mask$rheMac3 + Genome_branches_mask$ponAbe2 + Genome_branches_mask$gorGor3 + Genome_branches_mask$panTro4 + Genome_branches_mask$hg19 + Genome_branches_mask$PanHomo+ Genome_branches_mask$PanHomoGor + Genome_branches_mask$PanHomoGorPon)

# I include the inner branches, just in case

PanHomo_rb <- Genome_branches_mask$PanHomo   / (Genome_branches_mask$rheMac3 + Genome_branches_mask$ponAbe2 + Genome_branches_mask$gorGor3 + Genome_branches_mask$panTro4 + Genome_branches_mask$hg19 + Genome_branches_mask$PanHomo+ Genome_branches_mask$PanHomoGor + Genome_branches_mask$PanHomoGorPon)
PanHomoGor_rb <- Genome_branches_mask$PanHomoGor   / (Genome_branches_mask$rheMac3 + Genome_branches_mask$ponAbe2 + Genome_branches_mask$gorGor3 + Genome_branches_mask$panTro4 + Genome_branches_mask$hg19 + Genome_branches_mask$PanHomo+ Genome_branches_mask$PanHomoGor + Genome_branches_mask$PanHomoGorPon)
PanHomoGorPon_rb <- Genome_branches_mask$PanHomoGorPon   / (Genome_branches_mask$rheMac3 + Genome_branches_mask$ponAbe2 + Genome_branches_mask$gorGor3 + Genome_branches_mask$panTro4 + Genome_branches_mask$hg19 + Genome_branches_mask$PanHomo+ Genome_branches_mask$PanHomoGor + Genome_branches_mask$PanHomoGorPon)


# make a new dataframe of relative branch lenghts:

relBraches_mask <- as.data.frame(cbind(panTro4_rb, hg19_rb, ponAbe2_rb, gorGor3_rb, rheMac3_rb, PanHomo_rb, PanHomoGor_rb, PanHomoGorPon_rb))

# combine with the original table

branch_full_mask <- cbind(Genome_branches_mask, relBraches_mask)


library(ggplot2)

ggplot(branch_full_mask) +
  geom_density(data =branch_full_mask, aes(x=hg19_rb),  colour="black", fill="darkblue") +
  geom_density(data = branch_full_mask, aes(x=panTro4_rb),  colour="black", fill="green", alpha= 0.3) +
  geom_density(data =branch_full_mask, aes(x=gorGor3_rb),  colour="black", fill="blue", alpha= 0.5) +
  geom_density(data =branch_full_mask, aes(x=ponAbe2_rb),  colour="black", fill="red", alpha= 0.5) +
  geom_density(data = branch_full_mask, aes(x=rheMac3_rb),  colour="black", fill="yellow", alpha= 0.5) +
  theme_bw() + labs(x = "Branch length", y="Density", title = "Reference Branch Length Distribution (Masked)") +
  scale_x_continuous(limits = c(0, 0.5))
  
  
```

![alt text](https://github.com/wodanaz/adaptiPhy/blob/master/fig_rel_branch.png?raw=true)


```R

# remove trees with very small substitution rate. Change the value of the parameter (0.001) according to the distribution figure. The idea is to remove the peak of sites with a substitution rate near to zero (these are highly conserved sites).

NoMissData_mask <- subset(subset(subset(subset(subset(branch_full_mask, hg19>0.001), panTro4 > 0.001), gorGor3 > 0.001), ponAbe2 > 0.001), rheMac3 > 0.001)


library(ggplot2)

ggplot(NoMissData_mask) +
  geom_density(data =NoMissData_mask, aes(x=hg19_rb),  colour="black", fill="darkblue") +
  geom_density(data = NoMissData_mask, aes(x=panTro4_rb),  colour="black", fill="green", alpha= 0.3) +
  geom_density(data =NoMissData_mask, aes(x=gorGor3_rb),  colour="black", fill="blue", alpha= 0.5) +
  geom_density(data =NoMissData_mask, aes(x=ponAbe2_rb),  colour="black", fill="red", alpha= 0.5) +
  geom_density(data = NoMissData_mask, aes(x=rheMac3_rb),  colour="black", fill="yellow", alpha= 0.5) +
  theme_bw() + labs(x = "Branch length", y="Density", title = "Reference Branch Length Distribution (Masked)") +
  scale_x_continuous(limits = c(0, 0.5))

```

![alt text](https://github.com/wodanaz/adaptiPhy/blob/master/fig_rel_branch_2.png?raw=true)

```R
# compute mean, median and first and third quantile for the human branch:

mean_hg19 <- mean(NoMissData_mask$hg19_rb)
median_hg19 <- median(NoMissData_mask$hg19_rb)

quantile(NoMissData_mask$hg19_rb)
q1 <- matrix(summary(NoMissData_mask$hg19_rb))[2,]
q3 <- matrix(summary(NoMissData_mask$hg19_rb))[5,]


neutralset <- subset(NoMissData_mask, hg19_rb>q1 & hg19_rb<q3)
dim(neutralset)

# Finally, save this set of sequences. We assume these are non-functional and putatively neutral
neutralset.txt <- paste(neutralset$chromosome, "fa", sep=".")


write.table(neutralset.txt, file ="neutralset.txt", row.names=F, col.names=F, quote=F) 


```


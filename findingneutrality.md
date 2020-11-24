

# To create a subset list of putative neutral elements:

First, we need to compute the substitution rate among each node and tip of the three. Here, we used PhyloFit to do this.

```bash

mkdir -p MODELS_HKY85
nano domodel.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
for file in *.fa ; 
do root=`basename $file .fa`; 
phyloFit $file --tree "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))" --subst-mod HKY85 --out-root MODELS_HKY85/$root; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x



cd MODELS_HKY85

#For Reference
for filename in *.mod; do grep -H "TREE:" $filename; done > output.hky85.txt
cat output.hky85.txt | awk -F":" '{print $1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11  }' | \
awk -F"," '{print $1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10  }'  | \
awk -F")" '{print $1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10  }'  | \
awk '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $13 }'  > BranchLenghts.tab
# Let's stop here for a while.... we need to check the table and make sure it's looking good

# we need to further modify column 1 alone -> I don't like the dots, it would be nicer to have chromosome and location in separate columns

awk -F"." '{print $1 ":" $2  }' BranchLenghts.tab > chr_pos.tab

paste chr_pos.tab BranchLenghts.tab  | column -s '\t' -t | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }'  > Branches.tab

sed 1i"chromosome\trheMac3\tponAbe2\tgorGor3\tpanTro4\thg19\tPanHomo\tPanHomoGor\tPanHomoGorPon" Branches.tab > Branches.mask.data


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
dim(moreneutralset)

# Finally, save this set of sequences. We assume these are non-functional and putatively neutral


write.table(neutralset.txt, file ="neutralset.txt", row.names=F, col.names=F, quote=F) 


```


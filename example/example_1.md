
# A markdown pipeline of my analysis of selection using SARS-CoV-2 as a reference species


First, I downloaded files from email

Query_Align.tar.gz
REF.tar.gz


Uploaded data to Hardac


Tar in new directory query

```bash
tar -xzvf Query_Align.tar.gz
```

Tar in new directory ref


```bash
tar -xzvf REF.tar.gz
```
 
 
 
 Remove giberish and empty alignments from files and make a new queries list file
 
 ```bash
 ls -1 | awk -F"." '{ print $1 "." $2 }' > queries.list 
 ```
 
 
```bash

nano bfgenerator_global.py
import sys
import csv
import random
query=sys.argv[1]


with open(query) as f:
    querylist = f.read().splitlines() 







model = ['null','alt']
branches = ['hg38','panPan2','panTro5']
for lrt in model:
	for ape in branches:
		for i in querylist:
			k = random.randint(1,1000);
			f = open('%s.%s.%s.bf' % (i,ape,lrt), 'w')
			f.write('random_seed=%i;\n' % k)
			f.write('quer_seq_file= "/data/wraycompute/alejo/PS_tests/primate2/aish/query/%s.fa.prunned";\n' % i)
			f.write('ref_seq_file = "/data/wraycompute/alejo/PS_tests/primate2/aish/ref/%s.ref.fa";\n' % i)
			f.write('fit_repl_count = 20;\n')
			f.write('tree= "(rheMac8,(hg38,(panTro5,panPan2)))";\n')
			f.write('fgrnd_branch_name = "%s";\n' % ape)
			f.write('res_file = "res/%s.%s.%s.res";\n' % (i,ape,lrt))
			f.write('#include "%s4-fgrnd_spec.bf";\n' % lrt)





#exit nano ctrl+O ENTER ctrl+x
```


```bash
nano dobf.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=xxxxxx@xxxxxx.com
#SBATCH --mem=150
module load Anaconda/1.9.2-fasrc01
python bfgenerator_global.py queries.list
```

Submit to HPC (Slurm)

```bash

sbatch dobf.sh
```

Make lists of bfs for each branch and model

```bash
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
	echo $chr.*hg38.null.bf >> $chr.null.hg38.list;
	echo $chr.*hg38.alt.bf >> $chr.alt.hg38.list;
	echo $chr.*panTro5.null.bf >> $chr.null.panTro5.list;
	echo $chr.*panTro5.alt.bf >> $chr.alt.panTro5.list;	
	echo $chr.*panPan2.null.bf >> $chr.null.panPan2.list;
	echo $chr.*panPan2.alt.bf >> $chr.alt.panPan2.list;
done

```

```bash

rm shgenerator.py
nano shgenerator.py
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=xxxxxx@xxxxxx.com
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem-per-cpu=2000
import sys
import csv
import random

querylist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10','chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
model = ['null','alt']
branches = ['hg38','panPan2','panTro5']
for lrt in model:
	for ape in branches:
		for i in querylist:
			f = open('%s.%s.%s.sh' % (i,ape,lrt), 'w')
			f.write('#!/usr/bin/env bash\n')
			f.write('for file in `cat %s.%s.%s.list`; do root=`basename $file .%s.%s.bf`; HYPHYMP $file > HYPHY/$root.%s.%s.out;done\n' % (i,lrt,ape,ape,lrt,ape,lrt))









#exit nano ctrl+O ENTER ctrl+x

```

Run python

```bash

module load Anaconda/1.9.2-fasrc01
python shgenerator.py


```


```bash

mkdir HYPHY
mkdir res

module load hyphy
for file in chr*hg38*sh ; do sbatch $file ; done
for file in chr*panTro5*sh ; do sbatch $file ; done
for file in chr*panPan2*sh ; do sbatch $file ; done



```



```bash

# and for the query files ( do for both HKY85 and GTR substitution models )



cd /data/wraycompute/alejo/PS_tests/primate2/aish/query
mkdir -p MODELS_HKY85
nano domodel.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=xxxxxx@xxxxxx.com
for file in `cat query.list` ; 
do phyloFit $file.fa.prunned --tree "(rheMac8,(hg38,(panTro5,panPan2)))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85/$file; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x



cd ..



cd /data/wraycompute/alejo/PS_tests/primate2/aish/ref
mkdir -p MODELS_HKY85
rm domodel.sh
nano domodel.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=xxxxxx@xxxxxx.com
for file in `cat ref.list` ; 
do phyloFit $file.ref.fa --tree "(rheMac8,(hg38,(panTro5,panPan2)))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85/$file; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x

```



```bash

#For Reference



for filename in *mod; do grep -H "TREE:" $filename; done > output.hky85.txt

cat output.hky85.txt | sed -r 's/\(+/ /g' |  sed -r 's/\)/ /g'  | sed -r 's/:/ /g' |  sed -r 's/,/ /g' | sed -r 's/;//g' | awk '{ print $1  "\t" $4  "\t" $6  "\t" $8 "\t" $10  }' > BranchLenghts.tab

awk -F"." '{print $1 "\t" $2 }' BranchLenghts.tab  > chr_pos.tab

paste chr_pos.tab BranchLenghts.tab  | column -s '\t' -t | awk '{print $1  ":" $2 "\t"  $4 "\t" $5 "\t" $6 "\t" $7 }'  > Branches.tab

sed 1i"genome_location\trheMac8\thg38\tpanTro5\tpanPan2" Branches.tab > R.hky85.Branches.tab





#For query


for filename in *mod; do grep -H "TREE:" $filename; done > output.hky85.txt

cat output.hky85.txt | sed -r 's/\(+/ /g' |  sed -r 's/\)/ /g'  | sed -r 's/:/ /g' |  sed -r 's/,/ /g' | sed -r 's/;//g' | awk '{ print $1  "\t" $4  "\t" $6  "\t" $8 "\t" $10  }' > BranchLenghts.tab

awk -F"." '{print $1 "\t" $2 }' BranchLenghts.tab  > chr_pos.tab

paste chr_pos.tab BranchLenghts.tab  | column -s '\t' -t | awk '{print $1  ":" $2 "\t"  $4 "\t" $5 "\t" $6 "\t" $7 }'  > Branches.tab

sed 1i"genome_location\trheMac8\thg38\tpanTro5\tpanPan2" Branches.tab > Q.hky85.Branches.tab






```





```bash 

cd ..
cd ..


cp query/MODELS_HKY85/Q.hky85.Branches.tab .
cp ref/MODELS_HKY85/R.hky85.Branches.tab .


wc -l Q.hky85.Branches.tab
wc -l R.hky85.Branches.tab



sed -i -e "1d" Q.hky85.Branches.tab
sed -i -e "1d" R.hky85.Branches.tab



sort Q.hky85.Branches.tab -k1,1 -V > Q.hky85.Branches.sorted.tab
sort R.hky85.Branches.tab -k1,1 -V > R.hky85.Branches.sorted.tab


# genome_location	rheMac8	hg38	panTro5	panPan2



join -t $'\t' -j 1 -o  1.1 1.2 2.2 1.3 2.3 1.4 2.4 1.5 2.5  Q.hky85.Branches.sorted.tab R.hky85.Branches.sorted.tab   > join_trans.tab

awk '{  rate1 = $4 / $5 ;   print $1   "\t"  rate1  "\t" "hg38"  }'             join_trans.tab  > 	PhyloFit.hg38.tab
awk '{  rate2 = $6 / $7 ;   print $1   "\t"  rate2  "\t" "panTro5" }'          join_trans.tab  > 	PhyloFit.panTro5.tab
awk '{  rate3 = $8 / $9 ;   print $1   "\t"  rate3  "\t" "panPan2" }'        join_trans.tab  > 	PhyloFit.panPan2.tab


cat PhyloFit.hg38.tab  PhyloFit.panTro5.tab  PhyloFit.panPan2.tab  |  sort -k1 -k2,2n  -V | sed 1i"genome_location\tzeta\tspecies" > aisha.phyloFit.data



```





### After HyPhy had finished running. We can retrieve data from the res/ directory

```bash 
cd res
```
# we need to model the likelihood ratio test results in a chi-square distribution to obtain p-values for each pair of null v alternative models...
# to do this, we can use the program phase8.rb .... but maybe I can implement it in R



```bash
for file in *hg38.null.res; do echo $file >>null.hg38.log; done;
for file in *hg38.alt.res; do echo $file >>alt.hg38.log; done;

for file in *panTro5.null.res; do echo $file >>null.panTro5.log; done;
for file in *panTro5.alt.res; do echo $file >>alt.panTro5.log; done;

for file in *panPan2.null.res; do echo $file >>null.panPan2.log; done;
for file in *panPan2.alt.res; do echo $file >>alt.panPan2.log; done;

```


```bash

wc -l *log


for filename in `cat alt.hg38.log`; do grep -H "BEST LOG-L:" $filename >> alt.hg38.tab; done;
for filename in `cat null.hg38.log`; do grep -H "BEST LOG-L:" $filename >> null.hg38.tab; done;

for filename in `cat alt.panTro5.log`; do grep -H "BEST LOG-L:" $filename >> alt.panTro5.tab; done;
for filename in `cat null.panTro5.log`; do grep -H "BEST LOG-L:" $filename >> null.panTro5.tab; done;

for filename in `cat alt.panPan2.log`; do grep -H "BEST LOG-L:" $filename >> alt.panPan2.tab; done;
for filename in `cat null.panPan2.log`; do grep -H "BEST LOG-L:" $filename >> null.panPan2.tab; done;
```


# Consolidating tables for three branches

```bash
awk '{print $1 "\t" $3}' null.hg38.tab | sort -k1,1 -V  > nulls.hg38.tab
awk '{print $1 "\t" $3}' alt.hg38.tab | sort -k1,1 -V > alts.hg38.tab

awk '{print $1 "\t" $3}' null.panTro5.tab | sort -k1,1 -V  > nulls.panTro5.tab
awk '{print $1 "\t" $3}' alt.panTro5.tab | sort -k1,1 -V > alts.panTro5.tab

awk '{print $1 "\t" $3}' null.panPan2.tab | sort -k1,1 -V  > nulls.panPan2.tab
awk '{print $1 "\t" $3}' alt.panPan2.tab | sort -k1,1 -V > alts.panPan2.tab




awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.hg38.tab > col1.hg38.tab
paste col1.hg38.tab nulls.hg38.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $7   }' > lnulls.hg38.tab
paste lnulls.hg38.tab alts.hg38.tab | awk '{print $1 ":" $2  "\t" $4 "\t"  $6 "\t" "hg38" }'  |  sort -k1 -V  > likelihoods.hg38.tab
 
awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.panTro5.tab > col1.panTro5.tab
paste col1.panTro5.tab nulls.panTro5.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $7   }' > lnulls.panTro5.tab
paste lnulls.panTro5.tab alts.panTro5.tab | awk '{print $1 ":" $2  "\t" $4 "\t"  $6 "\t" "panTro5" }' |  sort -k1 -V  > likelihoods.panTro5.tab
 
awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.panPan2.tab > col1.panPan2.tab
paste col1.panPan2.tab nulls.panPan2.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $7   }' > lnulls.panPan2.tab
paste lnulls.panPan2.tab alts.panPan2.tab | awk '{print $1 ":" $2  "\t" $4 "\t"  $6 "\t" "panPan2" }' |  sort -k1 -V  > likelihoods.panPan2.tab



 
cat likelihoods.hg38.tab  likelihoods.panTro5.tab  likelihoods.panPan2.tab   |  sort -k1 -k2,2n -k3,3n | sed 1i"location\tlnull\tlalt\tspecies" > likelihoods.tab


```


 
# to go to R:


```bash
module load R
R

likelihoods = read.table("likelihoods.tab", header = TRUE)  # read tab file
LRT <- -2*(likelihoods$lnull - likelihoods$lalt)
pval <- 1-pchisq(LRT, 1)
l_pvals <- cbind(likelihoods, pval, -log(pval))
write.table(l_pvals, file ="likelihoods.pvals.tab", row.names=F, col.names=F, quote=F) 


# lets's exit the R environment
q()

```


###################################################
# 3 columns


```bash
awk '{ print $1 "\t" $5  "\t" $4 }' likelihoods.pvals.tab | sort -k1,1 -V  >  likelihoods.pvals2.tab
sed 1i"genome_location\tpval\tspecies" likelihoods.pvals2.tab > aisha.adaptiphy.data


cd ..
cd ..
cp test/res/aisha.adaptiphy.data .

sed -i -e "1d" aisha.adaptiphy.data

```






```bash



ls *data


#aisha.adaptiphy.data  aisha.phyloFit.data
module load R
R


```R
aisha_hyphy = as.data.frame(read.table("aisha.adaptiphy.data", header = F)) # read tab file 
colnames(aisha_hyphy) <- c("genome_location",  "pval",          "species")
aisha_zeta = as.data.frame(read.table("aisha.phyloFit.data", header = T)) # read tab file 
head(aisha_hyphy)
head(aisha_zeta)


aisha.selection.data <- merge(aisha_zeta, aisha_hyphy, by= c('genome_location', 'species'))
head(aisha.selection.data )


aisha.selection.data$QsubsRate <- NULL
aisha.selection.data$RsubsRate <- NULL


library(reshape)
aisha.selection.wide.data <- reshape(aisha.selection.data , idvar = c("genome_location"), timevar = "species", direction = "wide")
head(aisha.selection.wide.data)



write.table(aisha.selection.wide.data , file ="aisha.selection.data", row.names=F, col.names=T, quote=F) 


q()

```





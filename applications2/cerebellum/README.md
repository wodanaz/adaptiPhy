# Analysing evolution by positive selection in a single dataset of cerebellum atac seq


```bash
mkdir features
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY ; 
	do grep -w $chr cerebellum_atac.bed  | awk '{print $1 "\t" $2 - 1 "\t" $3 }' | sort -k1,1 -k2,2 -V >  features/$chr.feat.bed; 
done

```


To generate alignments for Cerebellum regions

```bash
mkdir query    
cd /data/wraycompute/alejo/PS_tests/primate/unmasked1
# for minus queries
rm do_queries.sh
nano do_queries.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=2000
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;
do msa_split $chr.primate.maf --refseq $chr.fa --gap-strip ANY -q --in-format MAF --features /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/features/$chr.feat.bed --for-features --out-root /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/query/$chr; done
```
Run it

```bash

sbatch do_queries.sh

```

Go to query directory and make a list of alignments


```bash

cd /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/query


for file in *fa; do echo $file >> all.list;done


wc -l all.list

head all.list
```

Need to clean alignments and prun deletions relative to the human genome

```bash
nano prun.py.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH --mem=200
module load Anaconda/1.9.2-fasrc01
python prunning.py 
```

Run it


```bash
sbatch prun.py.sh
```

Wait until it finishes running


```bash
for file in *prunned; do echo $file >> all.prunned.list;done
```

Filter out any alignment shorter than 250bps

```bash
rm filtering.py
nano filtering.py
import re # for Regular expressions
import sys
from Bio import AlignIO

#infile=sys.argv[1]

infile = 'all.prunned.list'

with open(infile) as f:
    reflist = f.read().splitlines() 



f.close()
missing = '[(+*)]'
ambiguous = '[N]'


badsimiosNs = []
badsimiosAsk = []
goodsimios = []

for i in reflist:
	myfile= open(i, 'r')
	mydata = myfile.read()
	fasta = AlignIO.read(i, "fasta")
	myfile.close()
	if fasta.get_alignment_length() > 250 :
		maxNs = len(re.findall(ambiguous, str(mydata), re.I))
		maxAsk = len(re.findall(missing, str(mydata), re.I))
		if maxNs > 2:
			print ')-B'
			badsimiosNs.append(i)
		elif maxAsk > 1:
			print ')-;'
			badsimiosAsk.append(i)
		else:
			print 'lol'
			goodsimios.append(i)
		




thefile3 = open('goodalignments.txt', 'w')
for item in goodsimios:
  thefile3.write("%s\n" % item)





thefile3.close()





thefile1 = open('ambiguous.txt', 'w')
for item in badsimiosNs:
  thefile1.write("%s\n" % item)




thefile1.close()

thefile2 = open('asterisks.txt', 'w')
for item in badsimiosAsk:
  thefile2.write("%s\n" % item)





thefile2.close()

```


Run this python script



```bash
nano dofilter.py.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH --mem=200
module load Anaconda/1.9.2-fasrc01
python filtering.py 

```

```bash
sbatch dofilter.py.sh
```
Check list lengths


```bash
  
wc -l ambiguous.txt asterisks.txt goodalignments.txt
```

Now, that the queries were created, we need to make a concatenated reference for each query using non-functional regions of NFRs.

```bash

cd ..
cat query/goodalignments.txt | awk -F"." '{ print $1 "\t" $2}' | awk -F"-" '{ print $1 "\t" $2 "\t" $3}' | sort -k1,1 -k2,2 -V > queries.bed  #100K
cat query/goodalignments.txt > queries_cerebelum.list


cp queries_cerebelum.list /data/wraycompute/alejo/PS_tests/primate/maskedv3/refmasked2

cd /data/wraycompute/alejo/PS_tests/primate/maskedv3/refmasked2

```

Now let's run a pyhon script that creates a dictionary using the list I just created.

```bash
module load Anaconda/1.9.2-fasrc01
rm DictGen_cerbellum.py
nano DictGen_cerbellum.py
import re # for Regular expressions
from Bio import AlignIO
import sys
import csv
import random
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment




keys= "queries_cerebelum.list"
values= "neutral.list"


with open(keys) as f:
    keylist = f.read().splitlines() 



with open(values) as f:
    valuelist = f.read().splitlines() 







dictionary = {}
for key in keylist:
	random.shuffle(valuelist)
	dictionary[key] = valuelist[0:10]
from pprint import pprint
fielddict_file = open("global.dict", "w")
pprint(dictionary, fielddict_file)
fielddict_file.close()


	
reference = []


for i,j in dictionary.iteritems():
	n= 0
	combined_seq = MultipleSeqAlignment([SeqRecord(Seq('', generic_dna),id="hg19"), 
						SeqRecord(Seq('', generic_dna),id="panTro4"),
						SeqRecord(Seq('', generic_dna),id="gorGor3"),
						SeqRecord(Seq('', generic_dna),id="rheMac3"),
						SeqRecord(Seq('', generic_dna),id="ponAbe2")])
	combined_seq.sort()
	for ref in j:
		n = n + 1
		seq_records = AlignIO.read(ref, 'fasta')
		seq_records.description = ""
		seq_records.sort()
		combined_seq =  combined_seq + seq_records
		combined_seq.description = ""
	with open('%s.ref' % i, 'w') as write_file:
		AlignIO.write(combined_seq, write_file, 'fasta')
	referencelist = open('reference.list', 'a')	
	referencelist.write('%s\t%i\n' % (i,n))





```


and run it:

```bash
rm dodicts.sh
nano dodicts.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem=150
module load Anaconda/1.9.2-fasrc01
python DictGen_cerbellum.py
```
```bash
sbatch dodicts.sh
```




Wait until it finishes running. And then, mv the alignments to a new ref directory in your working directory

```bash
cd /data/wraycompute/alejo/PS_tests/primate/maskedv3/refmasked2

for file in `cat queries_cerebelum.list`; do echo $file.ref >> ref.cerebellum.tab; done

mkdir /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/ref 

for file in `cat queries_cerebelum.list`; do mv $file.ref /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/ref ; done
```

After moving files, the description in the new fasta files must be changed

```bash
cd  /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/ref 

for file in *prunned.ref ; do echo $file >> all.ref.tab; done


rm do_refs.sh
nano do_refs.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem=1500
for file in `cat all.ref.tab`; do root=`basename $file .fa.prunned.ref`; awk '{if($1 ~ /^>/){split($1,a,"\t"); print a[1]}else{print}}' $file > $root.ref.fa; done
```
```bash
sbatch do_refs.sh
```


Create a new list and remove \*ref files to release some disk space


```bash
awk -F"." '{print $1 "." $2 }' all.ref.tab > ref.list



nano do_remove.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem=1500
for file in `cat all.ref.tab`; do rm $file; done
```
```bash
sbatch do_remove.sh
```


Now, this is almost ready to run AdaptiPhy. 

Copy the scripts where the models are defined to your test direcotry

alt4-fgrnd_spec.bf  null4-fgrnd_spec.bf


```bash
cd ..

mkdir test

cp ref/ref.list test
cp alt4-fgrnd_spec.bf  test
cp null4-fgrnd_spec.bf test
```

The following lines are used to create the batch file creator. This script, generates a batch file for each query for each adaptiPhy run.

```bash
rm bfgenerator_global.py
nano bfgenerator_global.py
import sys
import csv
import random
query=sys.argv[1]


with open(query) as f:
    querylist = f.read().splitlines() 







model = ['null','alt']
branches = ['hg19','panTro4' ]
for lrt in model:
	for ape in branches:
		for i in querylist:
			k = random.randint(1,1000);
			f = open('%s.%s.%s.bf' % (i,ape,lrt), 'w')
			f.write('random_seed=%i;\n' % k)
			f.write('quer_seq_file= "/data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/query/%s.fa.prunned";\n' % i)
			f.write('ref_seq_file = "/data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/ref/%s.ref.fa";\n' % i)
			f.write('fit_repl_count = 20;\n')
			f.write('tree= "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))";\n')
			f.write('fgrnd_branch_name = "%s";\n' % ape)
			f.write('res_file = "res/%s.%s.%s.res";\n' % (i,ape,lrt))
			f.write('#include "%s4-fgrnd_spec.bf";\n' % lrt)





#exit nano ctrl+O ENTER ctrl+x
```

```bash
nano dobf.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem=150
module load Anaconda/1.9.2-fasrc01
python bfgenerator_global.py ref.list
```
```bash
sbatch dobf.sh
```

now, we make a list of each batch file launchers per chromosome


```bash
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
	echo $chr.*hg19.null.bf >> $chr.null.hg19.list;
	echo $chr.*hg19.alt.bf >> $chr.alt.hg19.list;
	echo $chr.*panTro4.null.bf >> $chr.null.panTro4.list;
	echo $chr.*panTro4.alt.bf >> $chr.alt.panTro4.list;	
done
```

Next, we use a python script that generates batsh files for running adaptiphy


```bash
rm shgenerator.py
nano shgenerator.py
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem-per-cpu=2000
import sys
import csv
import random

querylist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10','chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
model = ['null','alt']
branches = [ 'hg19','panTro4' ]
for lrt in model:
	for ape in branches:
		for i in querylist:
			f = open('%s.%s.%s.sh' % (i,ape,lrt), 'w')
			f.write('#!/usr/bin/env bash\n')
			f.write('for file in `cat %s.%s.%s.list`; do root=`basename $file .%s.%s.bf`; HYPHYMP $file > HYPHY/$root.%s.%s.out;done\n' % (i,lrt,ape,ape,lrt,ape,lrt))









#exit nano ctrl+O ENTER ctrl+x
```


Run it:


```bash
module load Anaconda/1.9.2-fasrc01
python shgenerator.py
```



### Run Adaptiphy


```



mkdir HYPHY
mkdir res

module load hyphy
for file in chr*hg19*sh ; do sbatch $file ; done
for file in chr*panTro4*sh ; do sbatch $file ; done

```

### Now we are ready to run PhyloP to get the rates of substituion so, we can compute zeta, the value representing evolutionary ratio


For the queries:
```bash
cd /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/query
mkdir -p MODELS_HKY85
nano domodel.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
for file in `cat goodalignments.txt` ; 
do root=`basename $file .fa.prunned`; 
phyloFit $file --tree "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85/$root; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x

cd ..
```

For the references:

```bash
cd /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/ref
mkdir -p MODELS_HKY85
rm domodel.sh
nano domodel.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
for file in `cat ref.list` ; 
do phyloFit $file.ref.fa --tree "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85/$file; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x
```

Adaptiphy can take some time to finish, in the mean time we can process PhyloFit data.


For the query:

```bash
cd ../query/MODELS_HKY85/

for filename in *mod; do grep -H "TREE:" $filename; done > output.hky85.txt

cat output.hky85.txt | sed -r 's/\(+/ /g' |  sed -r 's/\)/ /g'  | sed -r 's/:/ /g' |  sed -r 's/,/ /g' | sed -r 's/;//g' | awk '{ print $1  "\t" $12  "\t" $10  "\t" $13 "\t" $8  "\t" $14 "\t" $6 "\t" $14 "\t" $4 }' | sed 's/.mod//g' > BranchLenghts.tab

sed 1i"genome_location\thg19\tpanTro4\tHCh\tgorGor3\tHchG\tponAbe2\tHChGO\trheMac3" BranchLenghts.tab > Q.hky85.Branches.tab


```

And, for the references:

```bash
cd ..
cd /ref/MODELS_HKY85/

for filename in *mod; do grep -H "TREE:" $filename; done > output.hky85.txt

cat output.hky85.txt | sed -r 's/\(+/ /g' |  sed -r 's/\)/ /g'  | sed -r 's/:/ /g' |  sed -r 's/,/ /g' | sed -r 's/;//g' | awk '{ print $1  "\t" $12  "\t" $10  "\t" $13 "\t" $8  "\t" $14 "\t" $6 "\t" $14 "\t" $4 }' | sed 's/.mod//g' > BranchLenghts.tab

sed 1i"genome_location\thg19\tpanTro4\tHCh\tgorGor3\tHchG\tponAbe2\tHChGO\trheMac3" BranchLenghts.tab > R.hky85.Branches.tab

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

awk '{print $1}' Q.hky85.Branches.sorted.tab | awk -F":" '{print $1 "\t" $2 }' | awk -F"-" '{print $1 "\t" $2 "\t" $3}' > Q.hky85.bed
awk '{print $1}' R.hky85.Branches.sorted.tab | awk -F":" '{print $1 "\t" $2 }'  | awk -F"-" '{print $1 "\t" $2 "\t" $3}' > R.hky85.bed


paste Q.hky85.bed Q.hky85.Branches.sorted.tab | column -s '\t' -t | awk '{print $1 "\t" $2 "\t" $3 "\t" $9 }' > Q.hky85.hg19.bed # human
paste R.hky85.bed R.hky85.Branches.sorted.tab | column -s '\t' -t | awk '{print $1 "\t" $2 "\t" $3 "\t" $9 }' > R.hky85.hg19.bed # human

paste Q.hky85.bed Q.hky85.Branches.sorted.tab | column -s '\t' -t | awk '{print $1 "\t" $2 "\t" $3 "\t" $8 }' > Q.hky85.panTro4.bed # chimp
paste R.hky85.bed R.hky85.Branches.sorted.tab | column -s '\t' -t | awk '{print $1 "\t" $2 "\t" $3 "\t" $8 }' > R.hky85.panTro4.bed # chimp
```

Compute Zeta

```bash
module load R
R
```

```R
Cerebellum_Q = as.data.frame(read.table("Q.hky85.Branches.sorted.tab", header = F)) # read tab file 
#colnames(K562_data) <- c('chromosome', 'pval.human', 'zeta.human', 'pval.chimp', 'zeta.chimp','phastcons')
Cerebellum_R= as.data.frame(read.table("R.hky85.Branches.sorted.tab", header = F)) # read tab file 
colnames(Cerebellum_Q) 
colnames(Cerebellum_R) 

Cerebellum.zeta <- merge(Cerebellum_Q, Cerebellum_R, by= "V1")

rate1 <- Cerebellum.zeta$V2.x / Cerebellum.zeta$V2.y # human
rate2 <- Cerebellum.zeta$V3.x / Cerebellum.zeta$V3.y # chimp


Cerebellum_selection <- data.frame(Cerebellum.zeta$V1 , rate1, rate2)
colnames(Cerebellum_selection) <- c('chromosome', 'zeta.human', 'zeta.chimp')


write.table(Cerebellum_selection, file ="PhyloFit.Cerebellum.data", row.names=F, col.names=T, quote=F) 


```
```bash
sed -ri 's/\./:/' PhyloFit.Cerebellum.data
```

If the part run by HyPhy has been completed, its time to consolidate tables and compute P-value


```bash


cd test/res

for file in *hg19.null.res; do echo $file >> null.hg19.log; done;
for file in *hg19.alt.res; do echo $file >>alt.hg19.log; done;
for file in *panTro4.null.res; do echo $file >>null.panTro4.log; done;
for file in *panTro4.alt.res; do echo $file >> alt.panTro4.log; done;



wc -l *log


nano domodelX1.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
for filename in `cat alt.panTro4.log`; do grep -H -m 1 "BEST LOG-L:" $filename >> alt.panTro4.tab; done;
for filename in `cat null.panTro4.log`; do grep -H -m 1 "BEST LOG-L:" $filename >> null.panTro4.tab; done;
for filename in `cat alt.hg19.log`; do grep -H -m 1 "BEST LOG-L:" $filename >> alt.hg19.tab; done;
for filename in `cat null.hg19.log`; do grep -H -m 1 "BEST LOG-L:" $filename >> null.hg19.tab; done;
```
```bash
sbatch domodelX1.sh
```

Consolidating tables for 5 branches:
```bash

awk '{print $1 "\t" $3}' null.hg19.tab | sort -k1,1 -V  | sed -r 's/.hg19.null.res:BEST//g' | sed -r 's/\./:/' > nulls.hg19.tab
awk '{print $1 "\t" $3}' alt.hg19.tab | sort -k1,1 -V  | sed -r 's/.hg19.alt.res:BEST//g' | sed -r 's/\./:/'> alts.hg19.tab

awk '{print $1 "\t" $3}' null.panTro4.tab | sort -k1,1 -V | sed -r 's/.panTro4.null.res:BEST//g' | sed -r 's/\./:/'  > nulls.panTro4.tab
awk '{print $1 "\t" $3}' alt.panTro4.tab | sort -k1,1 -V | sed -r 's/.panTro4.alt.res:BEST//g' | sed -r 's/\./:/' > alts.panTro4.tab



```

```bash
module load R
R
```

```R
Cerebellum_nulls_hg19 = as.data.frame(read.table("nulls.hg19.tab", header = F)) # read tab file 
Cerebellum_alts_hg19 = as.data.frame(read.table("alts.hg19.tab", header = F)) # read tab file 
Cerebellum_nulls_panTro4 = as.data.frame(read.table("nulls.panTro4.tab", header = F)) # read tab file 
Cerebellum_alts_panTro4 = as.data.frame(read.table("alts.panTro4.tab", header = F)) # read tab file 



colnames(Cerebellum_nulls_hg19) 
colnames(Cerebellum_alts_hg19) 

head(Cerebellum_nulls_hg19) 
head(Cerebellum_alts_hg19) 

Cerebellum.likelihoods.hg19 <- merge(Cerebellum_nulls_hg19, Cerebellum_alts_hg19, by= "V1")
LRT <- -2*(Cerebellum.likelihoods.hg19$V2.x - Cerebellum.likelihoods.hg19$V2.y)
Cerebellum.likelihoods.hg19$pval <- 1-pchisq(LRT, 1)
Cerebellum.likelihoods.hg19$V2.x <-NULL
Cerebellum.likelihoods.hg19$V2.y <-NULL

Cerebellum.likelihoods.panTro4 <- merge(Cerebellum_nulls_panTro4, Cerebellum_alts_panTro4, by= "V1")
LRT <- -2*(Cerebellum.likelihoods.panTro4$V2.x - Cerebellum.likelihoods.panTro4$V2.y)
Cerebellum.likelihoods.panTro4$pval <- 1-pchisq(LRT, 1)
Cerebellum.likelihoods.panTro4$V2.x <-NULL
Cerebellum.likelihoods.panTro4$V2.y <-NULL


Cerebellum.pvals <- merge(Cerebellum.likelihoods.hg19, Cerebellum.likelihoods.panTro4, by= "V1")
colnames(Cerebellum.pvals) <- c('chromosome', 'pval.human', 'pval.chimp')

write.table(Cerebellum.pvals, file ="adaptiphy.pvals.tab", row.names=F, col.names=T, quote=F) 


# lets's exit the R environment
q()
```

Copy resulting file into working directory

```bash
cd ..
cd ..

cp test/res/adaptiphy.pvals.tab .
```



# PhastCons from scratch:
Here a way to determine the degree of conservation among primates.

```bash

mkdir -p TREES
rm -f TREES/*      # in case old versions left over


nano dophastcons_train.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
for file in `cat ref.list` ; do
phastCons -i FASTA --estimate-trees TREES/$file query/$file.fa.prunned ref/MODELS_HKY85/$file.mod --no-post-probs
done;
```
```bash
sbatch dophastcons_train.sh 
```

For predicting:
```bash

mkdir SCORES
mkdir MOSTCONS
nano dophastcons_prediction.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
for file in `cat ref.list` ; do
phastCons -i FASTA query/$file.fa.prunned --most-conserved MOSTCONS/$file.bed TREES/$file.cons.mod,TREES/$file.noncons.mod > SCORES/$file.wig
done;
```
```bash
sbatch dophastcons_prediction.sh
```


For consolidating in a table:
```bash
cd SCORES
nano dophastcons_final.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH -n 1
for file in *wig ; do root=`basename $file .wig`;  sed -i -e "1d" $file  ; awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $file >  $root.average.wig ; done
for filename in *average.wig; do grep -H "" $filename ; done > output.wig.txt
sed -r 's/.average.wig:/\t/g'  output.wig.txt | sed -r '/\./s/\./:/' |  sort -k1 -V > cerebellum.phastCons.data
```
```bash
sbatch dophastcons_final.sh
```

# Obtain GREAT annotations

Upload a bed file of the tested queries (queries.bed) in the form: 

>chr	start	end	chr:start-end






```bash

awk '{print $1 "\t" $2 "\t" $3 "\t" $1 ":" $2 "-" $3 }' queries.bed  > cerebellum.bed2great.bed
```

Upload to GREAT at: http://great.stanford.edu/

And examine the closest single TSS within 10000 kb

Then download the resulting file


```bash
sed -ie "1d" cerebellum_great.txt 
sed -r 's/ \(\+/\t/g' cerebellum_great.txt | sed -r 's/ \(/\t/g'  | sed -r 's/\)/\t/g' | awk '{print $1 "\t" $2 "\t" $3 }' > cerebellum.great.data

grep "NONE" cerebellum.great.data -v > cerebellum.great.nonones.data 


```


Upload it to HARDAC


### Consolidate all results in a single table


```bash
module load R
R
```

```R

adaptiphy = as.data.frame(read.table("AdaptiPhy.Cerebellum.data", header = T)) # read tab file 
head(adaptiphy)
colnames(adaptiphy) <- c("genome_location",  "pval.human",  "pval.chimp" )


PhyloFit = as.data.frame(read.table("PhyloFit.Cerebellum.data", header = T)) # read tab file 
head(PhyloFit)
colnames(PhyloFit) <- c("genome_location",  "zeta.human",  "zeta.chimp" )
head(PhyloFit)

cerebellum_adaptiphy <- merge(PhyloFit, adaptiphy, by= c('genome_location'))
head(cerebellum_adaptiphy)


PhastCons = as.data.frame(read.table("cerebellum.phastCons.data", header = F)) # read tab file 
head(PhastCons)
colnames(PhastCons) <- c("genome_location", "PhastCons"  )
head(PhastCons)
dim(PhastCons)


cerebellum_adaptiphy2 <- merge(cerebellum_adaptiphy, PhastCons, by= c('genome_location'))
head(cerebellum_adaptiphy2)
dim(cerebellum_adaptiphy2)


GREAT = as.data.frame(read.table("cerebellum.great.nonones.data", header = F)) # read tab file 
head(GREAT)
colnames(GREAT) <- c("genome_location", "gene_id" , "dTSS" )


cerebellum_adaptiphy3 <- merge(cerebellum_adaptiphy2, GREAT, by= c('genome_location'), all.x = T)
head(cerebellum_adaptiphy3)


write.table(cerebellum_adaptiphy3 , file ="cerebellum.selection.data", row.names=F, col.names=T, quote=F) 
```





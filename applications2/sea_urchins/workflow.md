
# Workflow for analysis of evolution by positive selection in sea urchins

To run a test of selection, we need a set of query regions (promoters, ATACseq peaks, ChIPseq peaks, etc) and a set of putatively neutral regions.

The neutral regions are premade and they have been generated using the following instructions:

https://github.com/wodanaz/adaptiPhy/blob/master/applications2/sea_urchins/finding_sea_urchin_neutral_elements.md

The list of neutral regions contains about 21K akignments longer that 200 bp and they are saved at: 

```bash


cat allcis.bed alltrans.bed > allpeaks.bed
wc -l allpeaks.bed
40189 allpeaks.bed


mkdir features
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21  ; 
	do grep -w $chr allpeaks.bed  | awk '{print $1 "\t" $2 - 1 "\t" $3 }' | sort -k1,1 -k2,2 -V >  features/$chr.feat.bed; 
done


mkdir query
cd /data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/unmasked_urchins_maf

# for all queries
rm do_queries.sh
nano do_queries.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=2000
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;
do msa_split $chr.He.maf --refseq $chr.He.fa --gap-strip ANY -q --in-format MAF --features /data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/features/$chr.feat.bed --for-features --out-root /data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/query/$chr; done

```

Run it...

```bash
sbatch do_queries.sh
```


```bash
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
python2 prunning.py 
```
Run it

```bash
sbatch prun.py.sh
```

Wait until it finishes running


```bash
for file in *prunned; do echo $file >> all.prunned.list;done
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


badseqsNs = []
badseqsAsk = []
goodseqs = []

for i in reflist:
	myfile= open(i, 'r')
	mydata = myfile.read()
	fasta = AlignIO.read(i, "fasta")
	myfile.close()
	if fasta.get_alignment_length() > 250 : 
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

```bash
cd ..
cat query/goodalignments.txt | awk -F"." '{ print $1 "\t" $2}' | awk -F"-" '{ print $1 "\t" $2 "\t" $3}' | sort -k1,1 -k2,2 -V > queries.bed  #100K
cat query/goodalignments.txt > queries_urchins.list


cp queries_seaurchins.list /data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/neutral_set

cd /data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/neutral_set

```



```bash
nano DictGen.py
import re # for Regular expressions
from Bio import AlignIO
import sys
import csv
import random
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment




keys= "goodalignments.txt"
values= "neutralset.txt"


with open(keys) as f:
    keylist = f.read().splitlines() 



with open(values) as f:
    valuelist = f.read().splitlines() 





for replicate in range(10):

	dictionary = {}
	for key in keylist:
		random.shuffle(valuelist)
		dictionary[key] = valuelist[0:20]
	from pprint import pprint
	fielddict_file = open("global.dict", "w")
	pprint(dictionary, fielddict_file)
	fielddict_file.close()
	reference = []
	for i,j in dictionary.iteritems():
		n= 0
		combined_seq = MultipleSeqAlignment([SeqRecord(Seq('', generic_dna),id="He"), 
							SeqRecord(Seq('', generic_dna),id="Ht"),
							SeqRecord(Seq('', generic_dna),id="Lv")])
		combined_seq.sort()
		for ref in j:
			n = n + 1
			seq_records = AlignIO.read(ref, 'fasta')
			seq_records.description = ""
			seq_records.sort()
			combined_seq =  combined_seq + seq_records
			combined_seq.description = ""
		with open('%s.%i.ref' % (i,replicate), 'w') as write_file:
			AlignIO.write(combined_seq, write_file, 'fasta')
		referencelist = open('reference.list', 'a')	
		referencelist.write('%s\t%i\t%i\n' % (i,n, replicate))


```

```bash
nano genDict.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
module load Anaconda/1.9.2-fasrc01
python DictGen.py

```
Run it:

```bash
sbatch genDict.sh
```


```bash

awk -F"." '{print $1 "." $2 }' neutral_set/goodalignments.txt  > ref.list

nano convert.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
for file in `cat ref.list`; do for i in {0..9} ; do awk '{if($1 ~ /^>/){split($1,a,"\t"); print a[1]}else{print}}' $file.fa.prunned.$i.ref > /data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/ref/$file.$i.ref.fa; done; done


cd ..

mkdir test
cp query/ref.list test


cd test

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
branches = ['He','Ht' ]

for replicate in range(10):
	for lrt in model:
		for urchin in branches:
			for i in querylist:
				k = random.randint(1,1000);
				f = open('%s.%s.%s.%i.bf' % (i,urchin,lrt,replicate), 'w')
				f.write('random_seed=%i;\n' % k)
				f.write('quer_seq_file= "/data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/query/%s.fa.prunned";\n' % i)
				f.write('ref_seq_file = "/data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/ref/%s.%i.ref.fa";\n' % (i,replicate))
				f.write('fit_repl_count = 20;\n')
				f.write('tree= "((He,Ht),Lv)";\n')
				f.write('fgrnd_branch_name = "%s";\n' % urchin)
				f.write('res_file = "res/%s.%s.%s.%i.res";\n' % (i,urchin,lrt,replicate))
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




nano lists.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem=150
for i in {0..9} ; do echo *He.alt.$i.bf >> alt.He.$i.list ; done
for i in {0..9} ; do echo *He.null.$i.bf >> null.He.$i.list ; done
for i in {0..9} ; do echo *Ht.alt.$i.bf >> alt.Ht.$i.list ; done
for i in {0..9} ; do echo *Ht.null.$i.bf >> null.Ht.$i.list ; done
```


The following script will generate launchers for each branch and query


```bash

module load Anaconda3
conda activate aleconda


nano shgenerator.py
import sys
import csv
import random





model = ['null','alt']
branches = ['He','Ht' ]
for replicate in range(10):
	for lct in model:
		for urchin in branches:
			f = open('%s.%s.%i.sh' % (urchin,lct, replicate), 'w')
			f.write('#!/usr/bin/env bash\n')
			f.write('for file in `cat %s.%s.%i.list`; do root=`basename $file .%s.%s.%i.bf`; hyphy $file > HYPHY/$root.%s.%s.%i.out;done\n' % (lct,urchin,replicate,urchin,lct,replicate,urchin,lct,replicate))









#exit nano ctrl+O ENTER ctrl+x
```



```bash
Run and create new directiries to save data

python shgenerator.py 





mkdir HYPHY
mkdir res


for file in Ht*sh ; do sbatch $file ; done
for file in He*sh ; do sbatch $file ; done





```

```bash
cd ..
cd query



rm MODELS_HKY85_query -r
mkdir -p MODELS_HKY85_query
rm domodel_query.sh
nano domodel_query.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=your_email@email.com
for file in `cat ref.list` ; do
phyloFit $file.fa.prunned --tree "(Lv, (Ht,He))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85_query/$file; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x


cd ..

cd ref

rm MODELS_HKY85_ref -r
rm domodel_ref.sh
mkdir -p MODELS_HKY85_ref
nano domodel_ref.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=your_email@email.com
for file in `cat ref.list` ; do for i in {0..9} ; do
phyloFit $file.$i.ref.fa --tree "(Lv, (Ht,He))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85_ref/$file.$i; # HKY85 model, It runs fast and it also the model applied in HYPHY
done; done #exit nano ctrl+O ENTER ctrl+x

```

```bash
sbatch domodel_ref.sh
sbatch domodel_query.sh
```

Retrieving the MLEs and preform Likelihood Ratio Test (LRT)

```bash
for file in *He.null.*.res; do echo $file >>null.He.log; done;
for file in *He.alt.*.res; do echo $file >>alt.He.log; done;

for file in *Ht.null.*.res; do echo $file >>null.Ht.log; done;
for file in *Ht.alt.*.res; do echo $file >>alt.Ht.log; done;



wc -l *log


for filename in `cat alt.He.log`; do grep -H "BEST LOG-L:" $filename >> alt.He.tab; done;
for filename in `cat null.He.log`; do grep -H "BEST LOG-L:" $filename >> null.He.tab; done;

for filename in `cat alt.Ht.log`; do grep -H "BEST LOG-L:" $filename >> alt.Ht.tab; done;
for filename in `cat null.Ht.log`; do grep -H "BEST LOG-L:" $filename >> null.Ht.tab; done;


sed -r 's/.He.alt./\t/g' alt.He.tab | sed -r 's/.res:BEST LOG-L://g' | awk '{ print $1 "\t" $2 "\t" $3 "\t" "He" }'  |  sort -k1 -V  > alt_mle.He.tab
sed -r 's/.He.null./\t/g' null.He.tab | sed -r 's/.res:BEST LOG-L://g' | awk '{ print $1 "\t" $2 "\t" $3 "\t" "He" }'  |  sort -k1 -V  > null_mle.He.tab
 
sed -r 's/.Ht.alt./\t/g' alt.Ht.tab | sed -r 's/.res:BEST LOG-L://g' | awk '{ print $1 "\t" $2 "\t" $3 "\t" "Ht" }'  |  sort -k1 -V  > alt_mle.Ht.tab
sed -r 's/.Ht.null./\t/g' null.Ht.tab | sed -r 's/.res:BEST LOG-L://g' | awk '{ print $1 "\t" $2 "\t" $3 "\t" "Ht" }'  |  sort -k1 -V  > null_mle.Ht.tab
 

R
```

```R
He_null = read.table("null_mle.He.tab", header = F)  # read tab file
He_alt = read.table("alt_mle.He.tab", header = F)  # read tab file

He_likelihoods <-  merge(He_null, He_alt, by= c('V1', 'V2', 'V4'))

Ht_null = read.table("null_mle.Ht.tab", header = F)  # read tab file
Ht_alt = read.table("alt_mle.Ht.tab", header = F)  # read tab file

Ht_likelihoods <-  merge(Ht_null, Ht_alt, by= c('V1', 'V2', 'V4'))

likelihoods <- rbind(He_likelihoods, Ht_likelihoods)
colnames(likelihoods) <- c('genome_location' , 'repeat', 'species', 'lnull' , 'lalt')


LRT <- -2*(likelihoods$lnull - likelihoods$lalt)
pval <- 1-pchisq(LRT, 1)
l_pvals <- cbind(likelihoods, pval, -log(pval))
write.table(l_pvals, file ="likelihoods.pvals.tab", row.names=F, col.names=T, quote=F) 


# lets's exit the R environment
q()

```

```bash


cd ..
cd ..
cp test/res/likelihoods.pvals.tab urchins.adaptiphy.data

awk '{ print $1 "\t" $2 "\t" $6 "\t" $3}' urchins.adaptiphy.data | sed -r 's/\./:/' > urchins.adaptiphy.edited.data

```




# Compute Evolutionary Ratio (zeta)


```bash
cd /query/MODELS_HKY85_query

for filename in *mod; do grep -H "TREE:" $filename; done > output.hky85.txt

sed -r 's/\.mod:TREE: \(/\t/g' output.hky85.txt | sed -r 's/:/\t/g' | sed -r 's/,\(/\t/g' | sed -r 's/,/\t/g' | sed -r 's/\);//g' |  sed -r 's/\)//g' |  awk '{ print $1  "\t" $3  "\t" $5  "\t" $7 "\t" $8  }' |  sed -r 's/\./:/' | sed 1i"genome_location\tLv\tHt\tHe\tHt_He" > Q.hky85.Branches.tab

cd ..
cd ..

cd ref/MODELS_HKY85_ref

for filename in *mod; do grep -H "TREE:" $filename; done > output.hky85.txt

sed -r 's/\.mod:TREE: \(/\t/g' output.hky85.txt | sed -r 's/:/\t/g' |  sed -r 's/\./:/' |  sed -r 's/\./\t/' | sed -r 's/,\(/\t/g' | sed -r 's/,/\t/g' | sed -r 's/\);//g' |  sed -r 's/\)//g' | awk '{ print $1 "\t" $2  "\t" $4  "\t" $6  "\t" $8 "\t" $9  }' | sed 1i"genome_location\treplicate\tLv\tHt\tHe\tHt_He" > R.hky85.Branches.tab


cd ..
cd ..


cp query/MODELS_HKY85_query/Q.hky85.Branches.tab .
cp ref/MODELS_HKY85_ref/R.hky85.Branches.tab .


wc -l Q.hky85.Branches.tab
wc -l R.hky85.Branches.tab



sed -i -e "1d" Q.hky85.Branches.tab 
sed -i -e "1d" R.hky85.Branches.tab


sort -k1 -V  Q.hky85.Branches.tab > Q.hky85.sorted.tab 
sort -k1 -V  R.hky85.Branches.tab > R.hky85.sorted.tab 

join -t $'\t' -j 1 -o 1.1 2.2 1.2 2.3 1.3 2.4 1.4 2.5 1.5 2.6  Q.hky85.sorted.tab R.hky85.sorted.tab > join_trans.tab

head join_trans.tab
chr1:72806-73556	0	0.12833	0.138674	0.0823195	0.0300661	0.0623223	0.034622	0.12833	0.138674
chr1:72806-73556	1	0.12833	0.150726	0.0823195	0.0353043	0.0623223	0.0376183	0.12833	0.150726
chr1:72806-73556	2	0.12833	0.147326	0.0823195	0.0357258	0.0623223	0.0386072	0.12833	0.147326
chr1:72806-73556	3	0.12833	0.159169	0.0823195	0.0337531	0.0623223	0.0395536	0.12833	0.159169


awk '{  rate1 = $3 / $4 ;   print $1 "\t" $2  "\t" $3  "\t"  $4  "\t"  rate1  "\t" "Lv"  }' join_trans.tab  > PhyloFit.Lv.tab
awk '{  rate2 = $5 / $6 ;   print $1 "\t" $2  "\t" $5  "\t"  $6  "\t"  rate2  "\t" "Ht" }' join_trans.tab  > PhyloFit.Ht.tab
awk '{  rate3 = $7 / $8 ;   print $1 "\t" $2  "\t" $7  "\t"  $8  "\t"  rate3  "\t" "He" }' join_trans.tab  > PhyloFit.He.tab
awk '{  rate3 = $9 / $10 ;   print $1 "\t" $2  "\t" $7  "\t"  $8  "\t"  rate3  "\t" "Ht_He" }' join_trans.tab  > PhyloFit.Ht_He.tab



cat PhyloFit.Lv.tab PhyloFit.Ht.tab PhyloFit.He.tab PhyloFit.Ht_He.tab |  sort -k1 -k2,2n  -V | sed 1i"genome_location\trepeat\tQsubsRate\tRsubsRate\tzeta\turchin" > urchins.phyloFit.data



```


```R


urchins_hyphy = as.data.frame(read.table("urchins.adaptiphy.edited.data", header = T)) # read tab file 
colnames(urchins_hyphy) <- c("genome_location", "replicate",   "pval", "species")
urchins_zeta = as.data.frame(read.table("urchins.phyloFit.data", header = T)) # read tab file 
colnames(urchins_zeta) <- c("genome_location", "replicate",  "QsubsRate" , "RsubsRate"   ,  "zeta", "species")

head(urchins_hyphy)
head(urchins_zeta)


urchins.selection.data <- merge(urchins_zeta, urchins_hyphy, by= c('genome_location', 'replicate', 'species'))
head(urchins.selection.data)

library(reshape)
urchins.selection.wide.data <- reshape(urchins.selection.data, idvar = c("genome_location", "replicate"), timevar = "species", direction = "wide")
head(urchins.selection.wide.data)

write.table(urchins.selection.wide.data , file ="urchins.full_selection.data", row.names=F, col.names=T, quote=F) 


library(dplyr)

means_table <- urchins.selection.wide.data  %>% group_by(genome_location) %>% mutate(fdr.He =p.adjust(pval.He   , method="fdr") , fdr.Ht =p.adjust(pval.Ht   , method="fdr"))  %>% select(genome_location, replicate, zeta.He, pval.He, fdr.He, zeta.Ht, pval.Ht, fdr.Ht) %>% summarize(zeta.He = mean(zeta.He), pval.He = mean(pval.He),  fdr.He= mean(fdr.He), zeta.Ht = mean(zeta.Ht), pval.Ht = mean(pval.Ht),  fdr.Ht= mean(fdr.Ht))

write.table(means_table , file ="urchins.selection.data", row.names=F, col.names=T, quote=F) 

```

Measure local compositional complexity (LCC) of DNA sequence

```bash 

nano dolcc_content.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
python LCC_content.py


sbatch dolcc_content.sh

```

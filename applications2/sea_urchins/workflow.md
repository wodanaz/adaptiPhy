
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
cat query/goodalignments.txt > queries_cerebelum.list


cp queries_seaurchins.list /data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/neutral_set

cd /data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/neutral_set

```



```bash
ano DictGen.py
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
values= "neutral_Set.txt"


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
python DictGen.py

```
Run it:

```bash
sbatch genDict.sh
```


```bash

awk -F"." '{print $1 "." $2 }' goodalignments.txt  > ref.list


for file in `cat ref.list`; do for i in {0..9} ; do awk '{if($1 ~ /^>/){split($1,a,"\t"); print a[1]}else{print}}' $file.fa.$i.ref > $file.$i.ref.fa; done; done


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
				f = open('%s.%s.%s.%i.bf' % (i,virus,lrt,replicate), 'w')
				f.write('random_seed=%i;\n' % k)
				f.write('quer_seq_file= "/data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/query%s.fa";\n' % i)
				f.write('ref_seq_file = "/data/wraycompute/alejo/PS_tests/Sea_urchin_evolution/selection_test/ref/%s.%i.ref.fa";\n' % (i,replicate))
				f.write('fit_repl_count = 20;\n')
				f.write('tree= "((He,Ht),Lv)";\n')
				f.write('fgrnd_branch_name = "%s";\n' % urchin)
				f.write('res_file = "res/%s.%s.%s.%i.res";\n' % (i,urchin,lrt,replicate))
				f.write('#include "%s4-fgrnd_spec.bf";\n' % lrt)






#exit nano ctrl+O ENTER ctrl+x

```


```bash
python bfgenerator_global.py ref.list


for i in {0..9} ; do echo *He.alt.$i.bf >> alt.He.$i.list ; done
for i in {0..9} ; do echo *He.null.$i.bf >> null.He.$i.list ; done

for i in {0..9} ; do echo *Ht.alt.$i.bf >> alt.Ht.$i.list ; done
for i in {0..9} ; do echo *Ht.null.$i.bf >> null.Ht.$i.list ; done
```


The following script will generate launchers for each branch and query


```bash
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
			f.write('for file in `cat %s.%s.%i.list`; do root=`basename $file .%s.%s.%i.bf`; HYPHYMP $file > HYPHY/$root.%s.%s.%i.out;done\n' % (lct,urchin,replicate,urchin,lct,replicate,urchin,lct,replicate))









#exit nano ctrl+O ENTER ctrl+x
```



```bash
Run and create new directiries to save data

python shgenerator.py 





mkdir HYPHY
mkdir res

module load hyphy

for file in *sh ; do sbatch $file ; done
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
phyloFit $file.fa --tree "(Lv, (Ht,He))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85_query/$file; # HKY85 model, It runs fast and it also the model applied in HYPHY
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

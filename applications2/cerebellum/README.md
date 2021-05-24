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
	if fasta.get_alignment_length() > 200 :
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

cd query
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




keys= "queries_k562.list"
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





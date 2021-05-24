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

```bash

cd /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1


for file in *fa; do echo $file >> all.list;done


wc -l all.list

head all.list
```

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
for file in *prunned; do echo $file >> all.prunned.list;done
```



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
	if fasta.get_alignment_length() > 145 :
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











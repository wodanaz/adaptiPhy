# This pipeline shows how to run adaptiphy if the user desires to run local tests of selection using 100 Kb of the 
# flanking sequence of the site of interest

First, load the bedtools module according to your system (ignore if it is preloaded)

```bash
sortBed -i dhs_file.bed > dhs_file.sorted.bed

bedtools flank -i dhs_file.sorted.bed -g hg19.genome -b 50000 > flanks50K.bed    # flanks around dhs region window

sortBed -i flanks50K.bed > flanks50K.sorted.bed

bedtools merge -i flanks50K.sorted.bed > flanks50K.merged.bed

bedtools subtract -a flanks50K.merged.bed -b dhs_file.sorted.bed > flanks50K.clean.bed
```

to make local windows of 300bp around each dhs in the list
```bash
bedtools makewindows -b flanks50K.clean.bed -w 300 >  localsplit_100K.bed
```

To create a features file for each chromosome

```bash
mkdir features100K
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY ; 
	do grep -w $chr localsplit_100K.bed | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V > features100K/$chr.feat.bed; 
done
```

# In order to parse alignments for the local non-functional regions

```bash 
cd /directory/with/UCSC/maf_files/maskedv3/

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
msa_split $chr.masked.maf --refseq $chr.masked.fa --gap-strip ANY -q --in-format MAF --features /directory/where/you/run/test/features100K/$chr.feat.bed --for-features --out-root   /directory/where/you/run/test/local/ref/$chr; done;
```

To parse alignments for the DHS regions of interest


```bash
cd /directory/with/UCSC/maf_files/unmasked/
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;
do msa_split $chr.primate.maf --refseq $chr.fa --gap-strip ANY -q --in-format MAF --features /directory/where/you/run/test/features_query/$chr.feat.bed --for-features --out-root/directory/where/you/run/test/local/query/$chr;
done
```

To filter any NFR alignment containing any match to masked sequence


make a list of all the alignments generated

```bash
for file in *fa; do echo $file >> all.list;done


nano filtering.py
import re # for Regular expressions
import sys
from Bio import AlignIO

#infile=sys.argv[1]

infile = 'all.list'

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

Run script
```bash
python filtering.py 
```

Next, we need to format our search dictionary

Rformat and move to the working directory

```bash
cat ref/goodalignments.txt | awk -F"." '{ print $1 "\t" $2}' | awk -F"-" '{ print $1 "\t" $2 "\t" $3}' > goodalignments_100K_REF.300.bed 

```
Do the same for the query, it is important to redo the query bed list becasue after splitting the genome, the locations in the chromosome are sligthly moved


```bash
cat query/goodalignments.txt | awk -F"." '{ print $1 "\t" $2}' | awk -F"-" '{ print $1 "\t" $2 "\t" $3}' | sort -k1,1 -k2,2 -V > queries.bed  #100K
```


#### Sort the distribution of references


```bash
sortBed -i goodalignments_100K_REF.300.bed  > goodalignments_100K_REF.300.sort.bed 
```

### Create Dictionary of NFRs for each DHS

First, define regions:

```bash
bedtools window -a queries.bed -b goodalignments_100K_REF.300.sort.bed -w 50000 > DHSregions2ref100K.bed
```
# Make a 100K Dictionary
```bash
cat DHSregions2ref100K.bed | wc -l
cat DHSregions2ref100K.bed | head
cat DHSregions2ref100K.bed | awk '{print $1 "." $2 "-" $3 ".fa" "\t" $4 "." $5 "-" $6 ".fa" }'  > DHS2REF_100K.tab
awk '$1 != prev{printf "%s%s",ors,$1; ors=ORS; ofs="\t"} {printf "%s%s",ofs,$2; ofs="\t"; prev=$1} END{print ""}' DHS2REF_100K.tab >  DHS2REF_100K.fa.dict


cp DHS2REF_100K.fa.dict ref

cd ref


```

Now we need to concatenate reference files within 100K. 
To do this, run the following script:

```bash
nano merge_LocAlignments.py
import re # for Regular expressions
from Bio import AlignIO
import sys
import csv
import random
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
d = {}
infile='DHS2REF_100K.fa.dict'

with open(infile, 'rb') as csv_file:
     for row in csv.reader(csv_file, delimiter='\t'):
             d[row[0]] = row[1:]





csv_file.close()



reference = []

for i,j in d.iteritems():
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
	referencelist = open('reference.list', 'a')	
	referencelist.write('%s\t%i\n' % (i,n))
	if n > 3:
		with open('%s.ref' % i, 'w') as write_file:
			AlignIO.write(combined_seq, write_file, 'fasta')
	



```

Run previous script using:

```bash
python merge_LocAlignments.py
```


note: python Adds a description that can be removed with this simple Bash code:

```bash
for file in *fa.ref; do echo $file >> ref.tab; done
for file in `cat ref.tab`; do root=`basename $file .fa.ref`; awk '{if($1 ~ /^>/){split($1,a,"\t"); print a[1]}else{print}}' $file > $root.ref.fa; done


for file in *ref.fa; do echo $file >> ref2.list; done
awk -F"." '{print $1 "." $2 }' ref2.list > ref.list
```

This file generates the batch files that run HYPHY for each query and reference in the cluster, modify accordingly

```bash
nano bfgenerator_local.py
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
			f.write('quer_seq_file= "/path/to/query/files/%s.fa.prunned";\n' % i)
			f.write('ref_seq_file = "/path/to/reference/files//%s.ref.fa";\n' % i)
			f.write('fit_repl_count = 20;\n')
			f.write('tree= "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))";\n')
			f.write('fgrnd_branch_name = "%s";\n' % ape)
			f.write('res_file = "res/%s.%s.%s.res";\n' % (i,ape,lrt))
			f.write('#include "%s4-fgrnd_spec.bf";\n' % lrt)





#exit nano ctrl+O ENTER ctrl+x
```



Run python script

```bash
python bfgenerator_global.py ref.list
```


To Make lists of bf launchers per chromosome

```bash
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
	echo $chr.*hg19.null.bf >> $chr.null.hg19.list;
	echo $chr.*hg19.alt.bf >> $chr.alt.hg19.list;
	echo $chr.*panTro4.null.bf >> $chr.null.panTro4.list;
	echo $chr.*panTro4.alt.bf >> $chr.alt.panTro4.list;	
done
```

Create the following file"

```bash
	
nano shgenerator.py
#!/usr/bin/env python
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





``` bash
python shgenerator.py



mkdir HYPHY
mkdir res

module load hyphy



for file in chr*.list ; do wc -l $file ; done
for file in *alt.hg19.list ; do wc -l $file ; done
for file in *null.hg19.list ; do wc -l $file ; done
for file in *alt.panTro4.list ; do wc -l $file ; done
for file in *null.panTro4.list ; do wc -l $file ; done





mkdir HYPHY
mkdir res

module load hyphy
for file in chr*hg19*sh ; do sbatch $file ; done
for file in chr*panTro4*sh ; do sbatch $file ; done
```

After processing in HyPHy, process data just like global references



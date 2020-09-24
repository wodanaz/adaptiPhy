
# A markdown pipeline of my analysis of selection using SARS-CoV-2 as a reference species


First, I downloaded the current genomes alignment from the NCBI


#Upload data to Hardac

```bash
scp Coronavirus_genomes_final.fasta  xxxxx@your.hpc.system.edu:your/path/covid_19/alignment_final
scp SARS_CoV_2.fasta  xxxxx@your.hpc.system.edu:your/path//covid_19/alignment_final
scp SARS_CoV.fasta  xxxxx@your.hpc.system.edu:your/path//covid_19/alignment_final


grep ">" Coronavirus_genomes_final.fasta
>SARS_CoV_2
>Bat_CoV_RaTG13
>Pa_CoV_Guangdong
>Pa_CoV_Guangxi_P4L
>Bat_CoV_LYRa11
>SARS_CoV
>Bat_CoV_BM48



```



#### make a bed file of the reference genome
```bash
nano sars2.bed
chr	1	30068
```

```bash
module load bedtools2
bedtools makewindows -b  sars2.bed -w 299 -s 150 > split_sars2.bed


```

I made an additional segmentation of 500 bp to test for recombination across trees 

```bash

module load bedtools2
bedtools makewindows -b  sars2.bed -w 499 -s 250 > split_sars2_500.bed


```

#### To create sliding windows alignments

```bash

awk '{ print $1 "\t" $2 - 1 "\t" $3 }' split_sars2.bed > split_map.bed

mkdir query
msa_split Coronavirus_genomes_final.fasta --refseq SARS_CoV_2.fasta --gap-strip ANY -q --in-format FASTA --features split_map.bed  --for-features --out-root query/chr 


cd query
for file in *fa; do echo $file >> all.list;done

sort -k1 -V all.list > all_sort.list


```


create the following python program using nano to check sequence size


```bash
nano filtering.py
import re # for Regular expressions
import sys
from Bio import AlignIO

#infile=sys.argv[1]

infile = 'all_sort.list'

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
	if fasta.get_alignment_length() > 100 :
		maxNs = len(re.findall(ambiguous, str(mydata), re.I))
		maxAsk = len(re.findall(missing, str(mydata), re.I))
		if maxNs > 5:
			print ')-B'
			badsimiosNs.append(i)
		elif maxAsk > 1:
			print ')-;'
			badsimiosAsk.append(i)
		else:
			print '(-8 good'
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


######### Run this python script
```



run it

```bash
module load Anaconda/1.9.2-fasrc01
python filtering.py 
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
values= "goodalignments.txt"


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
		combined_seq = MultipleSeqAlignment([SeqRecord(Seq('', generic_dna),id="SARS_CoV_2"), 
							SeqRecord(Seq('', generic_dna),id="Bat_CoV_RaTG13"),
							SeqRecord(Seq('', generic_dna),id="Pa_CoV_Guangdong"),
							SeqRecord(Seq('', generic_dna),id="Pa_CoV_Guangxi_P4L"),
							SeqRecord(Seq('', generic_dna),id="SARS_CoV"),
							SeqRecord(Seq('', generic_dna),id="Bat_CoV_LYRa11"),
							SeqRecord(Seq('', generic_dna),id="Bat_CoV_BM48")])
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


Concatenate 20 sequences randomly


```bash
module load Anaconda/1.9.2-fasrc01
python DictGen.py


awk -F"." '{print $1 "." $2 }' goodalignments.txt  > ref.list


for file in `cat ref.list`; do for i in {0..9} ; do awk '{if($1 ~ /^>/){split($1,a,"\t"); print a[1]}else{print}}' $file.fa.$i.ref > $file.$i.ref.fa; done; done


cd ..

mkdir test
cp query/ref.list test


cd test
```

in test directory

```
nano bfgenerator_global.py
import sys
import csv
import random
query=sys.argv[1]


with open(query) as f:
    querylist = f.read().splitlines() 





model = ['null','alt']
branches = ['SARS_CoV_2','Bat_CoV_RaTG13','Pa_CoV_Guangdong','Pa_CoV_Guangxi_P4L','Bat_CoV_LYRa11','SARS_CoV' ]

for replicate in range(10):
	for lrt in model:
		for virus in branches:
			for i in querylist:
				k = random.randint(1,1000);
				f = open('%s.%s.%s.%i.bf' % (i,virus,lrt,replicate), 'w')
				f.write('random_seed=%i;\n' % k)
				f.write('quer_seq_file= "your/path//covid_19/alignment_final/query/%s.fa";\n' % i)
				f.write('ref_seq_file = "your/path//covid_19/alignment_final/query/%s.%i.ref.fa";\n' % (i,replicate))
				f.write('fit_repl_count = 20;\n')
				f.write('tree= "(((((SARS_CoV_2,Bat_CoV_RaTG13),Pa_CoV_Guangdong),Pa_CoV_Guangxi_P4L),(Bat_CoV_LYRa11,SARS_CoV)),Bat_CoV_BM48)";\n')
				f.write('fgrnd_branch_name = "%s";\n' % virus)
				f.write('res_file = "res/%s.%s.%s.%i.res";\n' % (i,virus,lrt,replicate))
				f.write('#include "%s4-fgrnd_spec.bf";\n' % lrt)






#exit nano ctrl+O ENTER ctrl+x
```

Make list of batch file runs


```bash
python bfgenerator_global.py ref.list


for i in {0..9} ; do echo *SARS_CoV_2.alt.$i.bf >> alt.SARS_CoV_2.$i.list ; done
for i in {0..9} ; do echo *SARS_CoV_2.null.$i.bf >> null.SARS_CoV_2.$i.list ; done

for i in {0..9} ; do echo *Bat_CoV_RaTG13.alt.$i.bf >> alt.Bat_CoV_RaTG13.$i.list ; done
for i in {0..9} ; do echo *Bat_CoV_RaTG13.null.$i.bf >> null.Bat_CoV_RaTG13.$i.list ; done

for i in {0..9} ; do echo *Pa_CoV_Guangdong.alt.$i.bf >> alt.Pa_CoV_Guangdong.$i.list ; done
for i in {0..9} ; do echo *Pa_CoV_Guangdong.null.$i.bf >> null.Pa_CoV_Guangdong.$i.list ; done

for i in {0..9} ; do echo *Pa_CoV_Guangxi_P4L.alt.$i.bf >> alt.Pa_CoV_Guangxi_P4L.$i.list ; done
for i in {0..9} ; do echo *Pa_CoV_Guangxi_P4L.null.$i.bf >> null.Pa_CoV_Guangxi_P4L.$i.list ; done

for i in {0..9} ; do echo *Bat_CoV_LYRa11.alt.$i.bf >> alt.Bat_CoV_LYRa11.$i.list ; done
for i in {0..9} ; do echo *Bat_CoV_LYRa11.null.$i.bf >> null.Bat_CoV_LYRa11.$i.list ; done

for i in {0..9} ; do echo *SARS_CoV.alt.$i.bf >> alt.SARS_CoV.$i.list ; done
for i in {0..9} ; do echo *SARS_CoV.null.$i.bf >> null.SARS_CoV.$i.list ; done

```



```bash




nano shgenerator.py
import sys
import csv
import random





model = ['null','alt']
branches = ['SARS_CoV_2','Bat_CoV_RaTG13','Pa_CoV_Guangdong','Pa_CoV_Guangxi_P4L','Bat_CoV_LYRa11','SARS_CoV' ]
for replicate in range(10):
	for lct in model:
		for virus in branches:
			f = open('%s.%s.%i.sh' % (virus,lct, replicate), 'w')
			f.write('#!/usr/bin/env bash\n')
			f.write('for file in `cat %s.%s.%i.list`; do root=`basename $file .%s.%s.%i.bf`; HYPHYMP $file > HYPHY/$root.%s.%s.%i.out;done\n' % (lct,virus,replicate,virus,lct,replicate,virus,lct,replicate))









#exit nano ctrl+O ENTER ctrl+x

```

Run



```bash


python shgenerator.py 





mkdir HYPHY
mkdir res

module load hyphy

for file in *sh ; do sbatch $file ; done

```



```

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
phyloFit $file.fa --tree "(Bat_CoV_BM48,((Bat_CoV_LYRa11,SARS_CoV),(Pa_CoV_Guangxi_P4L,(Pa_CoV_Guangdong, (SARS_CoV_2,Bat_CoV_RaTG13)))))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85_query/$file; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x





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
phyloFit $file.$i.ref.fa --tree "(Bat_CoV_BM48,((Bat_CoV_LYRa11,SARS_CoV),(Pa_CoV_Guangxi_P4L,(Pa_CoV_Guangdong, (SARS_CoV_2,Bat_CoV_RaTG13)))))" -i FASTA --subst-mod HKY85 --out-root MODELS_HKY85_ref/$file.$i; # HKY85 model, It runs fast and it also the model applied in HYPHY
done; done #exit nano ctrl+O ENTER ctrl+x



sbatch domodel_ref.sh
sbatch domodel_query.sh



```

# After HyPhy had finished running. We can retrieve data from the res/ directory
cd res
#### we need to model the likelihood ratio test results in a chi-square distribution to obtain p-values for each pair of null v alternative models...
#### to do this, we can use the program phase8.rb .... but maybe I can implement it in R



```bash
for file in *SARS_CoV_2.null.*.res; do echo $file >>null.SARS_CoV_2.log; done;
for file in *SARS_CoV_2.alt.*.res; do echo $file >>alt.SARS_CoV_2.log; done;

for file in *Bat_CoV_RaTG13.null.*.res; do echo $file >>null.Bat_CoV_RaTG13.log; done;
for file in *Bat_CoV_RaTG13.alt.*.res; do echo $file >>alt.Bat_CoV_RaTG13.log; done;

for file in *Pa_CoV_Guangdong.null.*.res; do echo $file >>null.Pa_CoV_Guangdong.log; done;
for file in *Pa_CoV_Guangdong.alt.*.res; do echo $file >>alt.Pa_CoV_Guangdong.log; done;

for file in *Pa_CoV_Guangxi_P4L.null.*.res; do echo $file >>null.Pa_CoV_Guangxi_P4L.log; done;
for file in *Pa_CoV_Guangxi_P4L.alt.*.res; do echo $file >>alt.Pa_CoV_Guangxi_P4L.log; done;


for file in *Bat_CoV_LYRa11.null.*.res; do echo $file >>null.Bat_CoV_LYRa11.log; done;
for file in *Bat_CoV_LYRa11.alt.*.res; do echo $file >>alt.Bat_CoV_LYRa11.log; done;

for file in *SARS_CoV.null.*.res; do echo $file >>null.SARS_CoV.log; done;
for file in *SARS_CoV.alt.*.res; do echo $file >>alt.SARS_CoV.log; done;

```


```bash

wc -l *log


for filename in `cat alt.SARS_CoV_2.log`; do grep -H "BEST LOG-L:" $filename >> alt.SARS_CoV_2.tab; done;
for filename in `cat null.SARS_CoV_2.log`; do grep -H "BEST LOG-L:" $filename >> null.SARS_CoV_2.tab; done;

for filename in `cat alt.Bat_CoV_RaTG13.log`; do grep -H "BEST LOG-L:" $filename >> alt.Bat_CoV_RaTG13.tab; done;
for filename in `cat null.Bat_CoV_RaTG13.log`; do grep -H "BEST LOG-L:" $filename >> null.Bat_CoV_RaTG13.tab; done;

for filename in `cat alt.Pa_CoV_Guangdong.log`; do grep -H "BEST LOG-L:" $filename >> alt.Pa_CoV_Guangdong.tab; done;
for filename in `cat null.Pa_CoV_Guangdong.log`; do grep -H "BEST LOG-L:" $filename >> null.Pa_CoV_Guangdong.tab; done;

for filename in `cat alt.Pa_CoV_Guangxi_P4L.log`; do grep -H "BEST LOG-L:" $filename >> alt.Pa_CoV_Guangxi_P4L.tab; done;
for filename in `cat null.Pa_CoV_Guangxi_P4L.log`; do grep -H "BEST LOG-L:" $filename >> null.Pa_CoV_Guangxi_P4L.tab; done;

for filename in `cat alt.SARS_CoV.log`; do grep -H "BEST LOG-L:" $filename >> alt.SARS_CoV.tab; done;
for filename in `cat null.SARS_CoV.log`; do grep -H "BEST LOG-L:" $filename >> null.SARS_CoV.tab; done;

for filename in `cat alt.Bat_CoV_LYRa11.log`; do grep -H "BEST LOG-L:" $filename >> alt.Bat_CoV_LYRa11.tab; done;
for filename in `cat null.Bat_CoV_LYRa11.log`; do grep -H "BEST LOG-L:" $filename >> null.Bat_CoV_LYRa11.tab; done;


```

# Consolidating tables for all the branches with the exception of the outgroup

```bash
awk '{print $1 "\t" $3}' null.SARS_CoV_2.tab | sort -k1,1 -V  > nulls.SARS_CoV_2.tab
awk '{print $1 "\t" $3}' alt.SARS_CoV_2.tab | sort -k1,1 -V > alts.SARS_CoV_2.tab

awk '{print $1 "\t" $3}' null.Bat_CoV_RaTG13.tab | sort -k1,1 -V  > nulls.Bat_CoV_RaTG13.tab
awk '{print $1 "\t" $3}' alt.Bat_CoV_RaTG13.tab | sort -k1,1 -V > alts.Bat_CoV_RaTG13.tab

awk '{print $1 "\t" $3}' null.Pa_CoV_Guangdong.tab | sort -k1,1 -V  > nulls.Pa_CoV_Guangdong.tab
awk '{print $1 "\t" $3}' alt.Pa_CoV_Guangdong.tab | sort -k1,1 -V > alts.Pa_CoV_Guangdong.tab

awk '{print $1 "\t" $3}' null.Pa_CoV_Guangxi_P4L.tab | sort -k1,1 -V  > nulls.Pa_CoV_Guangxi_P4L.tab
awk '{print $1 "\t" $3}' alt.Pa_CoV_Guangxi_P4L.tab | sort -k1,1 -V > alts.Pa_CoV_Guangxi_P4L.tab

awk '{print $1 "\t" $3}' null.SARS_CoV.tab | sort -k1,1 -V  > nulls.SARS_CoV.tab
awk '{print $1 "\t" $3}' alt.SARS_CoV.tab | sort -k1,1 -V > alts.SARS_CoV.tab

awk '{print $1 "\t" $3}' null.Bat_CoV_LYRa11.tab | sort -k1,1 -V  > nulls.Bat_CoV_LYRa11.tab
awk '{print $1 "\t" $3}' alt.Bat_CoV_LYRa11.tab | sort -k1,1 -V > alts.Bat_CoV_LYRa11.tab


```

#### Now, I can make each likelihood table to be joined in a single file

```bash
awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.SARS_CoV_2.tab > col1.SARS_CoV_2.tab
paste col1.SARS_CoV_2.tab nulls.SARS_CoV_2.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6   }' > lnulls.SARS_CoV_2.tab
paste lnulls.SARS_CoV_2.tab alts.SARS_CoV_2.tab |awk '{print "chr" ":" $2  "\t" $4 "\t" $5 "\t" $7 "\t" "SARS_CoV_2" }'  |  sort -k1 -V  > likelihoods.SARS_CoV_2.tab
 
awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.Bat_CoV_RaTG13.tab > col1.Bat_CoV_RaTG13.tab
paste col1.Bat_CoV_RaTG13.tab nulls.Bat_CoV_RaTG13.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6   }' > lnulls.Bat_CoV_RaTG13.tab
paste lnulls.Bat_CoV_RaTG13.tab alts.Bat_CoV_RaTG13.tab |awk '{print "chr" ":" $2  "\t" $4 "\t" $5 "\t" $7 "\t" "Bat_CoV_RaTG13" }'  |  sort -k1 -V  > likelihoods.Bat_CoV_RaTG13.tab
 
awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.Pa_CoV_Guangdong.tab > col1.Pa_CoV_Guangdong.tab
paste col1.Pa_CoV_Guangdong.tab nulls.Pa_CoV_Guangdong.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6   }' > lnulls.Pa_CoV_Guangdong.tab
paste lnulls.Pa_CoV_Guangdong.tab alts.Pa_CoV_Guangdong.tab |awk '{print "chr" ":" $2  "\t" $4 "\t" $5 "\t" $7 "\t" "Pa_CoV_Guangdong" }'  |  sort -k1 -V  > likelihoods.Pa_CoV_Guangdong.tab
 
awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.Pa_CoV_Guangxi_P4L.tab > col1.Pa_CoV_Guangxi_P4L.tab
paste col1.Pa_CoV_Guangxi_P4L.tab nulls.Pa_CoV_Guangxi_P4L.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6   }' > lnulls.Pa_CoV_Guangxi_P4L.tab
paste lnulls.Pa_CoV_Guangxi_P4L.tab alts.Pa_CoV_Guangxi_P4L.tab |awk '{print "chr" ":" $2  "\t" $4 "\t" $5 "\t" $7 "\t" "Pa_CoV_Guangxi_P4L" }'  |  sort -k1 -V  > likelihoods.Pa_CoV_Guangxi_P4L.tab
 
awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.SARS_CoV.tab > col1.SARS_CoV.tab
paste col1.SARS_CoV.tab nulls.SARS_CoV.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6   }' > lnulls.SARS_CoV.tab
paste lnulls.SARS_CoV.tab alts.SARS_CoV.tab |awk '{print "chr" ":" $2  "\t" $4 "\t" $5 "\t" $7 "\t" "SARS_CoV" }'  |  sort -k1 -V  > likelihoods.SARS_CoV.tab
 
awk -F"." '{print $1 "\t" $2 "\t" $3 "\t" $5  }'  nulls.Bat_CoV_LYRa11.tab > col1.Bat_CoV_LYRa11.tab
paste col1.Bat_CoV_LYRa11.tab nulls.Bat_CoV_LYRa11.tab  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6   }' > lnulls.Bat_CoV_LYRa11.tab
paste lnulls.Bat_CoV_LYRa11.tab alts.Bat_CoV_LYRa11.tab |awk '{print "chr" ":" $2  "\t" $4 "\t" $5 "\t" $7 "\t" "Bat_CoV_LYRa11" }'  |  sort -k1 -V  > likelihoods.Bat_CoV_LYRa11.tab
 

 
cat likelihoods.SARS_CoV_2.tab  likelihoods.Bat_CoV_RaTG13.tab  likelihoods.Pa_CoV_Guangdong.tab  likelihoods.Pa_CoV_Guangxi_P4L.tab likelihoods.SARS_CoV.tab   likelihoods.Bat_CoV_LYRa11.tab  |  sort -k1 -k2,2n -k3,3n | sed 1i"location\treplicate\tlnull\tlalt\tvirus" > likelihoods.tab
```


# Now, compute P-values in R:


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



```bash
awk '{ print $1 "\t" $2 "\t" $6 "\t" $5 }' likelihoods.pvals.tab | sort -k1,1 -V  >  likelihoods.pvals2.tab
sed 1i"genome_location\treplicate\tpval\tvirus" likelihoods.pvals2.tab > sars2.adaptiphy.data


cd ..
cd ..
cp test/res/sars2.adaptiphy.data .

sed -i -e "1d" sars2.adaptiphy.data

```

# Example:

For Ref
cat chr.29551-29850.9.mod
ALPHABET: A C G T 
ORDER: 0
SUBST_MOD: HKY85
TRAINING_LNL: -23760.951000
BACKGROUND: 0.303787 0.188364 0.199861 0.307988 
RATE_MAT:
  -0.900858    0.147594    0.511939    0.241325 
   0.238034   -1.183540    0.156602    0.788903 
   0.778145    0.147594   -1.167064    0.241325 
   0.238034    0.482490    0.156602   -0.877126 
TREE: (Bat_CoV_BM48:0.0780206,((Bat_CoV_LYRa11:0.0452548,SARS_CoV:0.0454444):0.0690806,(Pa_CoV_Guangxi_P4L:0.0813149,(Pa_CoV_Guangdong:0.0475241,(SARS_CoV_2:0.0155872,Bat_CoV_RaTG13:0.0189104):0.0355897):0.0463466):0.0908438):0.0780206);


#### Let's stop here for a while.... we need to check the table and make sure it's looking good
#### we need to further modify column 1 alone -> I don't like the dots, it would be nicer to have chromosome and location in separate columns


```bash
for filename in *mod; do grep -H "TREE:" $filename; done > output.hky85.txt

cat output.hky85.txt | sed -r 's/\(+/ /g' |  sed -r 's/\)/ /g'  | sed -r 's/:/ /g' |  sed -r 's/,/ /g' | sed -r 's/;//g' | awk '{ print $1  "\t" $15  "\t" $17  "\t" $13 "\t" $11  "\t" $8 "\t" $6 }' > BranchLenghts.tab

awk -F"." '{print $1 "\t" $2 "\t" $3 }' BranchLenghts.tab  > chr_pos.tab

paste chr_pos.tab BranchLenghts.tab  | column -s '\t' -t | awk '{print $1  ":" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }'  > Branches.tab

sed 1i"genome_location\treplicate\tSARS_CoV_2\tBat_CoV_RaTG13\tPa_CoV_Guangdong\tPa_CoV_Guangxi_P4L\tSARS_CoV\tBat_CoV_LYRa11" Branches.tab > R.hky85.Branches.tab



```



```bash


for filename in *mod; do grep -H "TREE:" $filename; done > output.hky85.txt

cat output.hky85.txt | sed -r 's/\(+/ /g' |  sed -r 's/\)/ /g'  | sed -r 's/:/ /g' |  sed -r 's/,/ /g' | sed -r 's/;//g' | awk '{ print $1  "\t" $15  "\t" $17  "\t" $13 "\t" $11  "\t" $8 "\t" $6 }' > BranchLenghts.tab

awk -F"." '{print $1 "\t" $2  }' BranchLenghts.tab  > chr_pos.tab

paste chr_pos.tab BranchLenghts.tab  | column -s '\t' -t | awk '{print $1  ":" $2  "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 }'  > Branches.tab

sed 1i"genome_location\tSARS_CoV_2\tBat_CoV_RaTG13\tPa_CoV_Guangdong\tPa_CoV_Guangxi_P4L\tSARS_CoV\tBat_CoV_LYRa11" Branches.tab > Q.hky85.Branches.tab


```

```bash

cd ..
cd ..


cp query/MODELS_HKY85_query/Q.hky85.Branches.tab .
cp query/MODELS_HKY85_ref/R.hky85.Branches.tab .


wc -l Q.hky85.Branches.tab
wc -l R.hky85.Branches.tab



sed -i -e "1d" Q.hky85.Branches.tab 
sed -i -e "1d" R.hky85.Branches.tab


sort -k1 -V  Q.hky85.Branches.tab > Q.hky85.sorted.tab 
sort -k1 -V  R.hky85.Branches.tab > R.hky85.sorted.tab 

join -t $'\t' -j 1 -o 1.1 2.2 1.2 2.3 1.3 2.4 1.4 2.5 1.5 2.6 1.6 2.7 1.7 2.8 1.8 2.9   Q.hky85.sorted.tab R.hky85.sorted.tab   > join_trans.tab

#sed 1i"genome_location\tSARS_CoV_2\tBat_CoV_RaTG13\tPa_CoV_Guangdong\tPa_CoV_Guangxi_P4L\tSARS_CoV\tBat_CoV_LYRa11" Branches.tab > Q.hky85.Branches.tab


head join_trans.tab
chr:1-300	0	1.7484e-16	0.0169358	0.0112381	0.0179126	0.00113576	0.0456392	0.0329711	0.0703059	6.82968e-18	0.0507684	0.0284754	0.0460905		
chr:1-300	1	1.7484e-16	0.0170609	0.0112381	0.0184921	0.00113576	0.0632499	0.0329711	0.0852851	6.82968e-18	0.0333493	0.0284754	0.0515742


#SARS_CoV_2_A2a	RatG13	Node3	PangolinCoV	Node2	LYRa11	SARS_CoV	Node7

awk '{  rate1 = $3 / $4 ;   print $1 "\t" $2  "\t" $3  "\t"  $4  "\t"  rate1  "\t" "SARS_CoV_2"  }'             join_trans.tab  > 	PhyloFit.SARS_CoV_2.tab
awk '{  rate2 = $5 / $6 ;   print $1 "\t" $2  "\t" $5  "\t"  $6  "\t"  rate2  "\t" "Bat_CoV_RaTG13" }'          join_trans.tab  > 	PhyloFit.Bat_CoV_RaTG13.tab
awk '{  rate3 = $7 / $8 ;   print $1 "\t" $2  "\t" $7  "\t"  $8  "\t"  rate3  "\t" "Pa_CoV_Guangdong" }'        join_trans.tab  > 	PhyloFit.Pa_CoV_Guangdong.tab
awk '{  rate4 = $9 / $10 ;  print $1 "\t" $2  "\t" $9  "\t"  $10 "\t"  rate4  "\t" "Pa_CoV_Guangxi_P4L" }'      join_trans.tab  > 	PhyloFit.Pa_CoV_Guangxi_P4L.tab
awk '{  rate5 = $11 / $12 ; print $1 "\t" $2  "\t" $11 "\t"  $12 "\t"  rate5  "\t" "SARS_CoV"          }'       join_trans.tab  > 	PhyloFit.SARS_CoV.tab
awk '{  rate6 = $13 / $14 ; print $1 "\t" $2  "\t" $13 "\t"  $14 "\t"  rate6  "\t" "Bat_CoV_LYRa11" }'          join_trans.tab  > 	PhyloFit.Bat_CoV_LYRa11.tab



cat PhyloFit.SARS_CoV_2.tab  PhyloFit.Bat_CoV_RaTG13.tab  PhyloFit.Pa_CoV_Guangdong.tab  PhyloFit.Pa_CoV_Guangxi_P4L.tab PhyloFit.SARS_CoV.tab PhyloFit.Bat_CoV_LYRa11.tab |  sort -k1 -k2,2n  -V | sed 1i"genome_location\treplicate\tQsubsRate\tRsubsRate\tzeta\tvirus" > sars2.phyloFit.data

```

# Computing PhastCons

```bash

mkdir -p TREES
rm -f TREES/*      # in case old versions left over
nano dophastcons.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=your_email@email.com
for file in `cat ref.list` ; do for i in {0..9} ; do
phastCons -i FASTA --estimate-trees TREES/$file.$i $file.fa MODELS_HKY85_ref/$file.$i.mod --no-post-probs;
done;
done;



nano dophastcons2.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=your_email@email.com
ls TREES/*.cons.mod > cons.txt
phyloBoot --read-mods '*cons.txt' --output-average ave.cons.mod 
ls TREES/*.noncons.mod > noncons.txt
phyloBoot --read-mods '*noncons.txt' --output-average ave.noncons.mod 




mkdir -p ELEMENTS SCORES
rm -f ELEMENTS/* SCORES/*
rm dophastcons3.sh
nano dophastcons3.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=your_email@email.com
for file in `cat ref.list` ; do 
phastCons -i FASTA --most-conserved ELEMENTS/$file.bed --score $file.fa ave.cons.mod,ave.noncons.mod > SCORES/$file.wig
done;
        
cd SCORES

for file in *wig ; do root=`basename $file .wig`;  sed -i -e "1d" $file  ; awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $file >  $root.average.wig ; done
grep "" *.average.wig | sed -r 's/.average.wig:/  /g' | sed -r 's/chr./chr:/g' |  awk '{ print $1 "\t" $2 }' | sort -k1 -V > sars2.phastCons.data
        
cd ..
        
cat ELEMENTS/*.bed | sort -k1,1 -k2,2n > most-conserved.bed



cd ..

cp query/SCORES/sars2.phastCons.data .

```


# Consolidate all datasets


```bash

ls *data


sars2.adaptiphy.data  sars2.phastCons.data  sars2.phyloFit.data




module load R
R



covid19_hyphy = as.data.frame(read.table("sars2.adaptiphy.data", header = F)) # read tab file 
colnames(covid19_hyphy) <- c("genome_location", "replicate",   "pval",          "virus")
covid19_zeta = as.data.frame(read.table("sars2.phyloFit.data", header = T)) # read tab file 
head(covid19_hyphy)
head(covid19_zeta)


head(covid19_zeta)
head(covid19_hyphy)


covid19.selection.data <- merge(covid19_zeta, covid19_hyphy, by= c('genome_location', 'replicate', 'virus'))
head(covid19.selection.data)

library(reshape)
covid19.selection.wide.data <- reshape(covid19.selection.data, idvar = c("genome_location", "replicate"), timevar = "virus", direction = "wide")
head(covid19.selection.wide.data)



PhastCons = as.data.frame(read.table("sars2.phastCons.data", header = F)) # read tab file 


colnames(PhastCons) <- c('genome_location', 'phastCons')
head(PhastCons)


covid19.selection.data2 <- merge(covid19.selection.wide.data, PhastCons, by= 'genome_location' , all.x = T)
head(covid19.selection.data2)




write.table(covid19.selection.data2 , file ="sars2.selection.data", row.names=F, col.names=T, quote=F) 


q()

```

# To run PREQUEL to generate common ancestors

```bash
ls *.mod > mytree.txt
phyloBoot --read-mods '*mytree.txt' --output-average mytree.mod

cd ..

cp query/MODELS_HKY85_ref/mytree.mod .

cat mytree.mod


prequel Coronavirus_genomes_final.fasta mytree.mod anc --no-probs  --keep-gaps

```


Do this In local machine

```bash
nano wuCor1.bed
chr	1	29903


bedtools makewindows -b  wuCor1.bed -w 500 -s 150 > split_wuCor1.bed


```
```bash
awk '{ print $1 "\t" $2 "\t" $3 "\t" $6 }' clades_covid19_Annotations.tab | sed -e "1d" > covid19_snps.bed





sed -r 's/,/\t/g' 5000_SARS_CoV_2_samples.csv > 5000_SARS_CoV_2_samples.tab
sed -ir 's/\%//g' 5000_SARS_CoV_2_samples.tab
awk '{print "chr" "\t" $1 "\t" $1 "\t" $2 / 100 }'  5000_SARS_CoV_2_samples.tab > 5000_SARS_CoV_2_samples.bed


bedtools makewindows -b  wuCor1.bed -w 500 -s 150 > split_wuCor1.bed
bedtools intersect -a split_wuCor1.bed -b 5000_SARS_CoV_2_samples.bed -c | awk '{print  $1 "\t" $2 "\t" $3 "\t" $4 / 500 }'  > Allele_Frequency_Spectrum2.bed
cp  Allele_Frequency_Spectrum2.bed /home/alejo/Desktop/Covid-19/Final


# report the number of snps per window (S)

bedtools intersect -a split_wuCor1.bed -b covid19_snps.bed -c > covid19_snps_window.bed
```


```bash

# report the number of common polymorphisms
awk '{  if ( $6 >= 0.05 )   print $1 "\t" $2 "\t" $3 "\t" $6 }' clades_covid19_Annotations.tab | sed -e "1d" > covid19_common_snps.bed


bedtools intersect -a split_wuCor1.bed -b covid19_common_snps.bed -c > covid19_common_snps_window.bed


awk '{  if ( $6 < 0.05 )   print $1 "\t" $2 "\t" $3 "\t" $6 }' clades_covid19_Annotations.tab | sed -e "1d" > covid19_rare_snps.bed


bedtools intersect -a split_wuCor1.bed -b covid19_rare_snps.bed -c > covid19_rare_snps_window.bed


```

#!/usr/bin/env python
#SBATCH --mail-type=END
#SBATCH --mail-user=youremail@mail.com
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem-per-cpu=2000
import sys
import csv
import random

querylist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10','chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
model = ['null','alt']
branches = ['hg19','panTro4' ]
for lrt in model:
	for ape in branches:
		for i in querylist:
			k = random.randint(1,1000);
			f = open('%s.%s.%s.sh' % (i,ape,lrt), 'w')
			f.write('#!/usr/bin/env bash\n')
			f.write('#SBATCH --mail-type=END\n')
			f.write('#SBATCH --mail-user=youremail@email.com";\n')
			f.write('#SBATCH -N 1\n')
			f.write('#SBATCH --mem-per-cpu=100\n')
			f.write('for file in `cat %s.%s.%s.list`; do root=`basename $file .%s.%s.bf`; HYPHYMP $file > HYPHY/$root.%s.%s.out;done\n' % (i,lrt,ape,ape,lrt,ape,lrt))





#exit nano ctrl+O ENTER ctrl+x

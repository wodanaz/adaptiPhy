import re # for Regular expressions
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.lcc import lcc_simp, lcc_mult

#infile=sys.argv[1]

infile = 'goodalignments.txt'

with open(infile) as f:
    reflist = f.read().splitlines() 



f.close()





for i in reflist:
	myfile= open(i, 'r')
	mydata = myfile.read()
	fasta = AlignIO.read(i, "fasta")
	myfile.close()
	lcc_values_lv = lcc_simp(fasta[0].seq)
	lcc_values_He = lcc_simp(fasta[1].seq)
	lcc_values_Ht = lcc_simp(fasta[2].seq)
	LCC_content = open('LCC_content.tab', 'a')	
	LCC_content.write('%s\t%.2f\t%.2f\t%.2f\n' % (i,lcc_values_lv,lcc_values_He, lcc_values_Ht))

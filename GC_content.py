import re # for Regular expressions
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
	gc_values_human = GC(fasta[0].seq)
	gc_values_chimp = GC(fasta[3].seq)
	gc_values = sorted(GC(rec.seq) for rec in SeqIO.parse(i, "fasta"))
	GC_mean = np.mean(gc_values)
	GC_content = open('GC_content.tab', 'a')	
	GC_content.write('%s\t%.2f\t%.2f\t%.2f\n' % (i,gc_values_human,gc_values_chimp, GC_mean))








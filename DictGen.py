import re # for Regular expressions
from Bio import AlignIO
import sys
import csv
import random
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment




keys= "queries.list"
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








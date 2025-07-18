import re # for Regular expressions
import sys
from Bio import AlignIO

infile=sys.argv[1]


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
	if fasta.get_alignment_length() > 200: 
		maxNs = len(re.findall(ambiguous, str(mydata), re.I))
		maxAsk = len(re.findall(missing, str(mydata), re.I))
		if maxNs > 5:
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

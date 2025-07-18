import re # for Regular expressions
import sys
from Bio import AlignIO

infile=sys.argv[1]

with open(infile) as f:
    reflist = f.read().splitlines() 

def filter_alignment(alignment):
    for record in alignment:
        sequence = str(record.seq).upper()
        if len(sequence) < 200:
            return False
        n_count = sequence.count('N')
        if n_count > 5:
          return False
        amb_count = sequence.count('[(+*)]')
        if amb_count > 1:
          return False
    return True

goodseqs = []
badseqs = []

for i in reflist:
    alignment = AlignIO.read(i, "fasta")
    if filter_alignment(alignment):
        goodseqs.append(i)
        print(f"Filtered alignment {i} saved to goodalignments.txt")
    else:
        badseqs.append(i)
        print(f"Alignment {i} did not meet the filtering criteria. Saved to ambiguous.txt")


############

goodsequences = open('goodalignments.txt', 'w')
for item in goodseqs:
  goodsequences.write("%s\n" % item)
goodsequences.close()

badsequences = open('ambiguous.txt', 'w')
for item in badseqs:
  badsequences.write("%s\n" % item)
badsequences.close()

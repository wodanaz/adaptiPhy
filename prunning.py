nano prunning.py
import re # for Regular expressions
import sys, sets


infile = 'all.list'

with open(infile) as f:
    reflist = f.read().splitlines() 



f.close()


for i in reflist:
	myfile= open(i, 'r')
	data=myfile.read().splitlines()
	myfile.close()
	output=open('%s.prunned' % i, 'w')
	results={}
	seq=""
	compName=""
	for line in data:
		if line[0]==">":
			if compName!="":
				results[compName]=seq
			compName=line[1:]
			seq=""
		else:
			seq=seq+line
	results[compName]=seq
	nsize=len(results[results.keys()[0]])
	ns=len(results.keys())
	indels=sets.Set([])
	prunned={}
	for spec in results.keys():
		for i in range(nsize):
			base=results[spec][i]
			if base.upper() not in ["A","C","G","T"]:
				indels.add(i)
	nsize2=nsize-len(indels)
	for spec in results.keys():
		seq=""
		for i in range(nsize):
			if i not in list(indels):
				seq=seq+results[spec][i]
		prunned[spec]=seq
	for f in range(len(prunned.keys())):
		spec=prunned.keys()[f]
		if i>0: output.write("\n")
		output.write(">"+spec+"\n")
		output.write(prunned[spec])
	output.close()
	

		



######### Run this python script
### IMPORTANT  ------>>>> 
#To run this you need Python 2.7.6

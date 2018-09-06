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
			f.write('quer_seq_file= "query/%s.fa.prunned";\n' % i)
			f.write('ref_seq_file = "ref/%s.ref.fa";\n' % i)
			f.write('fit_repl_count = 20;\n')
			f.write('tree= "((((hg19,panTro4),gorGor3),ponAbe2),rheMac3)";\n')
			f.write('fgrnd_branch_name = "%s";\n' % ape)
			f.write('res_file = "res/%s.%s.%s.res";\n' % (i,ape,lrt))
			f.write('#include "%s3-fgrnd_spec.bf";\n' % lrt)





#exit nano ctrl+O ENTER ctrl+x

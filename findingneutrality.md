

To subset a list of putative neutral elements:


```bash

mkdir -p MODELS_HKY85
nano domodel.sh
#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
for file in *.fa ; 
do root=`basename $file .fa`; 
phyloFit $file --tree "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))" --subst-mod HKY85 --out-root MODELS_HKY85/$root; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x



cd MODELS_HKY85

#For Reference
for filename in *.mod; do grep -H "TREE:" $filename; done > output.hky85.txt
cat output.hky85.txt | awk -F":" '{print $1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11  }' | \
awk -F"," '{print $1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10  }'  | \
awk -F")" '{print $1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10  }'  | \
awk '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $13 }'  > BranchLenghts.tab
# Let's stop here for a while.... we need to check the table and make sure it's looking good

# we need to further modify column 1 alone -> I don't like the dots, it would be nicer to have chromosome and location in separate columns

awk -F"." '{print $1 ":" $2  }' BranchLenghts.tab > chr_pos.tab

paste chr_pos.tab BranchLenghts.tab  | column -s '\t' -t | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }'  > Branches.tab

sed 1i"chromosome\trheMac3\tponAbe2\tgorGor3\tpanTro4\thg19\tPanHomo\tPanHomoGor\tPanHomoGorPon" Branches.tab > Branches.mask.data


```



```R



```

#!/usr/bin/env bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type=END
#SBATCH --mail-user=youremail@mail.com
for file in `cat goodalignments.txt` ; 
do root=`basename $file .ref.fa`; 
phyloFit $file --tree "(rheMac3,(ponAbe2,(gorGor3,(panTro4,hg19))))" --subst-mod HKY85 --out-root MODELS_HKY85/$root; # HKY85 model, It runs fast and it also the model applied in HYPHY
done #exit nano ctrl+O ENTER ctrl+x

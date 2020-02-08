Here, I used the following empirical data to simulate sequence alignments to evolve data under different substitution rates in a bash environment


summary(substitutionRates)

      rheMac3            ponAbe2           gorGor3            panTro4        
   Min.   :0.001769   Min.   :0.00906   Min.   :0.004659   Min.   :0.005273  
   1st Qu.:0.011851   1st Qu.:0.01339   1st Qu.:0.007625   1st Qu.:0.007168  
   Median :0.015966   Median :0.01657   Median :0.009577   Median :0.007807  
   Mean   :0.016051   Mean   :0.01647   Mean   :0.009590   Mean   :0.007812  
   3rd Qu.:0.020318   3rd Qu.:0.01946   3rd Qu.:0.011517   3rd Qu.:0.008439  
   Max.   :0.030650   Max.   :0.02536   Max.   :0.015135   Max.   :0.011455  
      hg19             PanHomo           PanHomoGor       PanHomoGorPon     
 Min.   :0.005218   Min.   :0.000000   Min.   :0.002509   Min.   :0.001769  
 1st Qu.:0.007145   1st Qu.:0.000870   1st Qu.:0.005712   1st Qu.:0.011851  
 Median :0.007765   Median :0.001754   Median :0.008077   Median :0.015966  
 Mean   :0.007789   Mean   :0.001780   Mean   :0.008089   Mean   :0.016051  
 3rd Qu.:0.008422   3rd Qu.:0.002598   3rd Qu.:0.010391   3rd Qu.:0.020318  
 Max.   :0.011348   Max.   :0.004387   Max.   :0.014976   Max.   :0.030650 



# To simulate 1000 alignments evolving neutrally in the background and foreground For the query

```bash
i=1;
while [[ i -le 1000 ]] ;
do
awk -v min_hg=71 -v max_hg=89 -v min_pt=66 -v max_pt=89  -v min_ph=0 -v max_ph=34  -v min_gg=57 -v max_gg=134  -v min_hpg=34 -v max_hpg=127 -v min_pa=103 -v max_pa=224  -v min_hpgp=19 -v max_hpgp=301 -v min_rm=19 -v max_rm=301  'BEGIN{"date +%N"|getline rseed;srand(rseed);close("date +%N"); print (min_hg+rand()*(max_hg-min_hg+1))/10000 "\t" (min_pt+rand()*(max_pt-min_pt+1))/10000 "\t" (min_ph+rand()*(max_ph-min_ph+1))/10000 "\t" (min_gg+rand()*(max_gg-min_gg+1))/10000  "\t" (min_hpg+rand()*(max_hpg-min_hpg+1))/10000  "\t" (min_pa+rand()*(max_pa-min_pa+1))/10000  "\t" (min_hpgp+rand()*(max_hpgp-min_hpgp+1))/10000  "\t" (min_rm+rand()*(max_rm-min_rm+1))/10000 }' >> branch_length.neutral.query.tab
  i=$((i+1));
done;

awk '{print "((((hg19:" $1 ",panTro4:" $2 "):" $3 ",gorGor3:" $4 "):" $5 ",ponAbe2:" $6 "):" $7 ",rheMac3:" $8 ");" }' branch_length.neutral.query.tab > Query_Tree.list
```
####################################################################################################\
# scp ab620@hardac-xfer.genome.duke.edu:/data/wraycompute/alejo/PS_tests/primate/Simulation_new/Experiment2/ref/ref30000/MODELS_HKY85/R.hky85.Branches.tab .

# summary(test30000_phylofit$RPhyloFit)
#   Min.     1st Qu.   Median     Mean     3rd Qu.     Max. 
#   0.005218 0.007157  0.007785   0.007801 0.008432    0.011455 


# To simulate 1000 alignments evolving neutrally in the background and foreground for the reference

```bash
i=1;
while [[ i -le 1000 ]] ;
do
awk -v min_hg=77 -v max_hg=78 -v min_pt=77 -v max_pt=78  -v min_ph=16 -v max_ph=17  -v min_gg=95 -v max_gg=96  -v min_hpg=80 -v max_hpg=81 -v min_pa=164 -v max_pa=165  -v min_hpgp=165 -v max_hpgp=166 -v min_rm=159 -v max_rm=160  'BEGIN{"date +%N"|getline rseed;srand(rseed);close("date +%N"); print (min_hg+rand()*(max_hg-min_hg+1))/10000 "\t" (min_pt+rand()*(max_pt-min_pt+1))/10000 "\t" (min_ph+rand()*(max_ph-min_ph+1))/10000 "\t" (min_gg+rand()*(max_gg-min_gg+1))/10000  "\t" (min_hpg+rand()*(max_hpg-min_hpg+1))/10000  "\t" (min_pa+rand()*(max_pa-min_pa+1))/10000  "\t" (min_hpgp+rand()*(max_hpgp-min_hpgp+1))/10000  "\t" (min_rm+rand()*(max_rm-min_rm+1))/10000 }' >> branch_length.neutral.reference.tab
  i=$((i+1));
done;
awk '{print "((((hg19:" $1 ",panTro4:" $2 "):" $3 ",gorGor3:" $4 "):" $5 ",ponAbe2:" $6 "):" $7 ",rheMac3:" $8 ");" }' branch_length.neutral.reference.tab > Reference_Tree.list
```




```bash
i=1;
for tree in `cat Reference_Tree.list`; do echo $tree > reference.$i.tree ; echo reference.$i.tree >> Reference_Tree.list2 ; i=$((i+1)); done

for tree in `cat Reference_Tree.list2`; do  root=`basename $tree .tree`; ./source/seq-gen -mHKY, -l300 -of  < $tree > ref300/$root.fasta ; done
for tree in `cat Reference_Tree.list2`; do  root=`basename $tree .tree`; ./source/seq-gen -mHKY, -l900 -of  < $tree > ref900/$root.fasta ; done
for tree in `cat Reference_Tree.list2`; do  root=`basename $tree .tree`; ./source/seq-gen -mHKY, -l3000 -of  < $tree > ref3000/$root.fasta ; done
for tree in `cat Reference_Tree.list2`; do  root=`basename $tree .tree`; ./source/seq-gen -mHKY, -l9000 -of  < $tree > ref9000/$root.fasta ; done
for tree in `cat Reference_Tree.list2`; do  root=`basename $tree .tree`; ./source/seq-gen -mHKY, -l15000 -of  < $tree > ref15000/$root.fasta ; done
for tree in `cat Reference_Tree.list2`; do  root=`basename $tree .tree`; ./source/seq-gen -mHKY, -l30000 -of  < $tree > ref30000/$root.fasta ; done


i=1;
for tree in `cat Query_Tree.list`; do echo $tree > query.$i.tree ; echo query.$i.tree >> Query_Tree.list2 ; i=$((i+1)); done
for tree in `cat Query_Tree.list2`; do  root=`basename $tree .tree`; seq-gen -mHKY, -l300 -of  < $tree > $root.fasta ; done


awk '!/^>/ { printf "%s", $0; n = "\n" }  
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' query.1000.fasta2
bioawk -c fastx '{ gsub(/\n/,"",seq); print ">"$name; print $seq }' file.fasta
```


# simulating relaxation of constraint
```bash
awk '{ constrant2 = $2 * 0.1 ;  constrant3 = $3 * 0.1 ;  constrant4 = $4 * 0.1 ;  constrant5 = $5 * 0.1 ;  constrant6 = $6 * 0.1 ;  constrant7 = $7 * 0.1 ;  constrant8 = $8 * 0.1 ;   print  $1 "\t"  constrant2  "\t" constrant3 "\t" constrant4  "\t" constrant5 "\t" constrant6 "\t" constrant7  "\t" constrant8 }' neutral.tab >  relaxcons.hg19.tab
awk '{print "((((hg19:" $1 ",panTro4:" $2 "):" $3 ",gorGor3:" $4 "):" $5 ",ponAbe2:" $6 "):" $7 ",rheMac3:" $8 ");" }' relaxcons.hg19.tab > relaxcons.hg19.tree.list
i=1;
for tree in `cat relaxcons.hg19.tree.list`; do echo $tree > query.$i.tree ; echo query.$i.tree >> relaxcons.hg19.tree.list2 ; i=$((i+1)); done
for tree in `cat relaxcons.hg19.tree.list2`; do  root=`basename $tree .tree`; ./source/seq-gen -mHKY, -l300 -of  < $tree > $root.fasta ; done
```


# simulating positive selection
```bash
awk '{ acc = $1 * 7 ; print acc "\t"  $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7  "\t" $8 }' neutral.tab >  possel.hg19.tab
awk '{print "((((hg19:" $1 ",panTro4:" $2 "):" $3 ",gorGor3:" $4 "):" $5 ",ponAbe2:" $6 "):" $7 ",rheMac3:" $8 ");" }' possel.hg19.tab > possel.hg19.tree.list

i=1;
for tree in `cat possel.hg19.tree.list`; do echo $tree > query.$i.tree ; echo query.$i.tree >> possel.hg19.tree.list2 ; i=$((i+1)); done
for tree in `cat possel.hg19.tree.list2`; do  root=`basename $tree .tree`; ./source/seq-gen -mHKY, -l300 -of  < $tree > $root.fasta ; done
```

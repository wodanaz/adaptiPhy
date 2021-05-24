# Analysing evolution by positive selection in a single dataset of cerebellum atac seq


```bash
mkdir features
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY ; 
	do grep -w $chr cerebellum_atac.bed  | awk '{print $1 "\t" $2 - 1 "\t" $3 }' | sort -k1,1 -k2,2 -V >  features/$chr.feat.bed; 
done

```


To generate alignments for Cerebellum regions

```bash
mkdir query    
cd /data/wraycompute/alejo/PS_tests/primate/unmasked1
# for minus queries
rm do_queries.sh
nano do_queries.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=2000
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;
do msa_split $chr.primate.maf --refseq $chr.fa --gap-strip ANY -q --in-format MAF --features /data/wraycompute/alejo/PS_tests/primate/K562_HyPhy/features/$chr.feat.bed --for-features --out-root /data/wraycompute/alejo/PS_tests/primate/Cerebellum_v1/query/$chr; done
```


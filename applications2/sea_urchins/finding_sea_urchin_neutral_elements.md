# To find neutral elements, we need to mask the genome-wide alignments (maf) using all known functional sequences and repetitive DNA

```bash
cat genes.He.bed Hery_repeats.bed  | awk '{ print $1 "\t" $2 "\t" $3 }'  | sort -k1,1V -k2,2n >  functional.feat.bed

bedtools merge -i functional.feat.bed > functional.feat.merge.bed

mkdir featBED

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ; 
	do grep -w $chr functional.feat.merge.bed  | awk '{print $1 "\t" $2 "\t" $3 }' > featBED/$chr.feat.bed; 
done

```

I removed 4 genes because the following error:

```bash
terminate called after throwing an instance of 'std::out_of_range'
  what():  basic_string::replace: __pos (which is 78561070) > this->size() (which is 78556103)
Aborted
terminate called after throwing an instance of 'std::out_of_range'
  what():  basic_string::replace: __pos (which is 71775098) > this->size() (which is 71715399)
Aborted
terminate called after throwing an instance of 'std::out_of_range'
  what():  basic_string::replace: __pos (which is 56660123) > this->size() (which is 56654506)
Aborted
terminate called after throwing an instance of 'std::out_of_range'
  what():  basic_string::replace: __pos (which is 56861364) > this->size() (which is 56843197)
Aborted

```


```bash
mkdir masked_genome
for seaurchins in chr*.maf ; do
root=`basename $seaurchins .He.maf`;
maf_parse $seaurchins --features featBED/$root.feat.bed --mask-features He,Ht,Lv > masked_genome/$root.masked.maf;   
done


for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21; do
bedtools maskfasta -fi $chr.He.fa -bed featBED/$chr.feat.bed -fo masked_genome/$chr.masked.fa;
done 


```


To extract a set of non-functional elements that can be used as  putatively neutral loci

```bash
module load bedtools2

bedtools random -n 5000000 -l 300 -g He > random.bed

sort -k1,1 -k2,2n random.bed > random.sorted.bed

awk '{ print $1 "\t" $2 - 1 "\t" $3 }' random.sorted.bed > random.300.bed 

wc -l random.300.bed # in masked version
5000000 


mkdir features

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY ; 
	do grep -w $chr random.300.bed | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V > features/$chr.feat.bed; 
done

```





# Docker Instructions


If you have a working installation of Docker and want to run small jobs in a contained working environment.

First, sign in and pull the image from docker hub


```bash
docker login --username=yourusername
docker pull wodenaz/accelerated-evolution-version3:latest

```

Next, make sure you see the image in your profile

```bash

docker image ls

```

Finally, enter the image:

```bash
docker run -it wodenaz/accelerated-evolution-version3

```

Navigate trough the file system and you can even run commands:


```bash
bash dobf.sh

#############################################################################wodenaz/accelerated-evolution-version3#############################################
# 2.2. Create lists of files for each model and species


###
# 2.2.A. Create a list of files to apply LRT method


for file in chr*bf ; do echo $file >> hyphy.list; done


### run HYPHY

for file in `cat hyphy.list`; do root=`basename $file .bf`; HYPHYMP $file > HYPHY/$root.out;done


###
# 2.2.B. or, create a list of batch files per chromosome to run in different nodes

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
	echo $chr.*hg19.null.bf >> $chr.null.hg19.list;
	echo $chr.*hg19.alt.bf >> $chr.alt.hg19.list;
	echo $chr.*panTro4.null.bf >> $chr.null.panTro4.list;
	echo $chr.*panTro4.alt.bf >> $chr.alt.panTro4.list;	
done


##########################################################################################################################
# 2.3. Create launcher files for each chromosome to run batches in parallel


		
python shgenerator.py


### Run HYPHY in parallel launchers

for file in chr20*sh ; do bash $file ; done


#########################################

###########################################3
# Step 3


# Once the HYPHY runs have finished, we need to create a table with pvalue and zeta

cd res


for file in *hg19.null.res; do echo $file >> null.hg19.log; done;
for file in *hg19.alt.res; do echo $file >>alt.hg19.log; done;
for file in *panTro4.null.res; do echo $file >>null.panTro4.log; done;
for file in *panTro4.alt.res; do echo $file >> alt.panTro4.log; done;




nano domodelX.sh
#!/usr/bin/env bash
for filename in `cat alt.panTro4.log`; do grep -H "BEST LOG-L:" $filename >> alt.panTro4.tab; done;
for filename in `cat null.panTro4.log`; do grep -H "BEST LOG-L:" $filename >> null.panTro4.tab; done;
for filename in `cat alt.hg19.log`; do grep -H "BEST LOG-L:" $filename >> alt.hg19.tab; done;
for filename in `cat null.hg19.log`; do grep -H "BEST LOG-L:" $filename >> null.hg19.tab; done;


bash domodelX.sh

#########




# Consolidating tables for 5 branches
awk '{print $1 "\t" $3}' null.hg19.tab | sort -k1,1 -V  > nulls.hg19.tab
awk '{print $1 "\t" $3}' alt.hg19.tab | sort -k1,1 -V > alts.hg19.tab

awk '{print $1 "\t" $3}' null.panTro4.tab | sort -k1,1 -V  > nulls.panTro4.tab
awk '{print $1 "\t" $3}' alt.panTro4.tab | sort -k1,1 -V > alts.panTro4.tab


#####
# 3 COLUMNS of data, I can pull down other parameteres later
awk -F"." '{print $1 ":" $2 "\t" $3 }'  nulls.hg19.tab > col1.hg19.tab
paste col1.hg19.tab nulls.hg19.tab  | awk '{print $1 "\t" $2 "\t" $4   }' > lnulls.hg19.tab
paste lnulls.hg19.tab alts.hg19.tab | awk '{print $1 "\t" $2  "\t" $3 "\t" $5  }' |  sort -k1,1 -V > likelihoods.hg19.tab
 
awk -F"." '{print $1 ":" $2 "\t" $3 }'  nulls.panTro4.tab > col1.panTro4.tab
paste col1.panTro4.tab nulls.panTro4.tab  | awk '{print $1 "\t" $2 "\t" $4   }'> lnulls.panTro4.tab
paste lnulls.panTro4.tab alts.panTro4.tab | awk '{print $1 "\t" $2  "\t" $3 "\t" $5  }' |  sort -k1,1 -V > likelihoods.panTro4.tab
 


cat likelihoods.hg19.tab likelihoods.panTro4.tab  |  sort -k1,1 -V | sed 1i"chromosome\tbranch\tlnull\tlalt" > likelihoods.tab



# to go to R:


nano LRT.R

likelihoods = read.table("likelihoods.tab", header = TRUE)  # read tab file
LRT <- -2*(likelihoods$lnull - likelihoods$lalt)
pval <- 1-pchisq(LRT, 1)
l_pvals <- cbind(likelihoods, pval, -log(pval))
write.table(l_pvals, file ="likelihoods.pvals.tab", row.names=F, col.names=F, quote=F) 



RScript LRT.R



#####
# To compute zeta, we need to run PhyloP from PHAST. 

# if the list of prunned alignments is already created in reference and query directories. please run:

bash domodel.sh 



cd ..
cd ..

cp query/MODELS_HKY85/Q.hky85.Branches.tab .
cp ref/MODELS_HKY85/R.hky85.Branches.tab .


wc -l Q.hky85.Branches.tab
wc -l R.hky85.Branches.tab



sed -i -e "1d" Q.hky85.Branches.tab
sed -i -e "1d" R.hky85.Branches.tab



sort Q.hky85.Branches.tab -k1,1 -V > Q.hky85.Branches.sorted.tab
sort R.hky85.Branches.tab -k1,1 -V > R.hky85.Branches.sorted.tab




####


nano zeta_cal.R

rate_Q = as.data.frame(read.table("Q.hky85.Branches.sorted.tab", header = F)) # read tab file 
#colnames(fibroblast_data) <- c('chromosome', 'pval.human', 'zeta.human', 'pval.chimp', 'zeta.chimp','phastcons')
rate_R= as.data.frame(read.table("R.hky85.Branches.sorted.tab", header = F)) # read tab file 
colnames(rate_Q) 
colnames(rate_R) 

zeta <- merge(rate_Q, rate_R, by= "V1")

rate1 <- zeta$V6.x / zeta$V6.y # human
rate2 <- zeta$V5.x / zeta$V5.y # chimp
#rate3 <- fibroblast.zeta$V7.x / fibroblast.zeta$V7.y # node6
#rate4 <- fibroblast.zeta$V4.x / fibroblast.zeta$V4.y # gorilla
#rate5 <- fibroblast.zeta$V8.x / fibroblast.zeta$V8.y # node4
#rate6 <- fibroblast.zeta$V3.x / fibroblast.zeta$V3.y # orangutan


ncHAE_zeta <- data.frame(zeta$V1 , rate1, rate2)
colnames(ncHAE_zeta) <- c('chromosome', 'zeta.human', 'zeta.chimp')


write.table(ncHAE_zeta, file ="PhyloFit.ncHAE.data", row.names=F, col.names=T, quote=F) 



### 

Rscript zeta_cal.R



####



awk '{ print $1 "\t" $2 "\t" $5 }' likelihoods.pvals.tab | grep 'hg19' > ncHAE.human.data
awk '{ print $1 "\t" $2 "\t" $5 }' likelihoods.pvals.tab | grep 'panTro4' > ncHAE.chimp.data



####

R



human_1 = as.data.frame(read.table("ncHAE.human.data", header = F)) # read tab file 
chimp_2 = as.data.frame(read.table("ncHAE.chimp.data", header = F)) # read tab file


 
head(human_1)
head(chimp_2)


human_1$V2 <- NULL

chimp_2$V2 <- NULL



ncHAE.selection.data1 <- merge(human_1, chimp_2, by= "V1")


colnames(ncHAE.selection.data1) <-  c('chromosome', 'pval.human', 'pval.chimp')


rates = as.data.frame(read.table("PhyloFit.ncHAE.data", header = T)) # read tab file   
head(rates)
dim(rates)

ncHAE.selection.data <- merge(ncHAE.selection.data1, rates, by= c("chromosome"))

write.table(ncHAE.selection.data, file ="ncHAE.selection.data", row.names=F, col.names=T, quote=F) 



#q()


# FLAGGING WITH THE GENOME BLACK LIST 



```



Finally, exit close Docker and delete images if necessary


```bash

exit

docker rm -f $(docker ps -aq)
docker system prune -a

```







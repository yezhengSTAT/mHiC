#!/bin/bash

## mHi-C
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu

'''
Script to call each step of mHi-C
April 2016
'''
## ************************************************
## step 0 - Download raw data - Example shown here.
## ************************************************

echo "Start step 0 - downloading!"

id="SRR1658591"
sraDir="/projects/sratoolkit.2.8.2-1-centos_linux64"
path="/projects/fastqFiles"
mkdir -p $path

## tar -zxvf sratoolkit.tar.gz
wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX764/SRX764954/SRR1658591/SRR1658591.sra -O $path/$id.sra
$sraDir/bin/fastq-dump -F --split-files $path/$id.sra -O $path

## step 1-3: Can be run in parallel.

## ******************
## step 1: Alignment
## ******************
name="IMR90_rep1"
ref="/projects/ReferenceGenome/hg19.fasta"
bwaDir="/projects/Softwares/bwa-0.5.9"
samtoolsDir="/projects/Softwares/samtools-1.3"
fastqDir="/projects/fastqFiles"
resultsDir="/projects/IMR90"
bin="/projects/bin"
cutsite="AAGCTAGCTT" # for HindIII
seqLength=25
resolution=40000

## compile cutsite to trim chimeric reads
g++ -std=c++0x -o $bin/cutsite_trimming_mHiC $bin/cutsite_trimming_mHiC.cpp

echo "Start step 1 - alignment!"
bash s1_bwaAlignment.sh "$name" "$ref" "$bwaDir" "$samtoolsDir" "$fastqDir" "$resultsDir/s1" "$bin" 8 "$cutsite" "$seqLength" "$resultsDir/mHiC.summary_w${resolution}_s1"


## **************************
## step 2: Read ends pairing
## **************************
name="IMR90_rep1"
resultsDir="/projects/IMR90"
resolution=40000

echo "Start step 2 - joining read ends!"
python s2_joinEnd.py -r1 ${resultsDir}/s1/${name}_1.sam -r2 ${resultsDir}/s1/${name}_2.sam -o ${resultsDir}/s2/${name}.sam -sf $resultsDir/mHiC.summary_w${resolution}_s2


## *********************************
## step 3: Valid fragment filtering
## *********************************
name="IMR90_rep1"
resultsDir="/projects/IMR90"
refrag="HindIII_resfrag_hg19.bed" #restriction fragment file
resolution=40000
lowerBound=$((resolution * 2))
refragL=50 #$((seqLength * 2))
refragU=500

echo "Start step 3 - categorize read pairs!"
python s3_categorizePairs.py -f ${bin}/${refrag} -r ${resultsDir}/s2/${name}.sam -o ${resultsDir}/s3 -l $refragL -u $refragU -d $lowerBound -m "window" -b $resolution -sf $resultsDir/mHiC.summary_w${resolution}_s3

## In case, chrM is not needed in downstream analysis
# awk -v OFS="\t" '$2 != "chrM" && $7!="chrM" {print $0}' $validP >$validP.noChrM
# rm -rf $validP
# mv $validP.noChrM $validP

## ***************************************
## step 4 - Remove duplicates and binning.
## ***************************************
name="IMR90_rep1"
resultsDir="/projects/IMR90"
resolution=40000
bin="/projects/bin"
validP="${resultsDir}/s3/w${resolution}/${name}.validPairs"
validI="${resultsDir}/s4/w${resolution}/${name}.validPairs"
mappFile="${bin}/human-hg19.HindIII.w${resolution}"
minMap=0.5 #min mappability threshold
minCount=1 #min contact counts allowed
maxIter=150

echo "Start step 4 - duplicates removal and binning!"
bash s4_bin.sh "$validP" "$validI" "$bin" "$mappFile" "$minMap" "$minCount" "$maxIter" "$resultsDir/mHiC.summary_w${resolution}_s4"


## **********************
## step 5 - Build prior.
## **********************
name="IMR90_rep1"
resultsDir="/projects/IMR90"
resolution=40000
validI="${resultsDir}/s4/w${resolution}/${name}.validPairs"
splineBin=200
priorName="uniPrior"

echo "Starts step 5 - prior construction based on uni-reads only!"
python s5_prior.py -f $validI.binPair.Marginal -i $validI.binPairCount.uni.afterICE -o ${resultsDir}/s5 -b $splineBin -l $priorName


## ************************************************************************************
## step 6 - Generative model to assign probability to multi-reads potential alignments.
## ************************************************************************************
name="IMR90_rep1"
resultsDir="/projects/IMR90"
resolution=40000
prior="${resultsDir}/s5/s5_w${resolution}_splineResults"
multi="${resultsDir}/s4/${name}.validPairs.MULTI.binPair.multi"
multiKeys="$resultsDir/s4/${name}.validPairs.MULTI.binPair.multiKeys" 
uni="$resultsDir/s4/${name}.validPairs.binPairCount.uni"
filename="${name}.validPairs.binPair.multi"
threshold=0.5

echo "Starts step 6 - assign probability to multi-reads potential alignment positions !"
awk -v OFS="_" '{print $2, $3, $4, $5}' $multi | sort -u >$multiKeys
python s6_em.py -p $prior -u $uni -m $multi -mk $multiKeys -t $threshold -o "${resultsDir}/s6" -f $filename

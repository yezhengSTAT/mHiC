#!/bin/bash

## mHi-C
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu

'''
Script to call each step of mHi-C
April 2016
'''

## Step 0 - Download raw data - Example shown here.

echo "Start step 0 - downloading!"

id="SRR1658591"
sraDir=$(pwd)/sratoolkit.2.8.2-1-centos_linux64
path=$(pwd)/raw
mkdir -p $path

## tar -zxvf sratoolkit.tar.gz
wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX764/SRX764954/SRR1658591/SRR1658591.sra -O $path/$id.sra
$sraDir/bin/fastq-dump -F --split-files $path/$id.sra -O $path


## Step 1-3: Can be run in parallel.

name=$1 #Can be replicate, e.g "rep1"
mkdir $name
resultsDir=$(pwd)/$name

sequencesDir=$(pwd)/Sequences # Given
cutsite="GATCGATC" # for MboI
seqLength=25

resolution=5000
validP="${resultsDir}/s3/w${resolution}/${name}.validPairs"
validI="${resultsDir}/s4/w${resolution}/${name}.validPairs"
refrag="MboI_resfrag_hg19.bed" #restriction fragment file
minMap=0.5
minCount=1
lowerBound=20000 #$((resolution * 2))
refragL=50 #$((seqLength * 2))
refragU=500

g++ -std=c++0x -o $bin/cutsite_trimming_mHiC $bin/cutsite_trimming_mHiC.cpp

echo "Prepare Output Directory and Summary File!"

echo "Start step 1 - alignment!"
bash s1_bwaAlignment.sh $name $ref $bwaDir $samtoolsDir $sequencesDir $resultsDir $bin 1 $cutsite $seqLength $resultsDir/mHiC.summary_w${resolution}

echo "Start step 2 - joining read ends!"
python s2_joinEnd.py -r1 ${resultsDir}/s1/${name}_1.sam -r2 ${resultsDir}/s1/${name}_2.sam -o ${resultsDir}/s2/${name}.sam -sf $resultsDir/mHiC.summary_w${resolution}

echo "Start step 3 - categorize read pairs!"
python s3_categorizePairs.py -f ${bin}/$refrag -r ${resultsDir}/s2/${name}.sam -o ${resultsDir}/s3/w$resolution -l $refragL -u $refragU -d $lowerBound -m "window" -b $resolution -sf $resultsDir/mHiC.summary_w${resolution}

## In case, chrM is not needed in downstream analysis
# awk -v OFS="\t" '$2 != "chrM" && $7!="chrM" {print $0}' $validP >$validP.noChrM
# rm -rf $validP
# mv $validP.noChrM $validP

## Step 4 - Remove duplicates and binning. Also ICE normalization on uni-reads are process in preparation for step 5 prior building.

echo "Start step 4 - duplicates removal and binning!"
maxIter=150
bin=$(pwd)/bin

cellLine="GM12878"
resolution=5000
validP="${resultsDir}/s4/w${resolution}/${name}.validPairs"
validI="${resultsDir}/s4/w${resolution}/${name}.validPairs"
refrag="MboI_resfrag_hg19.bed" #restriction fragment
minMap=0.5 #min mappability threshold
minCount=1 #min contact counts allowed
splineBin=200 # fit-hi-c spline bin number
refragL=50
refragU=500

mappability="null"


if [ ! -d "$resultsDir/s4/w${resolution}" ]; then
    mkdir -p "$resultsDir/s4/w${resolution}"
fi


awk -v OFS="\t" '{print $1, $2, $6, $7, $11}' $validP.UNI >$validP.UNI.binPair

echo "processing!"
bash s4_bin.sh "$validP" "$validI" "$bin" "$mappability" "$minMap" "$minCount" "$maxIter" "$resultsDir/mHiC.summary_s4_w${resolution}"


## step 5 - build prior
echo "Starts step 5 - prior construction based on uni-reads only!"
if [ ! -d "${resultsDir}/s5" ]; then
    mkdir -p "${resultsDir}/s5"
fi


python s5_prior.py -l uniPrior -f $validI.binPair.Marginal -i $validI.binPairCount.uni.afterICE -o ${resultsDir}/s5 -b $splineBin



## step 6 - Generative model to assign probability to multi-reads potential alignments.

echo "Starts step 6 - assign probability to multi-reads potential alignment positions !"


priorPath=$resultsDir/s5/s5_w${resolution}_splineResults 
multiFilePath=$resultsDir/s4/${name}.validPairs.MULTI.binPair.multi 
multiKeysPath=$resultsDir/s4/${name}.validPairs.MULTI.binPair.multiKeys 
uniFilePath=$resultsDir/s4/${name}.validPairs.binPairCount.uni 
filename=${name}.validPairs.binPair.multi
threshold=0.5

awk -v OFS="_" '{print $2, $3, $4, $5}' $multiFilePath | sort -u >$multiKeysPath

python s6_em.py -p $priorPath -u $uniFilePath -m $multiFilePath -mk $multiKeysPath -t $threshold -o ${resultsDir}/s6 -f $filename


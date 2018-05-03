#!/bin/bash

## mHi-C
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu

'''
Script to remove duplicates and binning.
ICE normalization of uni-reads are conducted in preparation
for prior building in step 5.

'''



################################################################
# step 4 remove duplicates among Valid Interactions and binning#
################################################################
validP=$1
validI=$2
bin=$3
dir=${validI%/*}
mappability=$4
minMap=${5:-0.5}
minCount=${6:-1}
maxIter=${7:-100}
summaryFile=${8:-""}

if [ ! -d $dir/sorttmp ]; then
    mkdir -p $dir/sorttmp
fi

# awk -v OFS="\t" '{print $1, $2, $6, $7, $11}' $validP.UNI >$validP.UNI.binPair

# Remove PCR duplicates based on alignment chrom + position
sort -k2,2V -k3,3n -k7,7V -k8,8n -T $dir/sorttmp  $validP  | awk -v OFS="\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$7 || s1!=$3 || s2!=$8){print;c1=$2;c2=$7;s1=$3;s2=$8}' | sort -k1,1V -k2,2V -k6,6n -k7,7V -k11,11n -T $dir/sorttmp > $validI.nodup

rm -rf $dir/sorttmp


# Take only the read ID and bin pairs information, then remove bin pair duplicates of one read pair
awk -v OFS="\t" 'BEGIN{q=0;c1=0;c2=0;b1=0;b2=0}(q!=$1 || c1!=$2 || c2!=$7 || b1!=$6 || b2!=$11){print $1, $2, $6, $7, $11;q=$1;c1=$2;c2=$7;b1=$6;b2=$11}' $validI.nodup | awk -v OFS="\t" '{split($2, chr1, "chr"); split($4, chr2, "chr"); if(chr1[2] > chr2[2]) {print $1, $4, $5, $2, $3} else {print $0}}' > $validI.nodup.binPair


# Split valid interaction bin pair files into unique bin pair and multiple bin pair file
awk '{print $1}' $validI.nodup.binPair | sort | uniq -c >$validI.nodup.binPair.uniMultiSummary

awk '$1 >1 {print $2}' $validI.nodup.binPair.uniMultiSummary >$validI.nodup.binPair.multiList

awk 'FNR==NR{a[$1];next};!($1 in a)' $validI.nodup.binPair.multiList $validI.nodup.binPair >$validI.nodup.binPair.uni

awk 'FNR==NR{a[$1];next};($1 in a)' $validI.nodup.binPair.multiList $validI.nodup.binPair >$validI.nodup.binPair.multi

# Uni bin pair count
awk '{print $2, $3, $4, $5}' $validI.nodup.binPair.uni | sort | uniq -c | awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | sort -k1,1V -k2,2n -k3,3V -k4,4n >$validI.binPairCount.uni



# Summary of removing duplicates
if [ "${#summaryFile}" -eq "0" ]; then
    summaryFile=$dir/rmDuplicates.summary
fi
echo -e "\n" >>$summaryFile
echo -e "\n" >>$summaryFile
echo -e "Step: Remove duplicates alignments" >>$summaryFile

validI_all=$(cat $validP | wc -l)
validI_nodup=$(cat $validI.nodup | wc -l)
validI_rmdup=$(( validI_all - validI_nodup ))
validI_nodupBin=$(cat $validI.nodup.binPair | wc -l)
readCount_uniBin=$(cat $validI.nodup.binPair.uni | wc -l)
readCount_multiBin=$(cat $validI.nodup.binPair.multiList | wc -l)


echo -e "Valid interaction pairs including uni-read pairs and all kinds of multi-mapping possible pairs:" >>$summaryFile
echo -e "  Valid interaction pairs count:\t"$validI_all >>$summaryFile
echo -e "    Valid interaction pairs without duplicates:\t"$validI_nodup >>$summaryFile
echo -e "    Valid interaction pairs removed duplicates:\t"$validI_rmdup >>$summaryFile
echo -e "  Valid interaction bin pairs without duplicates:\t"$validI_nodupBin >>$summaryFile
echo -e "Valid read pairs count:" >>$summaryFile
echo -e "  Valid read pair count with unique bin pair:\t"$readCount_uniBin >>$summaryFile
echo -e "  Valid read pair count with multiple possible bin pairs:\t"$readCount_multiBin >>$summaryFile

# Exclude low count contacts
if [ "$minCount" -gt "1" ];then

    awk -v minC=$minCount -v OFS="\t" '$5>=minC {print $0}' $validI.binPairCount.uni | sort -k1,1V -k2,2n -k3,3V -k4,4n >$validI.binPairCount.uni.minCount$minCount

    # ICE normalization with filtering low mappability regions
    python $bin/ICE-with-sparseMatrix.py $validI.binPairCount.uni.minCount$minCount $bin/$mappFile l1 $validI.binPairCount.uni.afterICE $minMap $maxIter
else
    # ICE normalization with filtering low mappability regions
    python $bin/ICE-with-sparseMatrix.py $validI.binPairCount.uni $bin/$mappFile l1 $validI.binPairCount.uni.afterICE $minMap $maxIter
fi

# Uni bin marginal pair
awk '{print $1, $2}' $validI.binPairCount.uni.afterICE >$validI.marginal1
awk '{print $3, $4}' $validI.binPairCount.uni.afterICE >$validI.marginal2
cat $validI.marginal1 $validI.marginal2 | sort -k1,1V -k2,2n | awk '!a[$0]++' >$validI.binPair.Marginal
rm -rf $validI.marginal*

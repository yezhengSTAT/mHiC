#!/bin/bash

#######################
# step1 bwa alignment #
#######################
#parameters
name=$1
ref=$2 #/projects/ReferenceGenome/PlasmoDB-9.0_Pfalciparum3D7_Genome.fasta

bwaDir=$3 #/projects/Softwares/bwa-0.5.9
samtoolsDir=$4 #/projects/Softwares/samtools-1.3

fastqDir=$5
resultDir=$6/s1
bin=$7
core=$8
cutsite=$9 #AAGCTAGCTT for HindIII
seqLength=${10:-25} #>=25 enforced by mHiC
summaryFile=${11:-"mHiC.summary"}

## min length for chimeric reads is enforced to be 25bp
if [ "$seqLength" -lt "25" ]; then
    seqLength=25
fi


#Make sure the directories all exist
if [ ! -d "$resultDir" ]; then
    mkdir -p $resultDir
fi

#BWA alignment
for i in 1 2
do
    echo "Start alignment of read end $i"

    echo "Step1.1 - BWA alignment"
    $bwaDir/bwa aln -n 2 -o 1 -q 0 -t $core $ref $fastqDir/$name\_$i.fastq >$resultDir/$name\_$i.sai
    $bwaDir/bwa samse -n 99 $ref $resultDir/$name\_$i.sai $fastqDir/$name\_$i.fastq >$resultDir/$name\_$i.sam

    # step1.2 - filter get aligned sam & unmapped sam
    echo "Step1.2 - Filter and get unmapped alignment sam file."
    ## Filter out the unmapped for future chimeric reads rescuing.
    $samtoolsDir/samtools view -f 4  $resultDir/$name\_$i.sam >$resultDir/$name\_unmapped_$i.sam
    mv $resultDir/$name\_$i.sam $resultDir/$name\_$i\_raw.sam
    
    # step1.3 - trim and filter unmapped
    echo "Step1.3 - Trim unmapped reads until the restriction enzyme cutting site."
    $samtoolsDir/samtools fastq $resultDir/$name\_unmapped_$i.sam >$resultDir/$name\_unmapped_$i.fastq
    $bin/cutsite_trimming_Ye --fastq $resultDir/$name\_unmapped_$i.fastq --cutsite $cutsite --out $resultDir/$name\_unmapped_trim_$i.fastq --rmuntrim

    awk -v minLen=$seqLength 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= minLen) {print header, seq, qheader, qseq}}' < $resultDir/$name\_unmapped_trim_$i.fastq >$resultDir/$name\_unmapped_trim_filter_$i.fastq

    # step4 - align trimed read
    echo "Step1.4 - Rescue chimeric reads by re-aligning trimmed unmapped reads."
    $bwaDir/bwa aln -n 2 -o 1 -q 0 -t $core $ref $resultDir/$name\_unmapped_trim_filter_$i.fastq >$resultDir/$name\_unmapped_trim_filter_$i.sai
    $bwaDir/bwa samse -n 99 $ref  $resultDir/$name\_unmapped_trim_filter_$i.sai $resultDir/$name\_unmapped_trim_filter_$i.fastq >$resultDir/$name\_unmapped_trim_filter_$i.sam
    $samtoolsDir/samtools view $resultDir/$name\_unmapped_trim_filter_$i.sam >$resultDir/$name\_unmapped_trim_filter_noheader_$i.sam
    # step5 - merge two step alignment
    echo "step1.5 - Merge chimeric reads with mapped full-length reads."
    awk 'NR==FNR{a[$1]=$0;next;}a[$1]{$0=a[$1]}1' $resultDir/$name\_unmapped_trim_filter_noheader_$i.sam $resultDir/$name\_$i\_raw.sam >$resultDir/$name\_$i.sam
done

rawAlign1=$(awk '$6 !="*" && $1!="@SQ" && $1!="@PG" {print $1}' $resultDir/${name}_1_raw.sam | wc -l)
rawMulti1=$(awk '{split($20, tag, ":"); if(tag[1] == "XA" && $2!="4"){print $1}}' $resultDir/${name}_1_raw.sam | wc -l)
chimericAlign1=$(awk '$6 != "*" && $1!="@SQ" && $1!="@PG" {print $1}' $resultDir/$name\_unmapped_trim_filter_noheader_1.sam | wc -l)
chimericMulti1=$(awk '{split($20, tag, ":"); if(tag[1] == "XA" && $2!="4"){print $0}}' $resultDir/$name\_unmapped_trim_filter_noheader_1.sam | wc -l)                  

rawAlign2=$(awk '$6 !="*" && $1!="@SQ" && $1!="@PG" {print $1}' $resultDir/${name}_2_raw.sam | wc -l)
rawMulti2=$(awk '{split($20, tag, ":"); if(tag[1] == "XA" && $2!="4"){print $1}}' $resultDir/${name}_2_raw.sam | wc -l)
chimericAlign2=$(awk '$6 != "*" && $1!="@SQ" && $1!="@PG" {print $1}' $resultDir/$name\_unmapped_trim_filter_noheader_2.sam | wc -l)
chimericMulti2=$(awk '{split($20, tag, ":"); if(tag[1] == "XA" && $2!="4"){print $0}}' $resultDir/$name\_unmapped_trim_filter_noheader_2.sam | wc -l)                  

echo -e "First stage alignment aligned reads total 1: "$rawAlign1 >>$summaryFile
echo -e "First stage alignment aligned reads total 2: "$rawAlign2 >>$summaryFile
echo -e "First stage alignment aligned multi-reads 1: "$rawMulti1 >>$summaryFile
echo -e "First stage alignment aligned multi-reads 2: "$rawMulti2 >>$summaryFile
                                                                                              
echo -e "Second stage alignment aligned reads total 1: "$chimericAlign1 >>$summaryFile
echo -e "Second stage alignment aligned reads total 2: "$chimericAlign2 >>$summaryFile
echo -e "Second stage alignment aligned multi-reads 1: "$chimericMulti1 >>$summaryFile
echo -e "Second stage alignment aligned multi-reads 2: "$chimericMulti2 >>$summaryFile

#Remove redundant files
rm -rf $resultDir/*sai
rm -rf $resultDir/*unmapped*
rm -rf $resultDir/*raw*

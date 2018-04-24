# mHi-C: robust leveraging of multi-mapping reads in Hi-C analysis
Zheng, Ye, Ferhat Ay, and Sunduz Keles. "mHi-C: robust leveraging of multi-mapping reads in Hi-C analysis." bioRxiv (2018): 301705.

The pipeline is developed in Keles Research Group in University of Wisconsin - Madison and please contact Ye Zheng (yezheng@stat.wisc.edu) for any question and suggestion.

## What is mHi-C?
mHi-C is short for **m**ulti-mapping strategy for **Hi-C** data in order to make use of reads aligned to multiple positions. mHi-C pipeline was developed to incorporate multi-mapping reads starting from unaligned read files and produces a set of statistically significant contacts at a give resolution. Remarkably, each main step is organized into independent script with flexible user-defined parameters and implementation options. Therefore, analysis can be carried out from any step of the work-flow and easily fits in high performance computing enviroments for parallel computations.

![Contact matrices comparison between Uni-setting and Uni&Multi-setting](/p/keles/ENCODE-TE/volume1/YeZheng/github/mHiC/figures/mHiC_gif.gif "Contact matrices comparison")

## mHi-C main procedures

### Step 0 - Pipeline caller [mhic_step0-6.sh]
Caller for all the steps in mHi-C pipeline, starting from alignment to multi-reads alignment probability assignment. This is a demo script to run multiple steps at once. Parameters in the script should be customize it for you own use.

#### 0.1 Usage

```
bash mhic_step0-6.sh
```

### Step 1 - Alignment [s1_bwaAlignment.sh]
Default aligner is BWA but other aligner such as bowtie can also be used as long as it records multi-mapping reads related alignment information. Besides, chimeric reads are rescued in this step.

#### 1.0 Requirements
- BWA (>=0.5.9)
- samtools (>=1.3)

#### 1.1 Arguments

```
1. name         : Name of the input fastq file. Hi-C reads are paired-end reads stored in a paired fastq files, for example, IMR90_primaryReplicate_1.fastq and IMR90_primaryReplicate_2.fastq. The name should be "IMR90_primaryReplicate".
2. ref          : Reference genome file including the path to it. For example, "/projects/ReferenceGenome/PlasmoDB-9.0_Pfalciparum3D7_Genome.fasta".
3. bwaDir       : Path to BWA aligner. For example, "/projects/Softwares/bwa-0.5.9".
4. samtoolsDir  : Path to samtools. For example, "/projects/Softwares/samtools-1.3".
5. fastqDir     : Path to the fastq files.
6. resultsDir   : Path to the output results including the intermediate files. 
7. bin          : Path to the bin folder where cutsite_trimming_mHiC.cpp can be found.
8. core         : Number of threads for parallel alignment.
9. cutsite      : Restriction enzyme cutting site. For example, "AAGCTAGCTT" for HindIII and "GATCGATC" for MboI.
10. seqLength   : [Optional] The minimum read length for chimeric reads. >=25 is enforced by mHiC.
11. summaryFile : [Optional] Name for the alignment summary file. By default, "mHiC.summary" will be used as the summary file name.
```

#### 1.2 Usage

```
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

bash s1_bwaAlignment.sh "$name" "$ref" "$bwaDir" "$samtoolsDir" "$fastqDir" "$resultsDir/s1" "$bin" 8 "$cutsite" "$seqLength" "$resultsDir/mHiC.summary_w${resolution}_s1"

```

### Step 2 - Read ends pairing [s2_joinEnd.py] 
In step 1, two ends (_1.fastq and _2.fastq) are aligned separetely to the reference genome which can be paired by read ID. Thus paired-end reads can be formed and each paired-end read represent one interaction.

#### 2.0 Requirements
- python3 (>= 3.6)
- numpy (>= 1.13.1)
- scipy (>= 0.19.1)
- pysam (>= 0.12.0)


#### 2.1 Arguments

```
readEnd1	(-r1/--readEnd1)	    : Path to the aligned file of read end 1.
readEnd2	(-r2/--readEnd2)	    : Path to the aligned file of read end 2.
output		(-o/--output)		      : Path to the aligned paired-end file.
summary		(-s/--summary)		    : [Optional] Summarize alignment results, i.e. number of unmapped reads, singleton and mapped reads, etc. Default is true.
summaryFile	(-sf/--summaryFile)	: [Optional] Path to the summary file including the summary file name. Default is alignedFileName.alignSummary
verbose		(-v/--verbose)		    : [Optional] Verbose. Default is true.

```

#### 2.2 Usage
```
name="IMR90_rep1"
resultsDir="/projects/IMR90"
resolution=40000

python3 s2_joinEnd.py -r1 ${resultsDir}/s1/${name}_1.sam -r2 ${resultsDir}/s1/${name}_2.sam -o ${resultsDir}/s2/${name}.sam -sf $resultsDir/mHiC.summary_w${resolution}_s2
```


### Step 3 - Valid fragment filtering [s3_categorizePairs.py]
This step is to ensure valid read pairs are passed on to downstream analysis while excluding dangling end, self circle, religation, too short-range interactions as well as  invalid alignment that are far away from restriction enzyme cutting sites. Read pairs in each category are summarized.

#### 3.0 Requirements
- python3 (>= 3.6)
- numpy (>= 1.13.1)
- scipy (>= 0.19.1)
- pysam (>= 0.12.0)
- bx-python (>= 0.5.0)

#### 3.1 Arguments

```
fragment	(-f/--fragment)		: Restriction Enzyme(RE) fragment file.
read		(-r/--read)		      : Paired-ended read alignment file from step 2.
outdir		(-o/--outdir)		  : Output directory to save category results.
lower		(-l/--lower)		    : Lower bound of the read pair distances summation from read end alignment position to its assigned RE fragment cutting site. Default value is None indicating no restriction. Recommended 50.
upper		(-u/--upper)		    : Upper bound of the read pair distances summation from read end alignment position to its assigned RE fragment cutting site. Default value is None indicating no restriction. Recommended 500.
distance	(-d/--distance)		: Minimum distance between intrachromosomal read pair alignments positions. Default value is None indicating no restriction. Recommended 20k or two times the resolution.
binMethod	(-m/--method)		  : Binning methods: by fixed number of window size, i.e. window, or fixed number of RE fragment, i.e. fragment. By default, binning by window is utilized.
resolution	(-b/--bin)		  : Binning resolution. If binning by fixed number of window size, it can be 10000 or any other suitable bin size. If binning by fixed number of RE fragment, it can be 10 representing 10RE fragments or other suitable bin size.
summary		(-s/--summary)		: [Optional] Summarize the categories of read pairs. Default is true.
summaryFile	(-sf/--summaryFile)	: [Optional] Path to the summary file including the summary file name. Default is alignedFileName.readCategorySummary.
verbose		(-v/--verbose)		: [Optional] Verbose. Default is true.
```

#### 3.2 Usage

```
name="IMR90_rep1"
resultsDir="/projects/IMR90"
refrag="HindIII_resfrag_hg19.bed" #restriction fragment file
resolution=40000
lowerBound=$((resolution * 2))
refragL=50 #$((seqLength * 2))
refragU=500

python s3_categorizePairs.py -f ${bin}/${refrag} -r ${resultsDir}/s2/${name}.sam -o ${resultsDir}/s3 -l $refragL -u $refragU -d $lowerBound -m "window" -b $resolution -sf $resultsDir/mHiC.summary_w${resolution}_s3
```

#### 3.3 Input file - Restriction enzyme fragment (BED file)
```
chr1    0       16007   HIC_chr1_1      0       +
chr1    16007   24571   HIC_chr1_2      0       +
chr1    24571   27981   HIC_chr1_3      0       +
chr1    27981   30429   HIC_chr1_4      0       +
chr1    30429   32153   HIC_chr1_5      0       +
chr1    32153   32774   HIC_chr1_6      0       +
chr1    32774   37752   HIC_chr1_7      0       +
chr1    37752   38369   HIC_chr1_8      0       +
chr1    38369   38791   HIC_chr1_9      0       +
chr1    38791   39255   HIC_chr1_10     0       +
chr1    39255   43602   HIC_chr1_11     0       +
```


### Step 4 - Duplicates removal and genome binning [s4_bin.sh]
Remove the PCR duplicates and bin the genome by fixed window size.

#### 4.0 Requirements
- python3 (>= 3.6)
- numpy (>= 1.13.1)

#### 4.1 Arguments

```
1. validP	      : Path to the valid read pairs obtained from step3.
2. validI	      : Path to the output non-duplicated valid interactions. Users can take this chance to rename the interaction files otherwise users can set it to be the same as validP
3. bin		      : Path to the bin folder where ICE normalization script can be found.
4. mappFile	    : Path to the mappability file.
5. minMapp	    : Minimum mappability of the regions considered. Default is 0.5.
6. minCount	    : Minimum contact counts of the bin pairs considered. Default is 1.
7. maxInter	    : Maximum iteraction of ICE to normalize contact matrix. Default is 100.
8. summaryFile	: Summary file name. Default is rmDuplicates.summary.

```
#### 4.2 Usage

```
name="IMR90_rep1"
resultsDir="/projects/IMR90"
resolution=40000
bin="/projects/bin"
validP="${resultsDir}/s3/w${resolution}/${name}.validPairs"
validI="${resultsDir}/s4/w${resolution}/${name}.validPairs"
mappFile="${bin}/human-hg19.HindIII.w${resolution}"
minMap=0.5 #min mappability threshold
minCount=1 #min contact counts allowed

bash s4_bin.sh "$validP" "$validI" "$bin" "$mappFile" "$minMap" "$minCount" "$maxIter" "$resultsDir/mHiC.summary_w${resolution}_s4"
```

#### 4.3 Input file - mappability file

```
chr1    20000   10      0.003775        0.526966666667
chr1    60000   8       0.01415 0.358675
chr1    100000  19      0.013   0.3898
chr1    140000  12      0.0008  0.4595
chr1    180000  4       0.0008  0.450192340816
chr1    220000  7       0.000575        0.402765636176
chr1    260000  13      0.004275        0.386305422274
chr1    300000  1       0.0     0.450679526523
chr1    340000  16      0.00055 0.421625
chr1    380000  10      0.000225        0.3797
chr1    420000  14      0.000575        0.44275
chr1    460000  10      0.0007  0.467705942362
chr1    500000  0       0.0     0
chr1    540000  13      0.007025        0.501630772417
chr1    580000  13      0.01005 0.39495
chr1    620000  10      0.0029  0.38705
```

### Step 5 - mHi-C generative model prior building [s5_prior.py]
Build the prior for mHi-C model using uni-reads only.

#### 5.0 Requirements
- python3 (>= 3.6)
- numpy (>= 1.13.1)
- scipy (>= 0.19.1)
- sklearn (>= 0.19.1) 

#### 5.1 Arguments

```
fragment	(-f/--fragments)	    : Marginal bin list of all the interactions. Midpoints of each fragment is utilized.
interaction	(-i/--interactions)	: Interaction files of uni-reads normalized by ICE in step 4.
outdir		(-o/--outdir)		      : Path to save outputs.
splineBin	(-b/--noOfBins)	    	: Number of equal-occupancy bins. Default is 100.
priorName	(-l/--lib)		        : Name of file that save the prior quantifying the relationship between random contact probability and genomic distances	.
```

#### 5.2 Usage

```
name="IMR90_rep1"
resultsDir="/projects/IMR90"
resolution=40000
validI="${resultsDir}/s4/w${resolution}/${name}.validPairs"
splineBin=200
priorName="uniPrior"

python s5_prior.py -f $validI.binPair.Marginal -i $validI.binPairCount.uni.afterICE -o ${resultsDir}/s5 -b $splineBin -l $priorName
```

### Step 6 - mHi-C assigning multi-reads [s6_em.py]
In this step, allocation probabilities are assigned to each multi-mapping reads at each potential alignment position. s6_em_cython.pyx will be called to accelerate computation process.

#### 6.0 Requirements
- python3 (>= 3.6)
- pyximport

#### 6.1 Arguments

```
prior		(-p/--prior)        : Prior built in step 5.
uniCount	(-u/--uni)		    : Uni-reads bin-pair contact count file.
multiBinPair	(-m/--multi)	: Multi-reads contact file.
multiKeys	(-mk/--multikeys)	: Multi-reads contact bin keys.
threshold	(-t/--threshold)	: Multi-reads contact probability threshold to extact high quality multi-mapping contact.
filename	(-f/--filename)		: Multi-reads contact posterior probability assignment output file name. ".mHiC" will be added as suffix.
outdir		(-o/--outdir)		  : Output directory to save results.
verbose		(-v/--verbose)		: [Optional] Verbose. Default is true.
```

#### 6.2 Usage

```
name="IMR90_rep1"
resultsDir="/projects/IMR90"
resolution=40000
prior="${resultsDir}/s5/s5_w${resolution}_splineResults"
multi="${resultsDir}/s4/${name}.validPairs.MULTI.binPair.multi"
multiKeys="$resultsDir/s4/${name}.validPairs.MULTI.binPair.multiKeys" 
uni="$resultsDir/s4/${name}.validPairs.binPairCount.uni"
filename="${name}.validPairs.binPair.multi"
threshold=0.5

awk -v OFS="_" '{print $2, $3, $4, $5}' $multi | sort -u >$multiKeys

python s6_em.py -p $prior -u $uni -m $multi -mk $multiKeys -t $threshold -o "${resultsDir}/s6" -f $filename
```

#### 6.3 Input file - Uni-reads bin-pair contact count file
```
chr1    20000   chr3    196460000       1
chr1    20000   chrX    155260000       1
chr1    60000   chr2    199340000       1
chr1    60000   chr4    84780000        1
chr1    60000   chr4    164380000       1
chr1    60000   chr7    2340000 1
chr1    60000   chr7    11700000        2
chr1    60000   chr7    29940000        2
chr1    60000   chr7    34460000        1
chr1    60000   chr7    45340000        1
chr1    60000   chr8    98140000        1
chr1    60000   chr8    131860000       1
chr1    60000   chr11   107220000       1
chr1    60000   chr13   34220000        1
chr1    60000   chr13   102580000       1
```

Multi-reads contact file
```
HWI-ST216_0180:3:1204:18475:195018      chr10   60000   chrX    42260000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr11   96940000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr11   96940000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr11   97620000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr11   97620000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr12   103860000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr12   103860000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr12   48780000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr12   48780000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr12   74860000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr12   74860000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr18   23340000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr18   23340000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr18   20000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr18   20000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr19   23780000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr19   23780000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chr20   19260000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chr20   19260000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chrX    104900000
HWI-ST216_0283:7:2206:17601:52440       chr10   60000   chrX    104900000
HWI-ST216_0283:7:2102:4063:178540       chr10   60000   chrX    115100000
```

### Step 7 - Significant contacts detection [Fit-Hi-C]

Fit-Hi-C pipeline: https://github.com/ay-lab/fithic

### bin/
Other necessary scripts and supplementary data can be found there. 

* cutsite_trimming_mHiC.cpp : utilized to trim unmapped reads to save chimeric reads.
* ICE-with-sparseMatrix.py : scripts for ICE normalization.
* human-hg19.HindIII.w300000, human-hg19.HindIII.w40000, hg19.MboI.w5000 : mappability file for different restriction enzyme and resolution.
* HindIII_resfrag_hg19.bed : restriction enzyme file.

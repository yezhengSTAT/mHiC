# mHiC
Multi-mapping strategy for Hi-C data in order to make use of reads aligned to multiple positions.

## mhic_step0-6.sh
Caller for all the steps in mHi-C pipeline, starting from alignment to multi-reads alignment probability assignment. Customize it for you own use.

## s1_bwaAlignment.sh
Step 1 BWA alignment. Other aligner can also be used instead of BWA. Users need to write their own code to align Hi-C fastq files. Besides, chimeric reads are rescued in this step.

## s2_joinEnd.py 
Step 2 join the read ends by read ID to form a paired-end alignment file.

## s3_categorizePairs.py
Step 3 cateogrize aligned read pairs and pass on the valid read pairs for downsteam analysis.

## s4_bin.sh
Step 4 duplicates removal and bin the genome.

## s5_prior.py
Step 5 build the prior for mHi-C model using uni-reads only. myStats.py and myUtils.py will be called.

## s6_em.py
Step 6 mHi-C assign probability to all potential alignment positions for each multi-reads. s6_em_cython.pyx will be called and accelerate computation process.

## bin/
Other necessary scripts can be found there. 
cutsite_trimming_mHiC.cpp is utilized to trim unmapped reads to save chimeric reads.
ICE-with-sparseMatrix.py is the scripts for ICE normalization.
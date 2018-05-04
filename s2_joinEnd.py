#!/usr/bin/env python

## mHiC
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu

'''
Script to join two single-end alignment with multi-mapping reads into pair-ended bam file.
Only aligned pair-end reads are kept.
Alignment statistics are summarized by demand.
Dec, 2016
'''

import sys
import os
import re
import pysam
import argparse


def get_args():
    """Get arguments"""
    parser = argparse.ArgumentParser(description = '------------Usage Start------------',
                                     epilog = '------------Usage End------------')
    parser.add_argument('-r1', '--readEnd1', help = 'Path to the aligned file of read end 1.')
    parser.add_argument('-r2', '--readEnd2', help = 'Path to the aligned file of read end 2.')
    parser.add_argument('-o', '--output', help = 'Output file.', default = None)
    parser.add_argument('-s', '--summary', help = '(Optional) Summary of alignment results. Default is true.', default = True)
    parser.add_argument('-sf', '--summaryFile', help = '(Optional) Summary file path . Default is infileName.alignSummary.', default = None)
    parser.add_argument('-v', '--verbose', help = '(Optional) Verbose. Default is true.', default = True)

    args = parser.parse_args()
    if args.readEnd1 is None or args.readEnd2 is None or args.output is None:
        parser.print_help()
        sys.exit()
    if not os.path.exists(args.readEnd1):
        print("Single-ended read 1 file does not exist!")
        sys.exit()
    if not os.path.exists(args.readEnd2):
        print("Single-ended read 2 file does not exist!")
        sys.exit()
    if args.output is None:
        print("Please specify the output file to store the pair-ended alignment files after join the ends.")
        sys.exit()
    if not os.path.exists(args.output.rsplit("/", 1)[0]):
        os.makedirs(args.output.rsplit("/", 1)[0])
    return args

def is_uniread(read):
    if not read.is_unmapped and read.mapping_quality > 24 and not read.has_tag('XA'):
        return True
    else:
        return False

def is_lowQualMulti(read):
    if not read.is_unmapped and read.mapping_quality < 25 and not read.has_tag("XA"):
        return True
    else:
        return False
    
def sam_flag(r1, r2, file1, file2):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mergeSAM.py) by Nicolas Servant, Eric Viara
    '''

    f1 = r1.flag
    f2 = r2.flag

    if r1.is_unmapped == False:
        r1_chrom = file1.getrname(r1.tid)
    else:
        r1_chrom="*"
    if r2.is_unmapped == False:
        r2_chrom = file2.getrname(r2.tid)
    else:
        r2_chrom="*"


  ##Relevant bitwise flags (flag in an 11-bit binary number)
  ##1 The read is one of a pair
  ##2 The alignment is one end of a proper paired-end alignment
  ##4 The read has no reported alignments
  ##8 The read is one of a pair and has no reported alignments
  ##16 The alignment is to the reverse reference strand
  ##32 The other mate in the paired-end alignment is aligned to the reverse reference strand
  ##64 The read is the first (#1) mate in a pair
  ##128 The read is the second (#2) mate in a pair
  
  ##The reads were mapped as single-end data, so should expect flags of 
  ##0 (map to the '+' strand) or 16 (map to the '-' strand)
  ##Output example: a paired-end read that aligns to the reverse strand 
  ##and is the first mate in the pair will have flag 83 (= 64 + 16 + 2 + 1)
  
    if f1 & 0x4:
        f1 = f1 | 0x8

    if f2 & 0x4:
        f2 = f2 | 0x8
    
    if (not (f1 & 0x4) and not (f2 & 0x4)):
    ##The flag should now indicate this is paired-end data
        f1 = f1 | 0x1
        f1 = f1 | 0x2
        f2 = f2 | 0x1
        f2 = f2 | 0x2
  
    
  ##Indicate if the pair is on the reverse strand
    if f1 & 0x10:
        f2 = f2 | 0x20
  
    if f2 & 0x10:
        f1 = f1 | 0x20
  
  ##Is this first or the second pair?
    f1 = f1 | 0x40
    f2 = f2 | 0x80
  
    ##Insert the modified bitwise flags into the reads
    r1.flag = f1
    r2.flag = f2

    ##Determine the RNEXT and PNEXT values (i.e. the positional values of a read's pair)
    #RNEXT
    if r1_chrom == r2_chrom:
        r1.rnext = r1.tid
        r2.rnext = r1.tid
    else:
        r1.rnext = r2.tid
        r2.rnext = r1.tid
   
   #PNEXT
    r1.pnext = r2.pos
    r2.pnext = r1.pos
 
    return(r1, r2)


if __name__ == "__main__":

    args = get_args()
    if args.verbose:
        print("Read 1 file = ", args.readEnd1)
        print("Read 2 file = ", args.readEnd2)
        print("Output file = ", args.output)
        print("Summary of alignment results = ", args.summary)
        print("Summary file path = ", args.summaryFile)
        print("Verbose = ", args.verbose)

    ## Initial count variable for summary
    singleton_count = 0
    unmapped_count = 0
    mapped_count = 0
    uni_count = 0
    multi_count = 0
    multi_lowqual_count = 0
    total_count = 0

    if args.verbose:
        print("------------Begin joining ends from read 1 file and read 2 file------------")
    with pysam.Samfile(args.readEnd1, "rb") as readEnd1, pysam.Samfile(args.readEnd2, "rb") as readEnd2:

        if args.output.endswith(".bam"):
            outfile = pysam.AlignmentFile(args.output, "wb", template = readEnd1)
        elif args.output.endswith(".sam"):
            outfile = pysam.AlignmentFile(args.output, "wh", template = readEnd1)
        else:
            print("Output file format should be either bam or sam.")
            sys.exit()
            
        for r1, r2 in zip(readEnd1.fetch(until_eof = True), readEnd2.fetch(until_eof = True)):
            total_count += 1

            if total_count % 1000000 == 0 and args.verbose:
                print("# of reads processed: ", total_count)

            if r1.qname == r2.qname:
                if r1.is_unmapped == True and r2.is_unmapped == True:
                    unmapped_count += 1
                    continue
                elif r1.is_unmapped == False and r2.is_unmapped == False:
                    mapped_count += 1
                    if is_uniread(r1) == True and is_uniread(r2) == True:
                        uni_count += 1
                        (r1, r2) = sam_flag(r1, r2, readEnd1, readEnd2)
                        outfile.write(r1)
                        outfile.write(r2)
                    elif is_lowQualMulti(r1) == True or is_lowQualMulti(r2) == True:
                        """
                        Need to check if it is low quality multi-mapping read, i.e. if score = 0 but no XA tag, namely multi-read which can be aligned to too many positions beyound upper bound and they are considered as low quality multi-mapping read -  they will be discarded eventually.
                        """
                        multi_lowqual_count += 1
                    else:
                        """
                        the rest aligned reads pair with at least one end having XA tag are multi-mapping read pairs.
                        """
                        multi_count += 1
                        (r1, r2) = sam_flag(r1, r2, readEnd1, readEnd2)
                        outfile.write(r1)
                        outfile.write(r2)

                else:
                    singleton_count += 1
            else:
                print("Read id in read 1 and read 2 file does not match. Please check the input read files and sort them correctly.")
                sys.exit()

    if args.summary:
        if args.summaryFile is None:
            summaryFile = re.sub("\.bam$|\.sam$", ".alignSummary", args.output)
            statSummary = open(summaryFile, "w+")

        else:
            summaryFile = args.summaryFile
            statSummary = open(summaryFile, "a")

        statSummary.write("Step: Alignment and join both ends." + "\n" + "\n")
        statSummary.write("Total number of read pairs:\t" + str(total_count) + "\n")
        statSummary.write("  Total number of unmapped read pairs:\t" + str(unmapped_count) + "\n")
        statSummary.write("  Total number of singleton read pairs:\t" + str(singleton_count) + "\n")
        statSummary.write("  Total number of mapped read pairs:\t" + str(mapped_count) + "\n")
        statSummary.write("    Total number of uniquely mapping read pair for both ends:\t" + str(uni_count) + "\n")
        statSummary.write("    Total number of multi-mapping for at least one end of the read pairs:\t" + str(multi_count + multi_lowqual_count) + "\n")
        statSummary.write("    Total number of multi-mapping for at least one end of the read pairs and is of high quality:\t" + str(multi_count) + "\n")
        statSummary.write("    Total number of multi-mapping for at least one end of the read pairs and is of low quality:\t" + str(multi_lowqual_count) + "\n")
        statSummary.close()

    readEnd1.close()
    readEnd2.close()
    outfile.close()
        
        
        
            
                

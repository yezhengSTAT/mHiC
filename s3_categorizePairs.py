#!/usr/bin/env python

## mHiC
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu

'''
Script to identify Restriction Enzyme fragment where the each read end comes from.
Categorize each read pairs: SC, DE, Religation, Dump, valid read pairs.
Summarize categories.
Only valid read pairs with fragment assignment are saved in output.
Dec, 2016
'''

import sys
import os
import re
import pysam
import argparse
from bx.intervals.intersection import Intersecter, Interval

def get_args():
    '''Get arguments'''
    parser = argparse.ArgumentParser(description = '------------Usage Start------------',
                                     epilog = '------------Usage End------------')
    parser.add_argument('-f', '--fragment', help = 'Restriction Enzyme(RE) fragment file.', default = None)
    parser.add_argument('-r', '--read', help = 'Paired-ended read alignment file.', default = None)
    parser.add_argument('-o', '--outdir', help = 'Output directory to save category results.', default = None)
    parser.add_argument('-l', '--lower', help = 'Lower bound of the read pair distances summation from read end alignment position to its assigned RE fragment cutting site. Default value is None indicating no restriction.', default = None)
    parser.add_argument('-u', '--upper', help = 'Upper bound of the read pair distances summation from read end alignment position to its assigned RE fragment cutting site. Default value is None indicating no restriction.', default = None)
    parser.add_argument('-d', '--distance', help = 'Minimum distance between intrachromosomal read pair alignments positions. Default value is None indicating no restriction.', default = None)
    parser.add_argument('-m', '--method', help = "Binning methods: by fixed number of window size, i.e. window, or fixed number of RE fragment, i.e. fragment.", default = "fragment")
    parser.add_argument('-b', '--bin', help = "Binning resolution. If binning by fixed number of window size, it can be 10000 or any other suitable bin size. If binning by fixed number of RE fragment, it can be 10 representing 10RE fragments or other suitable bin size.", default = None)
    parser.add_argument('-s', '--summary', help = '(Optional) Summary of categories. Default is true.', default = True)
    parser.add_argument('-sf', '--summaryFile', help = '(Optional) Summary file path. Default is alignedFileName.readCategorySummary.', default = None)
    parser.add_argument('-v', '--verbose', help = '(Optional) Verbose. Default is true.', default = True)

    args = parser.parse_args()
    if args.fragment is None or args.read is None or args.outdir is None:
        parser.print_help()
        sys.exit()
    if not os.path.exists(args.fragment):
        print("Resctriction Enzyme fragment file does not exist!")
        sys.exit()
    if not os.path.exists(args.read):
        print("Paired-ended read alignmeng file does not exist!")
        sys.exit()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    methodList = ['fragment', 'window', 'frag', 'win', 'f', 'w']
    if not any(args.method in x for x in methodList):
        print("Binning method should either be fragment or window.")
        sys.exit()
    if args.bin is None:
        print("Please provide number representing binning resolution.")
        sys.exit()
    else:
        args.bin = int(args.bin)
    return args


def get_read_strand(read):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''


    """
    Conversion of read position to naive strand representation
    Parameters
    ----------
    read : list
        list of aligned reads
    """
    strand = "+"
    if read.is_reverse:
        strand = "-"
    return strand


def isIntraChrom(read1, read2):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''
    
    """
    Return true is the reads pair is intrachromosomal
    
    read1 : [AlignedRead]
    read2 : [AlignedRead]
    """
    if read1.tid == read2.tid:
        return True
    else:
        return False


def get_read_pos(read):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''

    """
    Return the read position (zero-based) used for the intersection with
    the restriction fragment
    The 5' end is not a good choice for the reverse reads (which contain part
    of the restriction site, and thus overlap the next restriction fragment)
    Using the left-most position (5' for forward, 3' for reverse) or the
    middle of the read should work but the middle of the reads might be more
    safe
    Parameters
    -----------
    read : list
        list of aligned reads
    """
    if read.alen is not None:
        pos = read.pos + read.alen/2
    elif read.cigarstring is not None:
        alignedLen = re.findall(r'(\d+)[M, D]', read.cigarstring)
        len = sum([int(i) for i in alignedLen])
        pos = read.pos + len/2
    else:
        print("Warning: neither read aligned length nor cigar string!")
        sys.exit()

    return pos


def get_cis_dist(read1, read2):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''
    """
    Calculte the contact distance between two intrachromosomal reads
    read1 : [AlignedRead]
    read2 : [AlignedRead]
    """
    # Get oriented reads
    ##r1, r2 = get_ordered_reads(read1, read2)
    dist = None
    if not r1.is_unmapped and not r2.is_unmapped:         
        ## Contact distances can be calculated for intrachromosomal reads only
        if isIntraChrom(read1, read2):
            r1pos = get_read_pos(read1)
            r2pos = get_read_pos(read2)
            dist = abs(r1pos - r2pos)
            return dist
        
        
def get_read_start(read):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''
    
    """
    Return the 5' end of the read
    """
    if read.is_reverse:
        if read.alen is not None:
            pos = int(read.pos + read.alen)
        elif read.cigarstring is not None:
            alignedLen = re.findall(r'(\d+)[M, D]', read.cigarstring)
            len = sum([int(i) for i in alignedLen])
            pos = int(read.pos + len)
        else:
            print("Warning: neither read aligned length nor cigar string!")
            sys.exit()
            
    else:
        pos = read.pos
    return pos

def get_ordered_reads(read1, read2):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''
    """
    Reorient reads
    The sequencing is usually not oriented. Reorient the reads so that r1 is
    always before r2.
   
    read1 = [AlignedRead]
    read2 = [AlignedRead]
    """
    if read1.tid == read2.tid:
        if get_read_pos(read1) < get_read_pos(read2):
            r1 = read1
            r2 = read2
        else:
            r1 = read2
            r2 = read1
    else:
        if read1.tid < read2.tid:
            r1 = read1
            r2 = read2
        else:
            r1 = read2
            r2 = read1
                
    return r1, r2


def get_overlapping_restriction_fragment(resFrag, chrom, read):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''
    """
    Intersect a given read with the set of restriction fragments
    ##
    resFrag = the restriction fragments [hash]
    chrom = the chromosome to look at [character]
    read = the read to intersect [AlignedRead]
    """
    # Get read position (middle or 5' end)
    pos = int(get_read_pos(read))
    
    if chrom in resFrag:
        # Overlap with the position of the read (zero-based)
        resfrag = resFrag[chrom].find(pos, pos+1)
        if len(resfrag) > 1:
            print("Warning : " +  len(resfrag) + " restriction fragments found for " + read.qname + "- skipped")
            return None
        elif len(resfrag) == 0:
            print("Warning - no restriction fragments for " + read.qname + " at " + chrom + ":" + str(pos))
            return None
        else:
            return resfrag[0]
    else:
        print("Warning - no restriction fragments for " + read.qname + " at " + chrom + ":" + str(pos))
        return None

def load_restriction_fragment(in_file, minfragsize=None, maxfragsize=None, verbose=False):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''
    '''
    Function load_restriction_fragment cite the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara 
    '''
    
    """
    Read a BED file and store the intervals in a tree
    Intervals are zero-based objects. The output object is a hash table with
    one search tree per chromosome
    in_file = input file [character]
    verbose = verbose mode [logical]
    """
    resFrag = {}
    if verbose:
        print("## Loading Restriction File Intervals '" + in_file + "'...")

    bed_handle = open(in_file)
    nline = 0
    for line in bed_handle:
        nline +=1
        bedtab = line.split("\t")
        try:
            chromosome, start, end, name = bedtab[:4]
        except ValueError:
            print("Warning : wrong input format in line" + nline + ". Not a BED file !?")
            continue

        # BED files are zero-based as Intervals objects
        start = int(start)  # + 1
        end = int(end)
        midPoint = (start + end)/2
        fragl = abs(end - start)
        name = name.strip()

        ## Discard fragments outside the size range
        if minfragsize != None and int(fragl) < int(minfragsize):
            print("Warning : fragment "+ name + " [" +  fragl + "] outside of range. Discarded")
            continue
        if maxfragsize != None and int(fragl) > int(maxfragsize):
            print("Warning : fragment " + name + " [" + fragl + "] outside of range. Discarded")
            continue
       
        if chromosome in resFrag:
            tree = resFrag[chromosome]
            tree.add_interval(Interval(start, end, value={'name': name, 'midPoint': midPoint}))
        else:
            tree = Intersecter()
            tree.add_interval(Interval(start, end, value={'name': name, 'midPoint': midPoint}))
            resFrag[chromosome] = tree
    
    bed_handle.close()
    return resFrag



def get_overlapping_resFragBin(binTree, chrom, read):
    """
    Intersect a given read with the set of bin
    ##
    binTree = the bin tree [hash]
    chrom = the chromosome to look at [character]
    read = the read to intersect [AlignedRead]
    """
    # Get read position (middle or 5' end)
    pos = int(get_read_pos(read))
    
    if chrom in binTree.keys():
        # Overlap with the position of the read (zero-based)
        resFragBin = binTree[chrom].find(pos, pos+1)
        if len(resFragBin) > 1:
            print("Warning : " +  str(len(resFragBin)) + " bins found for " + read.qname + "- skipped")
            print(pos)
            print(resFragBin)
            return None
        elif len(resFragBin) == 0:
            print("Warning - no bin found for " + read.qname + " at " + chrom + ":" + str(pos))
            return None
        else:
            return resFragBin[0]
    else:
        print("Warning - no bin found for " + read.qname + " at " + chrom + ":" + str(pos))
        return None

def build_resFragBin_tree(in_file, resolution = None, verbose = False):
    """
    build restriction fragments bins tree based on the user defined resolution, e.g. 10fragment
    """
    binTree = {}
    if resolution is None:
        print("Please provide suitable resolution to bin contact! For example resolution = 10 for 10 RE fragment as one bin.")
        sys.exit()
    else:
        resolution = int(resolution)
        
    if verbose:
        print("## Building RE fragment bins tree from ordered'" + in_file + "'...")
    bed_handle = open(in_file)
    nline = 0
    flag = 0
    for line in bed_handle:
        nline += 1

        bedtab = line.split("\t")
        try:
            chromosome, start, end = bedtab[:3]
            start = int(start)
            end = int(end)
        except ValueError:
            print("Warning : wrong input format in line" + nline + ". Not a BED file !?")
            continue
            
        # First line
        if flag == 0:
            #Bin variables store current bin information
            startBin = start
            chromosomeBin = chromosome
            flag = 1
        elif chromosomePre != chromosome: #Start of another chromosome
            nline = 1
            if startBin != endPre:
                endBin = endPre
                midBin = round((startBin + endBin)/2)
                   
                if chromosomeBin in binTree.keys():
                    tree = binTree[chromosomeBin]
                    tree.add_interval(Interval(startBin, endBin, value={'midPoint': midBin}))
                else:
                    tree = Intersecter()
                    tree.add_interval(Interval(startBin, endBin, value={'midPoint': midBin}))
                    binTree[chromosomeBin] = tree
            startBin = start
            chromosomeBin = chromosome
        elif nline % resolution == 0:
            endBin = end
            midBin = round((startBin + endBin)/2)
                   
            if chromosomeBin in binTree.keys():
                tree = binTree[chromosomeBin]
                tree.add_interval(Interval(startBin, endBin, value={'midPoint': midBin}))
            else:
                tree = Intersecter()
                tree.add_interval(Interval(startBin, endBin, value={'midPoint': midBin}))
                binTree[chromosomeBin] = tree

            startBin = end
            chromosomeBin = chromosome

        # Update the Pre variable set
        startPre = start
        endPre = end
        chromosomePre = chromosome
        
    # for the last bin
    if nline % resolution != 0 and startBin != start and chromosomeBin == chromosome:
        endBin = end
        midBin = round((startBin + endBin)/2)
                   
        if chromosomeBin in binTree.keys():
            tree = binTree[chromosomeBin]
            tree.add_interval(Interval(startBin, endBin, value={'midPoint': midBin}))
        else:
            tree = Intersecter()
            tree.add_interval(Interval(startBin, endBin, value={'midPoint': midBin}))
            binTree[chromosomeBin] = tree

    bed_handle.close()
    return binTree




            

def are_contiguous_fragments(frag1, frag2, chr1, chr2):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''
    '''
    Compare fragment positions to check if they are contiguous
    '''
    ret = False
    if chr1 == chr2:
        if int(frag1.start) < int(frag2.start):
            d = int(frag2.start) - int(frag1.end)
        else:
            d = int(frag1.start) - int(frag2.end)
        
        if d == 0:
            ret = True
            
            return ret

def is_religation(read1, read2, frag1, frag2):
    '''
    This function cites the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara
    '''
    """
    Reads are expected to map adjacent fragments
    Check the orientation of reads -><-
    """
    ret=False
    if are_contiguous_fragments(frag1, frag2, read1.tid, read2.tid):
        ret=True
    return ret


def is_self_circle(read1, read2):
    '''
    Function cite the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara 
    '''
    """
    Both reads are expected to be on the same restriction fragments
    Check the orientation of reads <-->
    read1 : [AlignedRead]
    read2 : [AlignedRead]
    """
    ret = False
    # Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)
    # 1<- ->2 or 2<- ->1
    if get_read_strand(r1) == "-" and get_read_strand(r2) == "+":
        ret = True
    return ret


def is_dangling_end(read1, read2):
    '''
    Function cite the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara 
    '''
    """
    Both reads are expected to be on the same restriction fragments
    Check the orientation of reads -><-
    read1 : [AlignedRead]
    read2 : [AlignedRead]
    """
    ret = False
    # Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)
    # 1-> <-2 or 2-> <-1
    if get_read_strand(r1) == "+" and get_read_strand(r2) == "-":
        ret = True
    return ret


def get_interaction_type(read1, read1_chrom, resfrag1, read2, read2_chrom, resfrag2, verbose):
    '''
    Function cite the same function from Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara 
    '''
    """
    Returns the interaction type
    For a given reads pair and their related restriction fragment, classify
    the 3C products as :
    - Singleton Interaction (SI)
    - Valid Interaction (VI)
    - Self circle (SC)
    - Dangling end (DE)
    - Religation (RE)
    - Unknown (DUMP)

    ##
    read1 = the R1 read of the pair [AlignedRead]
    read1_chrom = the chromosome of R1 read [character]
    resfrag1 = restrictin fragment overlapping the R1 read [interval]
    read2 = the R2 read of the pair [AlignedRead]
    read2_chrom = the chromosome of R2 read [character]
    resfrag2 = restrictin fragment overlapping the R2 read [interval]
    verbose = verbose mode [logical]
    """

    # If returned InteractionType=None -> Same restriction fragment
    # and same strand = Dump
    interactionType = None
 
    if not r1.is_unmapped and not r2.is_unmapped and resfrag1 is not None and resfrag2 is not None:
        # same restriction fragment
        if resfrag1 == resfrag2:
            # Self_circle <- ->
            if is_self_circle(read1, read2):
                interactionType = "SC"
            # Dangling_end -> <-
            elif is_dangling_end(read1, read2):
                interactionType = "DE"
        elif is_religation(read1, read2, resfrag1, resfrag2):
            interactionType = "RE"
        else:
            interactionType = "VI"
    elif r1.is_unmapped or r2.is_unmapped:
        interactionType = "SI"

    return interactionType

def get_PE_fragment_size(read1, read2, resFrag1, resFrag2, interactionType):
    '''
    Function modified from function of Hi-C Pro (https://github.com/nservant/HiC-Pro/blob/master/scripts/mapped_2hic_fragments.py) by Nicolas Servant, Eric Viara 
    '''
    """
    Calculte the size of the DNA fragment library - (Ye: distances from read end alignment to its fragment cutting site.)
    read1 : [AlignedRead]
    read2 : [AlignedRead]
    resfrag1 = restrictin fragment overlapping the R1 read [interval]
    resfrag2 = restrictin fragment overlapping the R2 read [interval]
    interactionType : Type of interaction from get_interaction_type() [str]
    """

    fragmentsize = None

    # Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)
    if not r1.is_unmapped and not r2.is_unmapped:
        if r1 == read2:
            rfrag1 = resFrag2
            rfrag2 = resFrag1
        else:
            rfrag1 = resFrag1
            rfrag2 = resFrag2

        ## In this case use the read 5' end !
        r1pos = get_read_start(r1)
        r2pos = get_read_start(r2)

        if interactionType == "DE" or interactionType == "RE":
            fragmentsize = r2pos - r1pos
        elif interactionType == "SC":
            fragmentsize = (r1pos - rfrag1.start) + (rfrag2.end - r2pos)
        elif interactionType == "VI":
            dr1 = min(abs(rfrag1.end - r1pos), abs(r1pos - rfrag1.start))
            dr2 = min(abs(rfrag2.end - r2pos), abs(r2pos - rfrag2.start))
            # if get_read_strand(r1) == "+":
            #     dr1 = rfrag1.end - r1pos
            # else:
            #     dr1 = r1pos - rfrag1.start
            # if get_read_strand(r2) == "+":
            #     dr2 = rfrag2.end - r2pos
            # else:
            #     dr2 = r2pos - rfrag2.start
            fragmentsize = dr2 + dr1

    return fragmentsize

def get_valid_pair(r1, r2, r1_chrom, r2_chrom, r1_frag, r2_frag, r1_bin, r2_bin, cisDist):
    # print valid read pairs out
    or1, or2 = get_ordered_reads(r1, r2)

    if or1 == r1 and or2 == r2:
        or1_chrom = r1_chrom
        or2_chrom = r2_chrom
        or1_fragname = r1_frag.value['name'] if r1_frag is not None else "NA"
        or2_fragname = r2_frag.value['name'] if r2_frag is not None else "NA"
        # or1_fragPos = r1_frag.value['midPoint'] if r1_frag is not None else 0
        # or2_fragPos = r2_frag.value['midPoint'] if r2_frag is not None else 0
        or1_bin = r1_bin
        or2_bin = r2_bin
    elif or1 == r2 and or2 == r1:
        or1_chrom = r2_chrom
        or2_chrom = r1_chrom
        or1_fragname = r2_frag.value['name'] if r2_frag is not None else "NA"
        or2_fragname = r1_frag.value['name'] if r1_frag is not None else "NA"
        # or1_fragPos = r2_frag.value['midPoint'] if r2_frag is not None else 0
        # or2_fragPos = r1_frag.value['midPoint'] if r1_frag is not None else 0
        or1_bin = r2_bin
        or2_bin = r1_bin

    validInfo = "\t".join([or1.qname, or1_chrom, str(get_read_pos(or1)), get_read_strand(or1), or1_fragname, str(or1_bin),\
                           or2_chrom, str(get_read_pos(or2)), get_read_strand(or2), or2_fragname, str(or2_bin),\
                           "InterChrom" if cisDist is None else str(cisDist)])
    return validInfo

def get_multi_read(readFile, readName, tag):
    read = pysam.AlignedSegment()

    read.qname = readName
    readInfo = tag.split(',')
    read.tid = readFile.get_tid(readInfo[0])
    read.is_reverse = True if readInfo[1].startswith('-') else False
    read.pos = int(readInfo[1].strip('-')) - 1
    read.cigarstring = readInfo[2].strip()
    return read
    


if __name__ == "__main__":
    args = get_args()

    if args.verbose:
        print("RE fragment = ", args.fragment)
        print("Paired-ended read alignment file = ", args.read)
        print("Output directory = ", args.outdir)
        print("RE cutting site distance lower bound = ", args.lower)
        print("RE cutting site distance upper bound = ", args.upper)
        print("Minimum intrachromosomal contact distance = ", args.distance)
        print("Binning method = ", args.method)
        print("Binning resolution = ", args.bin)
        print("Summary = ", args.summary)
        print("SummaryFile = ", args.summaryFile)
        print("Verbose = ", args.verbose)

    # count for uni-mapping read pairs
    sc_count = 0
    de_count = 0
    re_count = 0
    ud_count = 0 #uni dump
    si_count = 0
    uv_count = 0
    uni_count = 0
    # count for multi-mapping read pairs/multi-mapping unique valid read pairs/multi-mapping multiple valid read pairs
    multi_count = 0
    md_count = 0 #multi dump
    multi_uv_count = 0
    multi_mv_count = 0
    multi_mv_uf_count = 0
    multi_mv_ub_count = 0
    
    fileName = re.sub("\.bam$|\.sam$", "", os.path.basename(args.read))
    outValid = open(args.outdir + "/" + fileName + ".validPairs", "w")
    outContactNum = open(args.outdir + "/" + fileName + ".uniMulti", "w")
    
    fragTree = load_restriction_fragment(in_file = args.fragment, verbose = args.verbose)
    print(fragTree.keys())

    fragMethodList = ['fragment', 'frag', 'f']
    if any(args.method in x for x in fragMethodList):
        binTree = build_resFragBin_tree(in_file = args.fragment, resolution = args.bin, verbose = args.verbose)
    else:
        binTree = None

    fragMethodList = ['window', 'win', 'w', 'windows']
    if any(args.method in x for x in fragMethodList):
        binWindow = int(args.bin)
    else:
        binWindow = None
        
    if args.verbose:
        print("------------Begin reading paired-ended file------------")

    if args.read.endswith(".bam"):
        readFile = pysam.Samfile(args.read, "rb")
    elif args.read.endswith(".sam"):
        readFile = pysam.Samfile(args.read, "r")
    else:
        print("Input paired-ended alignment file should either be in bam format or sam format.")
        sys.exit()

    if args.verbose:
        print("------------Categorize aligned read pairs------------")

    # new one with all kinds of potential contacts
    for read in readFile.fetch(until_eof = True):
        if read.is_read1:
            if not read.is_unmapped:
                r1 = read
            else:
                r1 = None
        elif read.is_read2:
            if not read.is_unmapped:
                r2 = read
            else:
                r2 = None

            if r1 is not None and r2 is not None:
                m1 = [r1]
                if r1.has_tag('XA'):
                    for xaTag in r1.get_tag('XA').strip(';').split(';'):
                        m1.append(get_multi_read(readFile, r1.qname, xaTag))
                m2 = [r2]
                if r2.has_tag('XA'):
                    for xaTag in r2.get_tag('XA').strip(';').split(';'):
                        m2.append(get_multi_read(readFile, r2.qname, xaTag))
                valid_count = 0
                fragPair = []
                binPair = []
                validInfoList = []
                for end1 in m1:
                    for end2 in m2:
                        
                        end1_chrom = readFile.getrname(end1.tid)
                        end2_chrom = readFile.getrname(end2.tid)
                        end1_frag = get_overlapping_restriction_fragment(fragTree, end1_chrom, end1)
                        if end1_frag is None:
                            interactionType = "DUMP"
                            break

                        end2_frag = get_overlapping_restriction_fragment(fragTree, end2_chrom, end2)
                        if end2_frag is None:
                            interactionType = "DUMP"
                            break
                        # get interaction type, alignmenet position to fragment cutting site distance sum, and intrachromosal contact distance
                        if end1_frag is not None and end2_frag is not None:
                            interactionType = get_interaction_type(end1, end1_chrom, end1_frag, end2, end2_chrom, end2_frag, args.verbose)
                            fragDist = get_PE_fragment_size(end1, end2, end1_frag, end2_frag, interactionType)
                            cisDist = get_cis_dist(end1, end2)
                            # If the alignment position to fragment cutting sites distance sum is out of range allowed
                            if (args.lower is not None and fragDist is not None and fragDist < int(args.lower)) or \
                               (args.upper is not None and fragDist is not None and fragDist > int(args.upper)):
                                interactionType = "DUMP"


                            # If intrachromosal contact distance is too near
                            if interactionType == "VI" and args.distance is not None and cisDist is not None and cisDist < int(args.distance):
                                interactionType = "DUMP"

                            if interactionType == "VI": #if this potential contact is valid - write out
                                valid_count += 1
                                # get the overlapped bin info
                                if binTree != None:
                                    end1_bin = get_overlapping_resFragBin(binTree, end1_chrom, end1)
                                    end2_bin = get_overlapping_resFragBin(binTree, end2_chrom, end2)
                                    mid1 = end1_bin.value["midPoint"]
                                    mid2 = end2_bin.value["midPoint"]
                                elif binWindow != None:
                                    mid1 = int(int(get_read_pos(end1))/binWindow) * binWindow + int(binWindow/2)
                                    mid2 = int(int(get_read_pos(end2))/binWindow) * binWindow + int(binWindow/2)
                                else:
                                    print("No binning method specified!")
                                    sys.exit()

                                if  end1_chrom != end2_chrom or mid1 != mid2:
                                    # Summarize valid pair information and prepare to write out
                                    validInfo = get_valid_pair(end1, end2, end1_chrom, end2_chrom, end1_frag, end2_frag, mid1, mid2, cisDist)
                                    validInfoList.append(validInfo)
                                    
                                    # valid interactions cover how many distinct fragment pairs
                                    fragPair.append([end1.qname, end1_frag.value["name"], end2_frag.value["name"]])
                                    # valid interactions cover how many distinct bins
                                    binPair.append([end1.qname, end1_chrom, mid1, end2_chrom, mid2])
                                else:
                                    interactionType = "DUMP"
                            
                # Identify read pair type: Uni or Multi
                if len(validInfoList) == 1:
                    outValid.write(validInfoList[0] + "\t" + "UNI" + "\n")
                elif len(validInfoList) > 1:
                    for eachInfo in validInfoList:
                        outValid.write(eachInfo + "\t" + "MULTI" + "\n")

                # Identify fragment pair type: Uni or Multi
                if len(fragPair) > 0:
                    fragPairUni = [list(x) for x in set(tuple(x) for x in fragPair)]
                else:
                    fragPairUni = []

                # Identify bin pair type: Uni or Multi
                if len(binPair) > 0:
                    binPairUni = [list(x) for x in set(tuple(x) for x in binPair)]
                else:
                    binPairUni = []

                # Count
                if len(m1) == 1 and len(m2) == 1: #uni-read pair
                    uni_count += 1
                    outContactNum.write('\t'.join([r1.qname, "Uni", "1"]) + "\n")
                    if interactionType == "VI":
                        uv_count += 1
                    elif interactionType == "DE":
                        de_count += 1
                    elif interactionType == "SC":
                        sc_count += 1
                    elif interactionType == "RE":
                        re_count += 1
                    else:
                        interactionType = "DUMP"
                        ud_count += 1
                else: #multi-read pair
                    multi_count +=1
                    if valid_count == 1: #indicating multi-mapping read pair has unique valid contact
                        multi_uv_count += 1
                        outContactNum.write('\t'.join([r1.qname, "MultiToUni", "1"]) + "\n")
                    elif valid_count == 0:
                        md_count += 1
                    else:
                        multi_mv_count += 1
                        outContactNum.write('\t'.join([r1.qname, "Multi", str(valid_count)]) + "\n")
                        if len(fragPairUni) == 1:
                            multi_mv_uf_count += 1
                        if len(binPairUni) == 1:
                            multi_mv_ub_count += 1
            else:
                interactionType = "DUMP"
                ud_count += 1
        else:
            print("Warning! Read is either read1 nor read2")

    
    if args.summary:
        if args.summaryFile is None:
            outSummary = open(args.outdir + "/" + fileName + ".readCategorySummary", "w")
        else:
            outSummary = open(args.summaryFile, "a")
        outSummary.write("\n" + "\n")
        outSummary.write("Step: Categorize read pair type.\n" + "\n")
        outSummary.write("Total number of read pairs:\t" + str(uni_count + multi_count) + "\n")
        outSummary.write("  Uniquely mapping read pairs:\t" + str(uni_count) + "\n")
        outSummary.write("    Uniquely mapping valid read pairs:\t" + str(uv_count) + "\n")
        outSummary.write("    Uniquely mapping dangling end read pairs:\t" + str(de_count) + "\n")
        outSummary.write("    Uniquely mapping self cycle read pairs:\t" +str(sc_count) + "\n")
        outSummary.write("    Uniquely mapping religation read pairs:\t" + str(re_count) + "\n")
        outSummary.write("    Uniquely mapping dumped read pairs:\t" + str(ud_count) + "\n")
        outSummary.write("  Multi-mapping read pairs:\t" + str(multi_count) + "\n")
        outSummary.write("    Multi-mapping unique valid read pairs:\t" + str(multi_uv_count) + "\n")
        outSummary.write("    Multi-mapping multiple valid read pairs:\t" + str(multi_mv_count) + "\n")
        outSummary.write("    Multi-mapping multiple valid read pair but unique fragment pairs:\t" + str(multi_mv_uf_count) + "\n")
        outSummary.write("    Multi-mapping multiple valid read pair but unique bin pairs:\t" + str(multi_mv_ub_count) + "\n")
        outSummary.write("    Multi-mapping dumped read pairs:\t" + str(md_count) + "\n")

        outSummary.close()
    outContactNum.close()
    outValid.close()


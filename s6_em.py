## mHi-C core model part
## used c/cython to significant increase the computation efficiency!
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu

'''
Script to assign probability to multi-reads candidates alignment 
positions.
April 2016
'''


import pyximport;pyximport.install()
from s6_em_cython import *
import argparse
import os
import sys
import re

def get_args():
    '''Get arguments'''
    parser = argparse.ArgumentParser(description = '------------Usage Start------------',
                                     epilog = '------------Usage End------------')
    parser.add_argument('-p', '--prior', help = 'Prior built by Fit-Hi-C.', default = None)
    parser.add_argument('-u', '--uni', help = 'Uniquely mapping bin contact count file.', default = None)
    parser.add_argument('-m', '--multi', help = 'Multi-mapping contact bin pair file.', default = None)
    parser.add_argument('-mk', '--multikeys', help = 'Multi-mapping contact bin pair keys.', default = None)
    parser.add_argument('-t', '--threshold', help = 'Multi-mapping contact bin pair probability threshold to extact high quality multi-mapping contact.', default = 0.5)
    parser.add_argument('-f', '--filename', help = 'Multi-mapping contact read pair posterior probability assignment output file name. ".mHiC" will be added as suffix.', default = None)
    parser.add_argument('-o', '--outdir', help = 'Output directory to save results.', default = None)
    #parser.add_argument('-s', '--summary', help = '(Optional) Summary of multi-mapping bin pair results. Default is true.', default = True)
    #parser.add_argument('-sf', '--summaryFile', help = '(Optional) Summary file path. Default is NOne.', default = None)
    parser.add_argument('-v', '--verbose', help = '(Optional) Verbose. Default is true.', default = True)

    args = parser.parse_args()

    if args.prior is None or args.uni is None or args.multi is None or args.outdir is None:
        parser.print_help()
        sys.exit()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if not os.path.exists(args.prior):
        print("Prior file does not exist!")
        sys.exit()

    if not os.path.exists(args.uni):
        print("Uniquely contact bin pair file does not exist!")
        sys.exit()
    if not os.path.exists(args.multi):
        print("Multi-mapping contact bin pair file does not exist!")
        sys.exit()
    args.threshold = float(args.threshold)
    if args.threshold < 0 or args.threshold >=1:
        print("Probability threshold to extact high quality multi-mapping contact should be between 0 and 1 excluding 1 as we use strict sign > threshold to do the filtering.")
        sys.exit()
    if args.filename is  None:
        args.filename = args.multi.strip().split('/')[-1]
        
    return args

if __name__ == '__main__':

    args = get_args()
        
    if args.verbose:
        print("Prior file is = ", args.prior)
        print("Uniquely mapping contact bin pair file = ", args.uni)
        print("Multi-mapping contact bin pair file = ", args.multi)
        print("Multi-mapping contact bin pair keys = ", args.multikeys)
        print("Probability threshold = ", args.threshold)
        print("Output directory = ", args.outdir)
        print("Multi-mapping probability assignment output = ", args.filename + ".mHiC")
        #print("Summary = ", args.summary)
        #print("SummaryFile = ", args.summaryFile)
        print("Verbose = ", args.verbose)


    priorPath = args.prior #"/projects/s5/splineResults"
    multiFilePath=args.multi #"/projects/s4/rep.validPairs.MULTI.binPair.multi"
    multiKeysPath=args.multikeys #"/projects/s4/rep.validPairs.MULTI.binPair.multiKeys"
    uniFilePath= args.uni #"/projects/s4/validPairs.binPairCount.uni"
    outFilePath= args.outdir + "/" + args.filename + ".mHiC" #"/projects/s6/outFile.multi"
    threshold = args.threshold #0.5


    main(priorPath, multiFilePath, multiKeysPath, uniFilePath, outFilePath, threshold)

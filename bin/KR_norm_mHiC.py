#!/usr/bin/env python
# Author: Arya Kaul
# Modifies by Ye Zheng (yezheng@stat.wisc.edu)
# May 2018

import sys
import numpy as np
import os
import re
import scipy.sparse as sps
import warnings
warnings.filterwarnings('ignore')
import argparse
import time
import resource

def get_args():
    '''Get arguments'''
    parser = argparse.ArgumentParser(description = '------------Usage Start------------',
                                     epilog = '------------Usage End------------')

    parser.add_argument('-f', '--file', help="Path to the interaction file to be read.", default = None)
    parser.add_argument('-o', '--outdir', help='Path to directory for outputted files.', default = None)
    parser.add_argument('-r', '--resolution', help='Resolution of the interaction, i.e., 40000.', type=int, default = None)
    parser.add_argument('-l' ,'--chrLens', help='Input file for chromosome lengths.', default = None)
    parser.add_argument('-c', '--chrNum', help='Chromosome number to read in following format: chrNum. If whole genome reading is required, please use "whole".', default = "whole")
    parser.add_argument('-tr', '--toRemove', help='Percentage of sparse row/columns to remove from matrix. Default is to remove top 10% sparse elements.', default = 10, type = int)
    parser.add_argument('-v', '--verbose', help = '(Optional) Verbose. Default is true.', default = True)

    args = parser.parse_args()
    if args.file is None or args.outdir is None or args.resolution is None or args.chrLens is None:
        parser.print_help()
        sys.exit()

    if not os.path.exists(args.file):
        print("Interaction file does not exist!")
        sys.exit()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if not os.path.exists(args.chrLens):
        print("Chromosome length file does not exist!")
        sys.exit()

    if args.toRemove <= 0 or args.toRemove >=100:
        print("Percentage of sparse rows/columns to be removed should be within (0, 100).")
        sys.exit()

    return args


def knightRuizAlg(A, tol=1e-6, f1 = False):

    ##FUNCTION DESCRIPTION
    # knighRuizAlg is an implementation of the matrix balancing algorithm
    #  developed by Knight and Ruiz. The goal is to take a matrix A and
    #  find a vector x such that, diag(x)*A*diag(x) returns a doubly
    #  stochastic matrix

    ##PARAMETERS
    #A is a given numpy array
    #tol is error tolerance
    #f1 boolean indicating if the intermediate convergance statistics
    # should also be outputted
    n = A.shape[0]
    e = np.ones((n,1), dtype = np.float64)
    res = []


    Delta = 3
    delta = 0.1
    x0 = np.copy(e)
    g = 0.9

    etamax = eta = 0.1
    stop_tol = tol*0.5
    x = np.copy(x0)

    rt = tol**2.0
    v = x * (A.dot(x))
    rk = 1.0 - v
#    rho_km1 = np.dot(rk.T, rk)[0, 0]
    rho_km1 = ((rk.transpose()).dot(rk))[0,0]
    rho_km2 = rho_km1
    rout = rold = rho_km1
    
    MVP = 0 #we'll count matrix vector products
    i = 0 #outer iteration count

    if f1:
        print("it in. it res\n"),

    while rout > rt: #outer iteration
        i += 1

        if i > 30:
            break

        k = 0
        y = np.copy(e)
        innertol = max(eta ** 2.0 * rout, rt)
        
        while rho_km1 > innertol: #inner iteration by CG
            k += 1
            if k == 1:
                Z = rk / v
                p = np.copy(Z)
                #rho_km1 = np.dot(rk.T, Z)
                rho_km1 = (rk.transpose()).dot(Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            if k > 10:
                break

            #update search direction efficiently
            w = x * A.dot(x * p) + v * p
            # alpha = rho_km1 / np.dot(p.T, w)[0,0]
            alpha = rho_km1 / (((p.transpose()).dot(w))[0,0])
            ap = alpha * p
            #test distance to boundary of cone
            ynew = y + ap
            
            if np.amin(ynew) <= delta:
                
                if delta == 0:
                    break

                ind = np.where(ap < 0.0)[0]
                gamma = np.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            if np.amax(ynew) >= Delta:
                ind = np.where(ynew > Delta)[0]
                gamma = np.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            y = np.copy(ynew)
            rk -= alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            #rho_km1 = np.dot(rk.T, Z)[0,0]
            rho_km1 = ((rk.transpose()).dot(Z))[0,0]
        x *= y
        v = x * (A.dot(x))
        rk = 1.0 - v
        #rho_km1 = np.dot(rk.T, rk)[0,0]
        rho_km1 = ((rk.transpose()).dot(rk))[0,0]
        rout = rho_km1
        MVP += k + 1
        
        #update inner iteration stopping criterion
        rat = rout/rold
        rold = rout
        res_norm = rout ** 0.5
        eta_o = eta
        eta = g * rat
        if g * eta_o ** 2.0 > 0.1:
            eta = max(eta, g * eta_o ** 2.0)
        eta = max(min(eta, etamax), stop_tol / res_norm)
        if f1:
            print("%03i %06i %03.3f %e %e \n") % \
                (i, k, res_norm, rt, rout), 
            res.append(res_norm)
    if f1:
        print("Matrix - vector products = %06i\n") % \
            (MVP),
    
    #X = np.diag(x[:,0])   
    #x = X.dot(A.dot(X))
    return [x,(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000)]

def chromLength(chrLens, resolution):
    '''
    This function is to read in chromosome length file and store each bin as a dictionary.
    '''
    if chrLens is not None:
        import math
        allFragsDic = {}
        revFragsDic = []
        chrLenDic = {}
        fragsIndex = 0
        with open(chrLens, 'r') as infile: 
            for line in infile:
                ## for each chromosome
                chrom, chrL = line.split()
                chrLenDic[str(chrom)] = int(chrL)
                if chrom not in allFragsDic:
                    allFragsDic[chrom]={}
                for i in range(int(math.ceil(1.0 * int(chrL) / resolution))):
                    mid = int( resolution / 2 ) + i * resolution
                    allFragsDic[chrom][mid] = fragsIndex
                    revFragsDic.append("=".join([str(chrom), str(mid)]))
                    fragsIndex += 1
        infile.close()

    print("## Total number of bins:" + str(fragsIndex))
    return allFragsDic, revFragsDic, chrLenDic

def readInteraction(file, chrNum, resolution, chrLenDic, allFragsDic):
    '''
    This function is to read in interaction files.
    '''
    
    import math
    
    hicFile = open(file, 'r')
    
    print("Reading in the interaction file...")
    startTime = time.time()
    nBin = 0
    halfBin = resolution/2

    if chrNum == 'whole':
        for key in chrLenDic:
            nBin += int(math.ceil(1.0 * int(chrLenDic[key]) / resolution))

    else:
        nBin = int(math.ceil(1.0 * int(chrLenDic[chrNum]) / resolution))
    
    #construct a sparse matrix of max resolution
    hic_mtx = sps.lil_matrix((nBin,nBin), dtype = np.int64) 

    #load values into the sparse matrix 
    for line in hicFile:

        if chrNum != "whole" and line.startswith(chrNum):
            fileLine = line.rstrip().split()
            i = int((fileLine[1]-halfBin)/resolution)
            j = int((fileLine[3]-halfBin)/resolution)
            k = float(fileLine[4])


            if (fileLine[0] == fileLine[2]):
                try:
                    hic_mtx[i,j] = k
                    hic_mtx[j,i] = k

                except:
      #              print fileLine
                    continue

        else:
            fileLine=line.rstrip().split()
            chr1 = fileLine[0]
            chr2 = fileLine[2]

            mid1 = int(fileLine[1])
            mid2 = int(fileLine[3])
            k = float(fileLine[4])

            i = allFragsDic[chr1][mid1]
            j = allFragsDic[chr2][mid2]
            hic_mtx[i,j] = k
            hic_mtx[j,i] = k


    #convert to csr!
    hic_mtx = hic_mtx.tocsr()
    R = sps.csr_matrix.sum(hic_mtx)
    #done loading!
    endTime = time.time()
    print("Reading in the interaction file took %f seconds" % (endTime - startTime))

    hicFile.close()
    return R, hic_mtx

def removeRowCSR(mat, i):
    if not isinstance(mat, sps.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])
    return mat
    
def dropcols_coo(M, idx_to_drop):
    idx_to_drop = np.unique(idx_to_drop)
    C = M.tocoo()
    keep = ~np.in1d(C.col, idx_to_drop)
    C.data, C.row, C.col = C.data[keep], C.row[keep], C.col[keep]
    C.col -= idx_to_drop.searchsorted(C.col) # decrement column indices
    C._shape = (C.shape[0], C.shape[1] - len(idx_to_drop))
    return C.tocsr()

def removeSparseCSR(mtx, revFragsDic = None, perc=0):
    iteration = 0
    toRemove = []
    ctr = 0

    if perc == 0:
        print("Remove rows/columns whose diagonal is 0.")
        diagonal = mtx.diagonal()
        for values in diagonal:
            if values == 0:
                toRemove.append(ctr)
            ctr += 1
    
    else:
        rowSums = mtx.sum(axis=0)
        rowSums = list(np.array(rowSums).reshape(-1,)) 
        rowSums = list(enumerate(rowSums))
        for value in rowSums:
            if int(value[1]) == 0:
                toRemove.append(value[0])
                rowSums.remove(value) 

        rowSums.sort(key=lambda tup: tup[1])
        size = len(rowSums)
        prop = perc/100.0
        rem = int(prop * size)
        while ctr < rem:
            toRemove.append(rowSums[ctr][0])
            ctr += 1
    list(set(toRemove))
    toRemove.sort()
    mtx = dropcols_coo(mtx, toRemove)
    for num in toRemove:
        if iteration != 0:
            num -= iteration
        mtx = removeRowCSR(mtx, num)
        revFragsDic.pop(num)
        iteration +=1
    return mtx, revFragsDic, toRemove

def computeBiasVector(x):
    one = np.ones((x.shape[0], 1))
    x = one/x
    sums = np.sum(x)
    avg = (1.0 * sums)/x.shape[0]
    bias = np.divide(x, avg)
    return bias

def addZeroBiases(removeList, biasVector):
    for pos in removeList:
        biasVector = np.insert(biasVector, pos, -1, axis = 0)
        
    return biasVector

def writeInteraction(matrix, filename, outdir, revFragsDic, chrNum, resolution):
    outFile = os.path.join(outdir, filename + ".KRnorm")
    with open(outFile, 'w') as matrixFile:
        row, col = matrix.nonzero()
        values = matrix.data
        if chrNum == 'whole':
            for i in range(len(row)):
                chr1, mid1 = revFragsDic[row[i]].split("=")
                chr2, mid2 = revFragsDic[col[i]].split("=")
                matrixFile.write(("%s\t%s\t%s\t%s\t%.4f\n") % (chr1, mid1, chr2, mid2, values[i]))
        else:
            for i in range(len(row)):
                mid1 = row[i] * resolution + resolution/2
                mid2 = col[i] * resolution + resolution/2
                matrixFile.write(("%s\t%s\t%s\t%s\t%.4f\n") % (chrNum, mid1, chrNum, mid2, values[i]))

def writeBias(bias, filename, outdir, revFragsDicAll, chrNum, resolution):

    biasFile = os.path.join(outdir, filename + ".KRnorm.bias")
    
    with open(biasFile + "All", 'w') as biasWzero:
        with open(biasFile, 'w') as biasWOzero:
            if chrNum == 'whole':
                for i in range(len(bias)):
                    chr, mid = revFragsDicAll[i].split("=")
                    if bias[i, 0] != -1:
                        biasWOzero.write("\t".join([str(chr), str(mid), str(bias[i, 0])]) + "\n")
                    biasWzero.write("\t".join([str(chr), str(mid), str(bias[i, 0])]) + "\n")

            else:
                for i in range(len(bias)):
                    mid = i * resolution + resolution/2
                    if bias[i, 0] != -1:
                        biasWOzero.write("\t".join([str(chrNum), str(mid), str(bias[i, 0])]) + "\n")
                    biasWzero.write("\t".join([str(chrNum), str(mid), str(bias[i, 0])]) + "\n")

                    
        biasWOzero.close()
    biasWzero.close()

def normMatrix(rawMatrix, R, removePerc, revFragsDic):

    '''
    This function remove sparse elements, normalize the contact matrix, calculate the bias
    '''

    ## remove sparse elements
    print("## Raw matrix shape :" + str(rawMatrix.shape[0]) + "*" + str(rawMatrix.shape[1]))
    st=time.time()
    hic_mtx, revFragsDic, toRemove = removeSparseCSR(rawMatrix, revFragsDic, removePerc)
    e=time.time()
    print("Removing sparse regions took %s seconds" % (e-st))

    print("## Matrix shape after removing sparse bins :" + str(hic_mtx.shape[0]) + "*" + str(hic_mtx.shape[1]))
    print("## Remove list length:" + str(len(toRemove)))
    
    initialSize = rawMatrix.shape[0]
    newSize = hic_mtx.shape[0]
    # rawMatrix = mtxAndRemoved[0]
    # removed = mtxAndRemoved[1]

    ## normalize matrix
    print("Generating Normalized Matrix.")
    st=time.time()
    result = knightRuizAlg(hic_mtx)
    colVec = result[0]
    if np.isnan(np.sum(colVec)):
        print("Too few rows/columns removed... try again")
        return None
    x = sps.diags(colVec.flatten(), 0, format='csr')
    
    ## calculate bias
    bias = computeBiasVector(colVec)
    biasWZeros = addZeroBiases(toRemove, bias)
        

    del(colVec)
    hic_norm = x.dot(hic_mtx.dot(x))
    n = hic_norm.shape[0]
    scalar = R / (2.0 * n)
    norm_mtx = hic_norm * scalar
    e=time.time()
    print("Normalization took %s seconds" % (e-st))
    print("Normalized Matrix Generated")
    print("## Matrix shape after normalization :" + str(norm_mtx.shape[0]) + "*" + str(norm_mtx.shape[1]))
    
    return norm_mtx, biasWZeros, revFragsDic


if __name__ == "__main__":
    args = get_args()

    if args.verbose:
        print("Interaction file = ", args.file)
        print("Output directory = ", args.outdir)
        print("Resolution = ", args.resolution)
        print("Chromosome length file = ", args.chrLens)
        print("Chromosome to process = ", args.chrNum)
        print("Percentage of sparse elements to be removed = ", args.toRemove)


    baseName = os.path.basename(args.file)
    biasFile = args.outdir + "/" + baseName + ".bias"

    ## 1. Read in chromosome length file and create fragment dictionary
    allFragsDic, revFragsDic, chrLenDic = chromLength(args.chrLens, args.resolution)    
    revFragsDicAll = list(revFragsDic)
    ## 2. read interaction file
    R, hic_mtx = readInteraction(args.file, args.chrNum, args.resolution, chrLenDic, allFragsDic)    

    ## 3. normalize the contact matrix
    norm_mtx, bias, revFragsDic = normMatrix(hic_mtx, R, args.toRemove, revFragsDic)

    print("## Length of revFragsDic before sparse element removal:" + str(len(revFragsDicAll)))
    print("## Length of revFragsDic after sparse element removal:" + str(len(revFragsDic)))
    ## 4. output normalized interaction
    writeInteraction(norm_mtx, baseName, args.outdir, revFragsDic, args.chrNum, args.resolution)

    ## 5. output bias
    writeBias(bias, baseName, args.outdir, revFragsDicAll, args.chrNum, args.resolution)

    print("KR normalization done!")

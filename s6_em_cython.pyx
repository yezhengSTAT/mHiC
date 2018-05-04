import bisect
import time
import re
import numpy as np
cimport cython
from libc.stdlib cimport malloc, free
from cpython cimport array
import array
import sys
import os

cdef struct Zstruct:
    long *alignPair
    double *alignProb
    int alignLen

cdef struct Zrevstruct:
    long *ind
    int alignLen
    int uniCount
    
@cython.profile(False)
cdef void make_Zstruct(Zstruct *ZstructC, dict Z, long ZstructN):
    cdef int Zn
    cdef int i, j
    cdef long[:] Zkeys
    cdef double[:]  Zvalues
    
    for i in xrange(ZstructN):
        Zkeys = array.array("l", Z[i].keys())
        Zvalues = array.array("d", Z[i].values())
        Zn = len(Zkeys)
        
        ZstructC[i].alignPair = <long *>malloc(Zn*cython.sizeof(long))
        ZstructC[i].alignProb = <double *>malloc(Zn*cython.sizeof(double))
        ZstructC[i].alignLen = Zn
        
        if not ZstructC[i].alignPair:
            raise MemoryError()
        if not ZstructC[i].alignProb:
            raise MemoryError()
        
        for j in xrange(Zn):
            ZstructC[i].alignPair[j] = Zkeys[j]
            ZstructC[i].alignProb[j] = Zvalues[j]

@cython.profile(False)
cdef void make_Zrevstruct(Zrevstruct *ZrevstructC, dict Zrev, long ZrevstructN):
    cdef int Zrevn
    cdef int i, j
    cdef long[:] Zrevkeys
    
    for i in xrange(ZrevstructN):
        Zrevkeys = array.array("l", Zrev[i])
        Zrevn = len(Zrevkeys)
        
        ZrevstructC[i].ind = <long *>malloc(Zrevn*cython.sizeof(long))
        ZrevstructC[i].alignLen = Zrevn
        ZrevstructC[i].uniCount = 0
        if not ZrevstructC[i].ind:
            raise MemoryError()
        for j in xrange(Zrevn):  
            ZrevstructC[i].ind[j] = Zrevkeys[j]
            
@cython.profile(False)
def read_spline_prior(str splineFilePath):
    """
    read in genomic distance prior
    """
    cdef dict spline = {}
    cdef double interProb
    cdef str line
    cdef list splineLine

    splineFile = open(splineFilePath, "r")
    for line in splineFile:
        splineLine = line.rstrip().split()
        spline[int(splineLine[0])] = float(splineLine[1])
    interProb = min(spline.values())/2
    splineFile.close()
    return (spline, interProb)


    
@cython.profile(False)
cdef double get_spline_probability(dict spline, double interProb, str pairString):
    """
    init pi based on chr and pos distance
    """
    cdef long pos1, pos2, dist, splineLen, xPoint
    cdef str chr1, chr2
    cdef list pair
    
    pair = pairString.rstrip().split("_")
    ch1 = pair[0]
    pos1 = int(float(pair[1]))
    ch2 = pair[2]
    pos2 = int(float(pair[3]))
    dist = abs(pos2 - pos1)
    
    
    if ch1 != ch2:
        return interProb
    elif dist in spline.keys():
        return spline[dist]
    else:
        splineLen = len(spline)
        xPoint = min(bisect.bisect_left(list(spline.keys()), dist), splineLen-1)
        return spline[list(spline.keys())[xPoint]]


@cython.profile(False)
def init_multi(dict spline, double interProb, str multiFilePath, str multiKeysPath):
    """
    initialize multiFile
    """
    
    
    cdef long index = -1, pairID
    cdef dict Z = {}, Zrev = {}, qnameList = {}, pairList = {}
    cdef str qname = None, qnamePre = None, pair, line
    cdef list pairKeys, piInit, pairLine, pairIndex

    
    multiFile = open(multiFilePath, "r")
    multiKeysFile = open(multiKeysPath, "r")
    
    # Create Zrev[pair] dictionary

    pairKeys = [ eachPair.rstrip() for eachPair in multiKeysFile ]
    

    pairIndex = list(range(0, len(pairKeys)))
    pairList = dict(zip(pairKeys, pairIndex))

    piInit = [get_spline_probability(spline, interProb, pairKeysEach) for pairKeysEach in pairKeys]
    cdef double[:] prior = array.array("d", piInit)
    cdef double[:] pi = array.array("d", piInit)
    
    
    Zrev = dict.fromkeys(pairIndex)
    for pairindex in pairIndex:
        Zrev[pairindex] = []


    multiKeysFile.close()
    
    # start For
    for line in multiFile:
        # Assume the multi file are sorted by the qname
        pairLine = line.rstrip().split()
        qname = pairLine[0]
        pair = "_".join([pairLine[1], pairLine[2], pairLine[3], pairLine[4]])
        
        if qname != qnamePre:
            index += 1
            Z[index] = {}
            qnameList[index] = qname # create index - query name connection

        pairID = pairList[pair]
        Z[index][pairID] = 0
        Zrev[pairID].append(index)
        qnamePre = qname
    # end For
    
    multiFile.close()
    
    
    return(pi, Z, Zrev, pairKeys, pairList, prior, qnameList, index)

@cython.profile(False)
cdef long init_uni(str uniFilePath, Zrevstruct *ZrevstructC, dict pairList, int index, int totalReadNumIndex):
    """init uni-read"""
    cdef str line, pair
    cdef list pairLine

    uniFile = open(uniFilePath, "r")
    #Start For
    for line in uniFile:
        index += 1
            
            
        pairLine = line.rstrip().split()
        pair = "_".join([pairLine[0], pairLine[1], pairLine[2], pairLine[3]])
        try:
            pairID = pairList[pair]
            ZrevstructC[pairID].uniCount = int(pairLine[4]) 
            totalReadNumIndex += int(pairLine[4])
        except:
            pass
    #End For
    uniFile.close()
    return (totalReadNumIndex)


@cython.cdivision(True)
cdef  update_z(double[:] pi, Zstruct *ZstructC, Zrevstruct *ZrevstructC, long multiIndThre, double threshold, double[:] sumPi):
    
    cdef long selectChange = 0, i, tmp, ind
    cdef double prob, sumPitmp = 0.0
    cdef long* alignPair
    cdef double* alignProb
    cdef int alignPairLen

    

    for ind in range(0, multiIndThre+1):
        alignPair = ZstructC[ind].alignPair
        alignProb = ZstructC[ind].alignProb
        
        alignPairLen = ZstructC[ind].alignLen
        sumPitmp = 0.0
        for i in xrange(alignPairLen):
            sumPitmp += pi[alignPair[i]]


        for i in xrange(alignPairLen):
            tmp = alignPair[i]
            prob = pi[tmp]/sumPitmp
            if(alignProb[i] - threshold) * (prob - threshold) < 0:
                selectChange += 1
            alignProb[i] = prob
        sumPi[ind] = sumPitmp
    return (selectChange, sumPi)


@cython.cdivision(True)
cdef double[:] update_pi(double[:] pi, double[:] prior, Zrevstruct *ZrevstructC, long totalReadNum, double * diff, double[:] sumPi):
    cdef double pi_new, alignProbSum = 0.0
    cdef  long i
    cdef double * alignProb
    cdef long * indList
    for i in xrange(len(pi)):
        alignProbSum = 0.0
        indList = ZrevstructC[i].ind
        for j in xrange(ZrevstructC[i].alignLen):
            alignProbSum += 1.0/sumPi[indList[j]]
        
        pi_new = alignProbSum * pi[i] + ZrevstructC[i].uniCount + totalReadNum * prior[i]
        diff[0] += (pi[i] - pi_new)*(pi[i] - pi_new)
        pi[i] = pi_new
    return(pi)


cdef void em(double[:] pi, double[:] prior, Zstruct *ZstructC, Zrevstruct *ZrevstructC, long multiIndThre, double threshold, int maxIter, double diffThre, int selectChangeThre, long totalReadNum):
    
    cdef  int eachIter, selectChange
    cdef double *diff = <double *>malloc(1*cython.sizeof(double))
    cdef double[:] sumPi = array.array("d", np.zeros(multiIndThre+1))
    
    start_time = time.time()
    
    for eachIter in range(0, maxIter):
        

        print("EM iteration :" + str(eachIter) + " Time:" + str(time.time() - start_time))
        start_time = time.time()
        
        # Update Z
        print("Update Z!")
        start_time_tmp = time.time()
    
        (selectChange, sumPi) = update_z(pi, ZstructC, ZrevstructC, multiIndThre, threshold, sumPi)
        
        print("update Z time:" + str(time.time() - start_time_tmp))
        print("Iteration: " + str(eachIter) + " Multi-bin pair selection change number: " + str(selectChange) + "\n")
        
        # Update pi
        print("Update Pi!")
        start_time_tmp = time.time()
        diff[0] = 0.0
        pi = update_pi(pi, prior, ZrevstructC, totalReadNum, diff, sumPi)

        
        print("update pi iteration time:" + str(time.time() - start_time_tmp))
    
        print("Iteration: " + str(eachIter) + " Difference of pi: "  + str(np.sqrt(diff[0])) + "\n")
        if np.sqrt(diff[0]) < diffThre and selectChange < selectChangeThre:
            break

cdef void write_out(str outFilePath, Zstruct *ZstructC, list pairKeys, dict qnameList, long multiIndThre):
    
    cdef long * alignPair
    cdef double * alignProb

    outFile = open(outFilePath, "w+")
    for ind in range(0, multiIndThre + 1):
        alignPair = ZstructC[ind].alignPair
        alignProb = ZstructC[ind].alignProb
        for j in xrange(ZstructC[ind].alignLen):
            position = "\t".join(pairKeys[alignPair[j]].strip().split("_"))
            outString = "\t".join([qnameList[ind], position, str(alignProb[j])])
            outFile.write(outString + "\n")
    outFile.close()

def main(str priorPath, str multiFilePath, str multiKeysPath, str uniFilePath, str outFilePath, double threshold):

    cdef int maxIter, selectChangeThre
    cdef double diffThre

    maxIter = 500
    diffThre = 0.01
    selectChangeThre = 1

    
    (spline, interProb) = read_spline_prior(priorPath)

    (pi, Z, Zrev, pairKeys, pairList, prior, qnameList, index) = init_multi(spline, interProb, multiFilePath, multiKeysPath)
    multiIndThre = index
    cdef long ZstructN = index + 1
    cdef long ZrevstructN = len(pairKeys)
    
    
    cdef Zstruct *ZstructC = <Zstruct *> malloc(ZstructN*sizeof(Zstruct))
    if not ZstructC:
        raise MemoryError()

    cdef Zrevstruct *ZrevstructC = <Zrevstruct *> malloc(ZrevstructN*sizeof(Zrevstruct))
    if not ZrevstructC:
        raise MemoryError()

    make_Zstruct(ZstructC, Z, ZstructN)
    make_Zrevstruct(ZrevstructC, Zrev, ZrevstructN)

    totalReadNum = init_uni(uniFilePath, ZrevstructC, pairList, index, index)

    
    em(pi, prior, ZstructC, ZrevstructC, multiIndThre, threshold, maxIter, diffThre, selectChangeThre, totalReadNum)

    write_out(outFilePath, ZstructC, pairKeys, qnameList, multiIndThre)
    

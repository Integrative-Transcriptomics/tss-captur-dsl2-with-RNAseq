"""Converter from .wig to .bedGraph file format

This script merges all wiggle files and then converts them into a bedGraph file. 
It calculates the full length of the genome in the wiggle file by using the appropriate fasta file.
"""

import sys
import MergeWigs

def WigToBedGraph(Wigglepath, chromSizesFilePath, outputPath):

    correctWalker = MergeWigs.merge_wigs(Wigglepath)

    BedGraphFile = open(f'{outputPath}.bedGraph', 'w')

    sizes = open(chromSizesFilePath, 'r')

    lastVal = 0
    lastReg = ""
    lastPos = 0
    first = True
    start = 0
    totalsum = 0

    sizeDict = {}

    for line in sizes:
        parts = line.split()
        sizeDict[parts[0]] = int(parts[1])

    trueCurrRegLength = 0
    for reg, pos, val in correctWalker:
        #handle minus of reverse wigs
        val = abs(val)
        if(first):
            lastPos = pos
            lastReg = reg
            lastVal = val
            start = pos
            first = False
            if(pos != 1):
                BedGraphFile.write(f'{lastReg}\t1\t{pos}\t0\n')
            trueCurrRegLength = sizeDict[reg]

        #If we are entering a new region
        if(reg != lastReg):
            #value might not have changed but still need to fill the current tracked intervall
            BedGraphFile.write(f'{lastReg}\t{start}\t{lastPos}\t{lastVal}\n')
            
            #If true chromosome size > lastPos, there are 0 at end of region which were cut off from .wig
            if(trueCurrRegLength > lastPos):
                BedGraphFile.write(f'{lastReg}\t{lastPos}\t{trueCurrRegLength}\t0\n')
            
            #if next chromosome does not start with position 1, fill missing poses with 0
            if(pos != 1):
                BedGraphFile.write(f'{reg}\t1\t{pos}\t0\n')

            end = pos
            lastVal = val
            lastReg = reg
            start = pos
            trueCurrRegLength = sizeDict[reg]
            
        #if we still are in the same region
        elif(val != lastVal):        
            BedGraphFile.write(f'{lastReg}\t{start}\t{pos}\t{lastVal}\n')

            end = pos
            lastVal = val
            lastReg = reg
            start = pos

        lastPos = pos
        
        totalsum += val

    #close final intervall of final chromosome
    BedGraphFile.write(f'{reg}\t{start}\t{pos}\t{lastVal}\n')

    #fill potentially end-cut off 0s from chromosome
    if(trueCurrRegLength > lastPos):
        BedGraphFile.write(f'{lastReg}\t{lastPos}\t{trueCurrRegLength}\t0\n')

    print(BedGraphFile)
    return BedGraphFile


if __name__ == '__main__':
    WigToBedGraph(sys.argv[1], sys.argv[2], sys.argv[3])
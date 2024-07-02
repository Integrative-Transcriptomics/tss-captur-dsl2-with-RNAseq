import glob
import itertools
import gff3_parser
import pyBigWig as bigwig
import numpy as np

def CalcBackgroundNoise(annotationPath, bigWigPath):

    bw = bigwig.open(bigWigPath)


    files = glob.glob(f"{annotationPath}/*.gff")
    print(f"annnot path bg: {annotationPath}")
    print(f"anoot globs bg: {files}")
    annot = gff3_parser.parse_gff3(files[0], verbose = False, parse_attributes = False)


    chromosome = annot.loc[1, "Seqid"]


    inverse = GetInverseOfGFF(annotationPath, bw.chroms(chromosome))

    #print(np.quantile(valueCopy, 1))
    #print(nonTranscr)

    quantiles = []
    for intervall in inverse:
        quantiles.append(np.quantile(bw.values(chromosome, intervall[0], intervall[1]), 0.5))
        
    return(np.quantile(quantiles, 0.5))



def GetInverseOfGFF(annotationPath, chromSize):


    files = glob.glob(f"{annotationPath}/*.gff")
    print(f"annnot path inverse: {annotationPath}")
    print(f"anoot globs inverse: {files}")
    annot = gff3_parser.parse_gff3(files[0], verbose = False, parse_attributes = False)  

    acceptedtypes = ["CDS", "gene"]

    #valueCopy = bw.values(chromosome, 1, bw.chroms(chromosome))
    transcr = []
    #print(bw.chroms(chromosome))

    #print(len(valueCopy))
    #might need to filter for certain regions here as well
    for start, end, feature, strand in zip(annot.loc[:, "Start"].astype(int).tolist(), annot.loc[:, "End"].astype(int).tolist(), annot.loc[:, "Type"],annot.loc[:, "Strand"]):
        if(strand != "+" or (not feature in acceptedtypes)):
            continue
        #valueCopy[start:end] = itertools.repeat("out", (end - start))
        transcr.append((start,end))
    #valueCopy = [i for i in valueCopy if i!= "out"]

    start = 0

    # Sort intervals by start position
    transcr.sort(key=lambda x: x[0])

    fused_transcr = []
    current_start, current_end = transcr[0]

    for start, end in transcr[1:]:
        # If there's no gap, extend the current interval
        if start <= current_end:
            current_end = max(current_end, end)
        # If there's a gap, add the current interval and start a new one
        else:
            fused_transcr.append((current_start, current_end))
            current_start, current_end = start, end

    # Add the last interval
    fused_transcr.append((current_start, current_end))

    inverse = []

    ##if first transcribed intervall starts after 1, pre intervall region is not trranscribed
    if(1 < fused_transcr[0][0]):
        inverse.append((1, fused_transcr[0][0]))

    for en in enumerate(fused_transcr):
        if(en[0] < (len(fused_transcr) -1)):
            #add intervall from end of current intervall till start of next one
            inverse.append((en[1][1], fused_transcr[en[0] + 1][0]))

    #if end of last transcribed intervall is before end of chromosome, add last untranscribed intervall till end
    if(chromSize > fused_transcr[-1][1]):
        inverse.append((fused_transcr[-1][1], chromSize))
    
    return inverse
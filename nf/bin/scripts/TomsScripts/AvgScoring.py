from collections import OrderedDict
from dataclasses import dataclass
import CalcbackgroundNoise
import pyBigWig as bigwig
import gff3_parser
import numpy as np
import pandas as pd

@dataclass
class ScoredTerm:
    initialScore = 0
    AvgScore = 0
    DerivScore = 0
    start = -1
    end = -1

def FindFirstUpstreamTSS(position, tssList):
    lastPos = -1
    for tss in tssList:
        if(tss >= position):
            return lastPos
        lastPos = tss
    return lastPos


def AvgScoreTerminators(termGff, bigwigpath, MasterTable, annotgff):
    #Params
    WindowSize = 25
    WindowOffsetFromEnd = 10

    rhotermdata = gff3_parser.parse_gff3(termGff, verbose = False, parse_attributes = False)
    bw = bigwig.open(bigwigpath)
    MasterTable = pd.read_csv(MasterTable, sep="\t")

    tssStarts = list(OrderedDict.fromkeys(MasterTable.loc[:, "SuperPos"].tolist()))

    chromosome = rhotermdata.loc[1, "Seqid"] + ".1"
    print(chromosome)
    rhoTermintervalls = zip(rhotermdata.loc[:, "Start"].astype(int), rhotermdata.loc[:, "End"].astype(int), rhotermdata.loc[:, "Score"])

    noiseLvL = CalcbackgroundNoise.CalcBackgroundNoise(annotgff, bigwigpath)

    TSSTermPairing = {}

    for start, end, score in rhoTermintervalls:
        myTss =  FindFirstUpstreamTSS(start, tssStarts)
        estGeneExprMed = np.quantile(bw.values(chromosome, myTss, end), 0.5)

        #if gene is not expressed, we cannot judge its terminators
        if(estGeneExprMed < noiseLvL):
            continue

        windStart = end + WindowOffsetFromEnd
        windEnd = min(windStart + WindowSize, bw.chroms(chromosome))

        #Curious as to which quantile is best here, im just guessing right now
        postTermExprQ = np.quantile(bw.values(chromosome, windStart , windEnd), 1)

        scored = ScoredTerm()
        #The way we score here is important, no idea whats best
        scored.AvgScore = 15*noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)
        scored.start = start
        scored.end = end
        scored.initialScore = score

        if(myTss in TSSTermPairing):
            TSSTermPairing[myTss].append(scored)
        else:
            TSSTermPairing[myTss] = [scored]
        
        
    return TSSTermPairing
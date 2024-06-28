from collections import OrderedDict
from dataclasses import dataclass
import re
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
    seqid = 'Region here'
    strand = 'NA'
    type = 'terminator type here'
    start = -1
    end = -1

def FindFirstUpstreamTSS(position, tssList):
    lastPos = -1
    for tss in tssList:
        if(tss >= position):
            return lastPos
        lastPos = tss
    return lastPos


def ScoreArea(WindowOffsetFromEnd, WindowSize, scoredTerm, bwFile, noiseLvL):
    windStart = scoredTerm.end + WindowOffsetFromEnd
    windEnd = min(windStart + WindowSize, bwFile.chroms(scoredTerm.seqid))
    postTermExprQ = np.quantile(bwFile.values(scoredTerm.seqid, windStart, windEnd), 1)

    return 15*noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)


def AvgScoreTerminators(gffRhoterm, gffNocornac, bigwigpath, MasterTable, annotgff):
    #Params
    WindowSize = 25
    WindowOffsetFromEnd = 10

    rhotermdata = gff3_parser.parse_gff3(gffRhoterm, verbose = False, parse_attributes = False)
    nocornacdata = gff3_parser.parse_gff3(gffNocornac, verbose = False, parse_attributes = False)

    all_term_data = pd.concat([rhotermdata, nocornacdata], ignore_index = True)

    bw = bigwig.open(bigwigpath)
    MasterTable = pd.read_csv(MasterTable, sep="\t")

    tssStarts = list(OrderedDict.fromkeys(MasterTable.loc[:, "SuperPos"].tolist()))

    chromosome = re.match('^[^\.]+', rhotermdata.loc[1, "Seqid"]).group(0)

    for chrom in bw.chroms():
        prefix = re.match('^[^\.]+', chrom).group(0)
        if(prefix == chromosome):
            chromosome = chrom
            break
    print(chromosome)

    rhoTermintervalls = zip(rhotermdata.loc[:, "Start"].astype(int), rhotermdata.loc[:, "End"].astype(int), rhotermdata.loc[:, "Score"], rhotermdata.loc[:, "Type"],rhotermdata.loc[: "Strand"])

    noiseLvL = CalcbackgroundNoise.CalcBackgroundNoise(annotgff, bigwigpath)

    TSSTermPairing = {}

    for start, end, score, type, strand in rhoTermintervalls:
        if(type != 'terminator' and type != 'RhoTerminator'):
            print(type)
            continue

        myTss =  FindFirstUpstreamTSS(start, tssStarts)
        print(f"Intervalls: {myTss} {end} {tssStarts}")
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
        scored.start = start
        scored.end = end
        scored.seqid = chromosome
        scored.strand = strand
        scored.initialScore = score
        scored.type = type

        #Use specific scoring windows depending on Rho dependent or intrinsic terminator

        #intrinsic terminator
        if(type == "terminator"):
            scored.AvgScore = ScoreArea(WindowOffsetFromEnd=WindowOffsetFromEnd, WindowSize=WindowSize, scoredTerm=scored, bwFile=bw, noiseLvL=noiseLvL)
        #Rho dep. terminator
        else:
             scored.AvgScore = ScoreArea(WindowOffsetFromEnd=WindowOffsetFromEnd, WindowSize=WindowSize, scoredTerm=scored, bwFile=bw, noiseLvL=noiseLvL)

        #scored.AvgScore = 15*noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)

        if(myTss in TSSTermPairing):
            TSSTermPairing[myTss].append(scored)
        else:
            TSSTermPairing[myTss] = [scored]
        
        
    return TSSTermPairing
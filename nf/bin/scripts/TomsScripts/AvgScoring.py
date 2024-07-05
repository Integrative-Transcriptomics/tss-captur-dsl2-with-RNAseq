from collections import OrderedDict
from dataclasses import dataclass
import glob
import re
import CalcbackgroundNoise
import pyBigWig as bigwig
import gff3_parser
import numpy as np
import pandas as pd

@dataclass
class ScoredTerm:
    initialScore = -1
    avgScore = -1
    derivScore = -1
    seqid = 'Region here'
    strand = 'NA'
    type = 'terminator type here'
    start = -1
    end = -1

def FindFirstUpstreamTSS(position, tssList, strand):
    lastPos = -1
    if strand == '+':
        for tss in tssList:
            if(tss >= position):
                return lastPos
            lastPos = tss
        return lastPos
    else:
        for tss in tssList:
            if(tss > position):
                return tss
    return -1


def ScoreArea(WindowOffsetFromEnd, WindowSize, startOfArea, scoredTerm, bwFile, noiseLvL):
    if(scoredTerm.strand == '-'):
        WindowOffsetFromEnd *= -1
        WindowSize *=-1
        windStart = startOfArea + WindowOffsetFromEnd
        windEnd = max(windStart + WindowSize, 1)
    else:
        windStart = startOfArea + WindowOffsetFromEnd
        windEnd = min(windStart + WindowSize, bwFile.chroms(scoredTerm.seqid))
    
    if(windStart > windEnd):
        b = windEnd
        windEnd = windStart
        windStart = b

    postTermExprQ = np.quantile(bwFile.values(scoredTerm.seqid, windStart, windEnd), 1)

    return noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)


def AvgScoreTerminators(gffRhoterm, gffNocornac, bigwigpath, MasterTable, annotgff):
    #Params
    window_size = 25
    scoring_window_offset_rho = 10
    scoring_window_offset_intrinsic = 10


    rhotermdata = gff3_parser.parse_gff3(gffRhoterm, verbose = False, parse_attributes = False)
    nocornacdata = gff3_parser.parse_gff3(gffNocornac, verbose = False, parse_attributes = False)

    all_term_data = pd.concat([rhotermdata, nocornacdata], ignore_index = True)


    print(f"shape of concated: {all_term_data.shape}" )
    print(f"shape of concated: {rhotermdata.shape}" )
    print(f"shape of concated: {nocornacdata.shape}" )
    #throw = TSSTermPairing[99999999999]

    bw = bigwig.open(bigwigpath)
    MasterTable = pd.read_csv(MasterTable, sep="\t")

    positive_tss = sorted(MasterTable[MasterTable['SuperStrand'] == '+']['SuperPos'].tolist())
    negative_TSS = sorted(MasterTable[MasterTable["SuperStrand"] == '-']['SuperPos'].tolist())

    chromosome = re.match('^[^\.]+', rhotermdata.loc[1, "Seqid"]).group(0)

    for chrom in bw.chroms():
        prefix = re.match('^[^\.]+', chrom).group(0)
        if(prefix == chromosome):
            chromosome = chrom
            break
    print(chromosome)
    
    all_terms_intervalls = zip(all_term_data.loc[:, "Start"].astype(int), all_term_data.loc[:, "End"].astype(int), all_term_data.loc[:, "Score"], all_term_data.loc[:, "Type"], all_term_data.loc[: ,"Strand"])

    print(f"FROM AVGSCORING: {annotgff}")

    noiseLvL = CalcbackgroundNoise.CalcBackgroundNoise(annotgff, bigwigpath)

    TSSTermPairing = {}
    print(f"noiseLVL avg: {noiseLvL}")
    x = 0
    for start, end, score, type, strand in all_terms_intervalls:
        x += 1
        if(type != 'terminator' and type != 'RhoTerminator'):
            continue
        
        if(strand == '+'):
            myTss =  FindFirstUpstreamTSS(start, positive_tss, '+')
            estGeneExprMed = np.quantile(bw.values(chromosome, myTss, end), 0.5)
        else:
            myTss =  FindFirstUpstreamTSS(end, negative_TSS, '-')
            #print(f"Start: {start} End: {end} tss: {myTss}")
            estGeneExprMed = np.quantile(bw.values(chromosome, start, myTss), 0.5)
        #print(f"Intervalls: {myTss} {end} {tssStarts}")

        scored = ScoredTerm()
        #if gene is not expressed, we cannot judge its terminators
        if(estGeneExprMed <= noiseLvL):
            scored = ScoredTerm()
            
            scored.start = start
            scored.end = end
            scored.seqid = chromosome
            scored.strand = strand
            scored.initialScore = score
            scored.type = type
            scored.avgScore = "NA"
            if(myTss in TSSTermPairing):
                TSSTermPairing[myTss].append(scored)
            else:
                TSSTermPairing[myTss] = [scored]

            continue
        else:
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
            scored.avgScore = ScoreArea(WindowOffsetFromEnd=scoring_window_offset_intrinsic, WindowSize=window_size, startOfArea= scored.end, scoredTerm=scored, bwFile=bw, noiseLvL=noiseLvL)
        #Rho dep. terminator
        else:
            scored.avgScore = ScoreArea(WindowOffsetFromEnd=scoring_window_offset_rho, WindowSize=window_size, startOfArea= scored.end ,scoredTerm=scored, bwFile=bw, noiseLvL=noiseLvL)

        #scored.AvgScore = 15*noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)

        if(myTss in TSSTermPairing):
            TSSTermPairing[myTss].append(scored)
        else:
            TSSTermPairing[myTss] = [scored]
    
    
    print(f"Length of all term intervalls: {x}")

    return TSSTermPairing
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


def ScoreArea(WindowOffsetFromEnd, WindowSize, startOfArea, scoredTerm, bwFile, noiseLvL, iqr):
    if(scoredTerm.strand == '-'):
        WindowSize *=-1
        WindowOffsetFromEnd *= -1
        windStart = max(startOfArea + WindowOffsetFromEnd,1)
        windEnd = max(windStart + WindowSize, 1)
    else:
        windStart = min(startOfArea + WindowOffsetFromEnd, bwFile.chroms(scoredTerm.seqid))
        windEnd = min(windStart + WindowSize, bwFile.chroms(scoredTerm.seqid))
    
    if(windStart > windEnd):
        b = windEnd
        windEnd = windStart
        windStart = b

    postTermExprQ = np.quantile(bwFile.values(scoredTerm.seqid, windStart, windEnd), 1)

    upper_bound = noiseLvL + iqr * 2

    if(postTermExprQ <= noiseLvL):
        return 1
    elif(postTermExprQ >= upper_bound):
        return 0
    else:
         return (1 - ((postTermExprQ - noiseLvL) / (upper_bound - noiseLvL)))


    #return noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)


def AvgScoreTerminators(gffRhoterm, gffNocornac, forward_bigwig_path, reverse_bigwig_path, master_table_path, annotgff):
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

    fwbw = bigwig.open(forward_bigwig_path)
    rvbw = bigwig.open(reverse_bigwig_path)

    MasterTable = pd.read_csv(master_table_path, sep="\t")

    positive_tss = sorted(MasterTable[MasterTable['SuperStrand'] == '+']['SuperPos'].tolist())
    negative_TSS = sorted(MasterTable[MasterTable["SuperStrand"] == '-']['SuperPos'].tolist())

    forward_chromosome = re.match('^[^\.]+', rhotermdata.loc[1, "Seqid"]).group(0)

    for chrom in fwbw.chroms():
        prefix = re.match('^[^\.]+', chrom).group(0)
        if(prefix == forward_chromosome):
            forward_chromosome = chrom
            break

    reverse_chromosome = re.match('^[^\.]+', rhotermdata.loc[1, "Seqid"]).group(0)

    for chrom in rvbw.chroms():
        prefix = re.match('^[^\.]+', chrom).group(0)
        if(prefix == reverse_chromosome):
            reverse_chromosome = chrom
            break

    print(forward_chromosome, reverse_chromosome)
    
    all_terms_intervalls = zip(all_term_data.loc[:, "Start"].astype(int), all_term_data.loc[:, "End"].astype(int), all_term_data.loc[:, "Score"], all_term_data.loc[:, "Type"], all_term_data.loc[: ,"Strand"])

    print(f"FROM AVGSCORING: {annotgff}")

    forward_noise_lvl, forward_IQR = CalcbackgroundNoise.InverseOfMasterTableNoise(annotgff, forward_bigwig_path, master_table_path)
    reverse_noise_lvl, reverse_IQR = CalcbackgroundNoise.InverseOfMasterTableNoise(annotgff, reverse_bigwig_path, master_table_path)


    TSSTermPairing = {}
    print(f"noiseLVL avg: {forward_noise_lvl}")
    for start, end, score, type, strand in all_terms_intervalls:
        if(type != 'terminator' and type != 'RhoTerminator'):
            continue

        scored = ScoredTerm()
        scored.start = start
        scored.end = end        
        scored.strand = strand
        scored.initialScore = score
        scored.type = type

        if(strand == '+'):
            bw = fwbw
            chromosome = forward_chromosome
            noiseLvL = forward_noise_lvl
            iqr = forward_IQR
            scoreStartArea = scored.end
            myTss =  FindFirstUpstreamTSS(start, positive_tss, '+')
            estGeneExprMed = np.quantile(bw.values(chromosome, myTss, end), 0.5)
        else:
            bw = rvbw
            reverse_chromosome = reverse_chromosome
            noiseLvL = reverse_noise_lvl
            iqr = reverse_IQR
            scoreStartArea = scored.start
            myTss =  FindFirstUpstreamTSS(end, negative_TSS, '-')
            #print(f"Start: {start} End: {end} tss: {myTss}")
            estGeneExprMed = np.quantile(bw.values(chromosome, start, myTss), 0.5)
        #print(f"Intervalls: {myTss} {end} {tssStarts}")

        scored.seqid = chromosome

        #if gene is not expressed, we cannot judge its terminators
        if(estGeneExprMed <= noiseLvL):
            scored.avgScore = "NA"

            if(myTss in TSSTermPairing):
                TSSTermPairing[myTss].append(scored)
            else:
                TSSTermPairing[myTss] = [scored]
            continue

        #Use specific scoring windows depending on Rho dependent or intrinsic terminator

        #intrinsic terminator
        if(type == "terminator"):
            scored.avgScore = ScoreArea(WindowOffsetFromEnd=scoring_window_offset_intrinsic, WindowSize=window_size, startOfArea= scoreStartArea, scoredTerm=scored, bwFile=bw, noiseLvL=noiseLvL, iqr=iqr)
        #Rho dep. terminator
        else:
            scored.avgScore = ScoreArea(WindowOffsetFromEnd=scoring_window_offset_rho, WindowSize=window_size, startOfArea= scoreStartArea,scoredTerm=scored, bwFile=bw, noiseLvL=noiseLvL, iqr=iqr)

        #scored.AvgScore = 15*noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)
        if(myTss in TSSTermPairing):
            TSSTermPairing[myTss].append(scored)
        else:
            TSSTermPairing[myTss] = [scored]

    onlynas= 0
    for key in TSSTermPairing:
        if(all(term.avgScore == 'NA' for term in TSSTermPairing[key])):
            print(key)
            onlynas += 1

    print(f"Total TSS: {len(TSSTermPairing.keys())} tss with only NA: {onlynas}" )

    return TSSTermPairing
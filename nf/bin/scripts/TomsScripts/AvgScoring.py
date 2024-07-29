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
    dropScore = -1
    seqid = 'Region here'
    strand = 'NA'
    type = 'terminator type here'
    start = -1
    end = -1
    avgScoreStart = -1
    avgScoreEnd = -1
    dropScorePreTermStart = -1
    dropScorePreTermEnd = -1
    dropScorePostTermStart = -1
    dropScorePostTermEnd = -1
    minExprDrop = -1

def FindFirstUpstreamTSS(position, minimum_gene_length, tssList, strand):
    lastPos = -1
    if strand == '+':
        for tss in tssList:
            if(tss >= position):
                return lastPos
            
            if(tss <= position - minimum_gene_length):
                lastPos = tss
        return lastPos
    else:
        for tss in tssList:
            if(tss > position + minimum_gene_length):
                return tss
    return -1


def DropScoreArea(window_gene_ratio, post_term_window_offset, post_term_window_size, scoredTerm, myTss, bwFile, noiseLvl):

    chrom = scoredTerm.seqid
    minimum_expr_drop_ratio = 0.25

    if(scoredTerm.strand == '+'):
        score_pre_term_start = int(scoredTerm.start - ((scoredTerm.start - myTss) * window_gene_ratio))
        score_pre_term_end = scoredTerm.start

        score_post_term_start = min(scoredTerm.end + post_term_window_offset, bwFile.chroms(chrom) -1)
        score_post_term_end = min(score_post_term_start + post_term_window_size, bwFile.chroms(chrom))

    else:
        score_pre_term_start = scoredTerm.end
        score_pre_term_end = int(scoredTerm.end + ((myTss - scoredTerm.end) * window_gene_ratio))

        score_post_term_end = max(scoredTerm.start - post_term_window_offset, 2)
        score_post_term_start = max(score_post_term_end - post_term_window_size, 1)

    print("bounds: ", score_pre_term_start, score_pre_term_end, bwFile.chroms(chrom), scoredTerm.strand, scoredTerm.start, myTss, scoredTerm.end)
    pre_term_expr = np.quantile(bwFile.values(chrom, score_pre_term_start, score_pre_term_end), 0.5)

    post_term_expr = np.quantile(bwFile.values(chrom, score_post_term_start, score_post_term_end), 0.5)

    expression_ratio = post_term_expr / pre_term_expr
    
    score = 1 if expression_ratio < minimum_expr_drop_ratio else 0

    return score, score_pre_term_start, score_pre_term_end, score_post_term_start, score_post_term_end, (pre_term_expr * minimum_expr_drop_ratio)



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

    score = 0

    if(postTermExprQ <= noiseLvL):
        score = 1
    elif(postTermExprQ >= upper_bound):
        score = 0
    else:
         score = (1 - ((postTermExprQ - noiseLvL) / (upper_bound - noiseLvL)))

    return score, windStart, windEnd

    #return noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)


def AvgScoreTerminators(gffRhoterm, gffNocornac, forward_bigwig_path, reverse_bigwig_path, master_table_path, annotgff):

    #Params avg scoring
    AVG_WINDOW_SIZE = 25
    AVGSCORING_WINDOW_OFFSET_RHO = 10
    AVGSCORING_WINDOW_OFFSET_INTRINSIC = 10

    #Params drop score
    DROPSCORING_PRE_TERM_WINDOW_RATIO = 0.5
    DROPSCORING_POST_TERM_WINDOW_SIZE = 50

    DROPSCORING_POST_TERM_OFFSET_RHO = 10
    DROPSCORING_POST_TERM_OFFSET_INTRINSIC = 10

    MINIMUM_GENE_LENGTH = 30
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
    
    all_terms_intervalls = zip(all_term_data.loc[:, "Start"].astype(int), all_term_data.loc[:, "End"].astype(int), all_term_data.loc[:, "Score"],
                                all_term_data.loc[:, "Type"], all_term_data.loc[: ,"Strand"])

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
            avgScoreStartArea = scored.end
            myTss = FindFirstUpstreamTSS(start, MINIMUM_GENE_LENGTH, positive_tss, '+')
            #No TSS found, no gene here, we assume
            if(myTss <= -1):
                continue
            estGeneExprMed = np.quantile(bw.values(chromosome, myTss, end), 0.5)
        else:
            bw = rvbw
            chromosome = reverse_chromosome
            noiseLvL = reverse_noise_lvl
            iqr = reverse_IQR
            avgScoreStartArea = scored.start
            myTss = FindFirstUpstreamTSS(end, MINIMUM_GENE_LENGTH, negative_TSS, '-')
            #No TSS found, no gene here, we assume
            if(myTss <= -1):
                continue
            #print(f"Start: {start} End: {end} tss: {myTss}")
            #No TSS found, no gene here, we assume
            if(myTss <= -1):
                continue
            estGeneExprMed = np.quantile(bw.values(chromosome, start, myTss), 0.5)

        

        #print(f"Intervalls: {myTss} {end} {tssStarts}")

        scored.seqid = chromosome

        #if gene is not expressed, we cannot judge its terminators
        if(estGeneExprMed <= noiseLvL):
            scored.avgScore = "NA"
            scored.dropScore = "NA"

            if(myTss in TSSTermPairing):
                TSSTermPairing[myTss].append(scored)
            else:
                TSSTermPairing[myTss] = [scored]
            continue

        #Use specific scoring windows depending on Rho dependent or intrinsic terminator

        #intrinsic terminator
        if(type == "terminator"):
            scored.avgScore, scored.avgScoreStart, scored.avgScoreEnd = ScoreArea(WindowOffsetFromEnd=AVGSCORING_WINDOW_OFFSET_INTRINSIC, WindowSize=AVG_WINDOW_SIZE, startOfArea= avgScoreStartArea, scoredTerm=scored, bwFile=bw, noiseLvL=noiseLvL, iqr=iqr)

            scored.dropScore, scored.dropScorePreTermStart, scored.dropScorePreTermEnd, scored.dropScorePostTermStart, scored.dropScorePostTermEnd, scored.minExprDrop = DropScoreArea(
                window_gene_ratio=DROPSCORING_PRE_TERM_WINDOW_RATIO, post_term_window_size= DROPSCORING_POST_TERM_WINDOW_SIZE, post_term_window_offset=DROPSCORING_POST_TERM_OFFSET_INTRINSIC,
                scoredTerm=scored, myTss=myTss, bwFile=bw, noiseLvl=noiseLvL
            )
        #Rho dep. terminator
        else:
            scored.avgScore, scored.avgScoreStart, scored.avgScoreEnd = ScoreArea(WindowOffsetFromEnd=AVGSCORING_WINDOW_OFFSET_RHO, WindowSize=AVG_WINDOW_SIZE, startOfArea= avgScoreStartArea,scoredTerm=scored, bwFile=bw, noiseLvL=noiseLvL, iqr=iqr)

            scored.dropScore, scored.dropScorePreTermStart, scored.dropScorePreTermEnd, scored.dropScorePostTermStart, scored.dropScorePostTermEnd, scored.minExprDrop = DropScoreArea(
                window_gene_ratio=DROPSCORING_PRE_TERM_WINDOW_RATIO, post_term_window_size= DROPSCORING_POST_TERM_WINDOW_SIZE, post_term_window_offset=DROPSCORING_POST_TERM_OFFSET_RHO,
                scoredTerm=scored, myTss=myTss, bwFile=bw, noiseLvl=noiseLvL
            )


        #scored.AvgScore = 15*noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)
        if(myTss in TSSTermPairing):
            TSSTermPairing[myTss].append(scored)
        else:
            TSSTermPairing[myTss] = [scored]

    return TSSTermPairing
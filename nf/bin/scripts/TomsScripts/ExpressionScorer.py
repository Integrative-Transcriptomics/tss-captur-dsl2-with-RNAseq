import csv
import re
import sys
import CalcbackgroundNoise
import pyBigWig as bigwig
import gff3_parser
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly

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

def score_area(start, end, bwFile, chromosome, noiseLvL, iqr):

    postTermExprQ = np.quantile(bwFile.values(chromosome, start, end), 1)

    upper_bound = noiseLvL + iqr * 2

    score = 0

    if(postTermExprQ <= noiseLvL):
        score = 1
    elif(postTermExprQ >= upper_bound):
        score = 0
    else:
         score = (1 - ((postTermExprQ - noiseLvL) / (upper_bound - noiseLvL)))

    return score

def avg_score_area(term_start, term_end, window_size, window_offset, chromosome , strand, bwFile, noiseLvL, iqr):

    if(strand == '-'):
        window_size *=-1
        window_offset *= -1
        windStart = max(term_start + window_offset,1)
        windEnd = max(windStart + window_size, 1)
    else:
        windStart = min(term_end + window_offset, bwFile.chroms(chromosome))
        windEnd = min(windStart + window_size, bwFile.chroms(chromosome))
    
    if(windStart > windEnd):
        b = windEnd
        windEnd = windStart
        windStart = b

    return score_area(windStart, windEnd, bwFile, chromosome, noiseLvL, iqr), windStart, windEnd
    #return noiseLvL / max(noiseLvL, postTermExprQ, 0.00001)


def deriv_score_area(term_start, term_end, fit_window_size, scoring_window_size, scoring_window_offset, chromosome , strand, bwFile, noiseLvL, iqr):    
    
    if(strand == '+'):
        fit_window_start = max(term_start - fit_window_size, 1)
        fit_window_end = min(term_end + fit_window_size, bwFile.chroms(chromosome) -1)

    else:
        fit_window_start = max(term_start - fit_window_size, 1)
        fit_window_end = min(term_end + fit_window_size, bwFile.chroms(chromosome) -1)

    Wendepunkte = []
    fit_values = np.polynomial.polynomial.Polynomial.fit(range(fit_window_start, fit_window_end), bwFile.values(chromosome, fit_window_start, fit_window_end), deg=5, domain = [])
    firstDeriv = fit_values.deriv()
    secondDeriv = firstDeriv.deriv()
    thirdDeriv = secondDeriv.deriv()
    x = np.linspace(fit_window_start, fit_window_end, 100)
    #y = [np.polyval(rawValues.coef[::-1], i) for i in x]
    # plt.plot(x, rawValues(x))
    # plt.plot(x, firstDeriv(x))
    # plt.plot(x, secondDeriv(x))
    # plt.plot(x, thirdDeriv(x))
    plt.plot(term_start, fit_values(term_start) ,'o', label= "Term start", color = "yellow")
    plt.plot(term_end, fit_values(term_end) ,'o' , label= "Term End", color = '#49257B')

    roots = poly.polyroots(secondDeriv.coef)
    secondDeriv_real_roots = roots
    secondDeriv_real_roots = roots.real[abs(roots.imag) < 1e-5]
    secondDeriv_real_roots = list(filter(lambda r: r >= fit_window_start and r <= fit_window_end, secondDeriv_real_roots))
    #print("2nd Deriv roots: ", secondDeriv_real_roots)

    #print(poly.polyroots(secondDeriv.coef[::-1]))
    for root in secondDeriv_real_roots:
        #if abs(thirdDeriv(root) > 1e-5):
            #print("drew", root)
            #print(secondDeriv(root))
        Wendepunkte.append(root)


    bestScore = -1
    score = -1

    if(len(Wendepunkte) <= 0):
        #print(f"Alarm, keine wendepunkte: {tss}" )
        return "NA", -1, -1
    else:
        INFO_wind_start = 0
        INFO_wind_end = 0

        for wp in Wendepunkte:
            wp = round(wp)
            if(strand == "+"):
                score_end = min(wp + scoring_window_size + scoring_window_offset,  bwFile.chroms(chromosome) -1)
                score_start = min(wp + scoring_window_offset, bwFile.chroms(chromosome) -2)
            else:
                score_start = max(wp - scoring_window_size - scoring_window_offset, 1)
                score_end = max(wp - scoring_window_offset, 2)


            score = score_area(score_start, score_end, bwFile, chromosome, noiseLvL, iqr)
            tmp_start = wp
            tmp_end = score_end
            if(score > bestScore):
                bestScore = score
                INFO_wind_start = tmp_start
                INFO_wind_end = tmp_end


    return bestScore if bestScore != -1 else "NA", INFO_wind_start, INFO_wind_end


def drop_score_area(term_start, term_end, window_gene_ratio, post_term_window_offset, post_term_window_size, min_expr_drop, tss, bwFile, chromosome, strand):

    if(strand == '+'):
        score_pre_term_start = int(term_start - ((term_start - tss) * window_gene_ratio))
        score_pre_term_end = term_start

        score_post_term_start = min(term_end + post_term_window_offset, bwFile.chroms(chromosome) -1)
        score_post_term_end = min(score_post_term_start + post_term_window_size, bwFile.chroms(chromosome))

    else:
        score_pre_term_start = term_end
        score_pre_term_end = int(term_end + ((tss - term_end) * window_gene_ratio))

        score_post_term_end = max(term_start - post_term_window_offset, 2)
        score_post_term_start = max(score_post_term_end - post_term_window_size, 1)

    #print("bounds: ", score_pre_term_start, score_pre_term_end, bwFile.chroms(chrom), scoredTerm.strand, scoredTerm.start, myTss, scoredTerm.end)
    pre_term_expr = np.quantile(bwFile.values(chromosome, score_pre_term_start, score_pre_term_end), 0.5)

    post_term_expr = np.quantile(bwFile.values(chromosome, score_post_term_start, score_post_term_end), 0.5)

    expression_ratio = post_term_expr / max(pre_term_expr, 0.00001)
    
    score = 1 if expression_ratio < min_expr_drop else 0

    return score, score_pre_term_start, score_pre_term_end, score_post_term_start, score_post_term_end, (pre_term_expr * min_expr_drop)


def score_genome_region(forward_bigwig_path, 
                      reverse_bigwig_path, annotation_path, 
                      gffRhoterm, gffNocornac, master_table_path):
    

    MINIMUM_GENE_LENGTH = 15

    #Params avg scoring
    AVG_SCORE_WINDOW_SIZE = 25
    AVGSCORING_WINDOW_OFFSET_RHO = 150
    AVGSCORING_WINDOW_OFFSET_INTRINSIC = 20

    #Params deriv scoring
    FIT_WINDOW_HALF_SIZE = 75
    DERIV_SCORE_WINDOW_SIZE = 25
    DERIV_SCORE_WINDOW_OFFSET_RHO = 150
    DERIV_SCORE_WINDOW_OFFSET_INTRINSIC = 20


    #Params drop scoring
    #the ration of the evaluated pre-terminator window as a ratio to the gene length
    PRE_TERM_WINDOW_LENGTH_RATIO = 0.25
    #The min. ratio to which the post term window needs to be expressed compared to the pre term region
    MIN_EXPRESSION_RATIO = 0.25
    POST_TERM_WINDOW_SIZE = 50

    #offsets from downstream point of terminator
    POST_TERM_WINDOW_OFFSET_INTRINSIC = 150
    POST_TERM_WINDOW_OFFSET_RHO = 20

    fwbw = bigwig.open(forward_bigwig_path)
    rvbw = bigwig.open(reverse_bigwig_path)

    forward_noise, forward_iqr = CalcbackgroundNoise.InverseOfMasterTableNoise(annotation_path, fwbw, master_table_path)
    reverse_noise, reverse_iqr = CalcbackgroundNoise.InverseOfMasterTableNoise(annotation_path, rvbw, master_table_path)


    rhotermdata = gff3_parser.parse_gff3(gffRhoterm, verbose = False, parse_attributes = False)
    nocornacdata = gff3_parser.parse_gff3(gffNocornac, verbose = False, parse_attributes = False)

    all_term_data = pd.concat([rhotermdata, nocornacdata], ignore_index = True)

    forward_chromosome = re.match(r'^[^\.]+', rhotermdata.loc[1, "Seqid"]).group(0)
    reverse_chromosome = re.match(r'^[^\.]+', rhotermdata.loc[1, "Seqid"]).group(0)

    #free memory with next garbage collection
    del rhotermdata
    del nocornacdata

    for chrom in fwbw.chroms():
        prefix = re.match(r'^[^\.]+', chrom).group(0)
        if(prefix == forward_chromosome):
            forward_chromosome = chrom
            break

    for chrom in rvbw.chroms():
        prefix = re.match(r'^[^\.]+', chrom).group(0)
        if(prefix == reverse_chromosome):
            reverse_chromosome = chrom
            break
    

    all_terms_intervalls = zip(all_term_data.loc[:, "Start"].astype(int), all_term_data.loc[:, "End"].astype(int), all_term_data.loc[:, "Score"],
                                all_term_data.loc[:, "Type"], all_term_data.loc[: ,"Strand"])
    

    MasterTable = pd.read_csv(master_table_path, sep="\t")

    positive_tss = sorted(MasterTable[MasterTable['SuperStrand'] == '+']['SuperPos'].tolist())
    negative_TSS = sorted(MasterTable[MasterTable["SuperStrand"] == '-']['SuperPos'].tolist())

    del MasterTable

    with open("AllTermScoringPlusInfo.csv", 'w', newline ='') as info, open("AllTermScoring.tsv", 'w', newline = "") as allTermScoring:
        infoWriter = csv.writer(info, delimiter=",")
        termWriter = csv.writer(allTermScoring, delimiter='\t')

        termWriter.writerow(["seqid", "myTSS", "type", "strand", "start", "end", "initalScore", "avgScore", "derivScore", "dropScore"])
        infoWriter.writerow(["seqid", "myTSS", "type", "strand", "start", "end", "initalScore", "avgScore",
                      "derivScore", "dropScore", "bgNoise", "derivScoreWindowStart", "derivScoreWindowEnd", "avgScoreWindowStart", "avgScoreWindowEnd",
                        "dropScorePreTermStart", "dropScorePreTermEnd", "dropScorePostTermStart","dropScorePostTermEnd", "minExprDrop"])
        for start, end, inital_score, region, strand in all_terms_intervalls:
            if(region != 'terminator' and region != 'RhoTerminator'):
                continue
            
            #Get TSS
            if(strand == '+'):
                bw = fwbw
                chromosome = forward_chromosome
                noiseLvL = forward_noise
                iqr = forward_iqr
                myTss = FindFirstUpstreamTSS(start, MINIMUM_GENE_LENGTH, positive_tss, '+')
                #No TSS found, no gene here, we assume
                if(myTss <= -1):
                    continue
                estGeneExprMed = np.quantile(bw.values(chromosome, myTss, end), 0.5)
            else:
                bw = rvbw
                chromosome = reverse_chromosome
                noiseLvL = reverse_noise
                iqr = reverse_iqr
                myTss = FindFirstUpstreamTSS(end, MINIMUM_GENE_LENGTH, negative_TSS, '-')
                print(f"Start: {start} End: {end} tss: {myTss}")
                #No TSS found, no gene here, we assume
                if(myTss <= -1):
                    continue
                print("chromosome: ", chromosome)
                print("length: ", bw.chroms(chromosome))

                values = bw.values(chromosome, start, myTss)
                estGeneExprMed = np.quantile(values, 0.5)

            #If gene is not expressed, don't score at all
            if(estGeneExprMed <= noiseLvL):
                termWriter.writerow([chromosome, myTss, region, strand, start, end, inital_score, "NA", "NA", "NA"])
                infoWriter.writerow([chromosome, myTss, region, strand, start, end, inital_score, "NA", "NA", "NA", noiseLvL, -1, -1, -1, -1, -1, -1, -1, -1, -1])
                continue
            
            if(region == "terminator"):
                avg_offset = AVGSCORING_WINDOW_OFFSET_INTRINSIC
                deriv_offset = DERIV_SCORE_WINDOW_OFFSET_INTRINSIC
                drop_post_offset = POST_TERM_WINDOW_OFFSET_INTRINSIC
            else:
                avg_offset = AVGSCORING_WINDOW_OFFSET_RHO
                deriv_offset = DERIV_SCORE_WINDOW_OFFSET_RHO
                drop_post_offset = POST_TERM_WINDOW_OFFSET_RHO

            #Get avg score
            avg_score, avg_start, avg_end = avg_score_area(term_start=start, term_end=end, window_size=AVG_SCORE_WINDOW_SIZE, window_offset=avg_offset, 
                           chromosome=chromosome, strand=strand, bwFile=bw, noiseLvL=noiseLvL, iqr=iqr)

            #get deriv score
            deriv_score, deriv_start, deriv_end = deriv_score_area(term_start=start,  term_end=end, fit_window_size=FIT_WINDOW_HALF_SIZE,
                                                                   scoring_window_size= DERIV_SCORE_WINDOW_SIZE, 
                                                                   scoring_window_offset=deriv_offset, chromosome=chromosome, 
                                                                   strand=strand, bwFile=bw,noiseLvL=noiseLvL, iqr=iqr)

            #get drop score
            drop_score, drop_pre_start, drop_pre_end, drop_post_start, drop_post_end, drop_min_drop = drop_score_area(
                term_start=start, term_end=end, window_gene_ratio= PRE_TERM_WINDOW_LENGTH_RATIO, post_term_window_offset= drop_post_offset,
                post_term_window_size=POST_TERM_WINDOW_SIZE, min_expr_drop=MIN_EXPRESSION_RATIO, tss=myTss, bwFile=bw, chromosome=chromosome, strand=strand)

            #write to info and non info tsv
            toWrite = [chromosome, myTss, region, strand, start, end, inital_score, avg_score, deriv_score, drop_score]

            termWriter.writerow(toWrite)

            toWrite.extend([noiseLvL, deriv_start, deriv_end, avg_start, avg_end, drop_pre_start, drop_pre_end, drop_post_start, drop_post_end, drop_min_drop])

            infoWriter.writerow(toWrite)


if __name__ == "__main__":
    fwbw = sys.argv[1]
    rvbw = sys.argv[2]
    annot = sys.argv[3]
    gffRhoterm = sys.argv[4]
    gffNocornac= sys.argv[5]
    masterTable = sys.argv[6]
    score_genome_region(forward_bigwig_path=fwbw, reverse_bigwig_path=rvbw, annotation_path=annot, gffRhoterm=gffRhoterm, gffNocornac=gffNocornac, master_table_path=masterTable)
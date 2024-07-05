import glob
import math
import re
import pyBigWig as bigwig
import numpy as np
import AvgScoring
import CalcbackgroundNoise
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy.polynomial.polynomial as poly
import sys
import csv


def ScoreArea(WindowOffsetFromEnd, WindowSize, startOfArea, scoredTerm, bwFile, noiseLvL):
    if(scoredTerm.strand == '-'):
        WindowSize *= -1
        WindowOffsetFromEnd *= 1
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

def DerivScroring(bigWigPath, annotationPath, gffRhoterm, gffNocornac, MasterTablePath):    
    #PARAMS
    SearchWindow = 75
    scoring_window_offset_rho = 150
    scoring_window_offset_intrinsic = 20

    scoring_window_size = 25    
    bw = bigwig.open(bigWigPath)

    noiseLVL = CalcbackgroundNoise.CalcBackgroundNoise(annotationPath, bigWigPath)
    TSSTermPairings = AvgScoring.AvgScoreTerminators(gffRhoterm, gffNocornac, bigWigPath, MasterTablePath, annotationPath)

    chromo = re.match('^[^\.]+', TSSTermPairings[list(TSSTermPairings.keys())[0]][0].seqid).group(0)
    for chrom in bw.chroms():
        prefix = re.match('^[^\.]+', chrom).group(0)
        if(prefix == chromo):
            chromo = chrom
            break
    
    inverseGFF = CalcbackgroundNoise.GetInverseOfGFF(annotationPath, bw.chroms(chromo))



    #testf = np.polynomial.Polynomial.fit(range(1,200,1), bw.values(chromo, 1, 200), deg = 10)

    # x = range(1,200,1)
    # plt.plot(x, testf(x))
    # roots = poly.polyroots(testf.coef)
    # real_roots =roots.real[abs(roots.imag) < 1e-5]
    # real_roots = list(filter(lambda r: r > 0, real_roots))
    # print(real_roots)
    #plt.figure()
    fig, ax = plt.subplots()

    data =[]
    data.append(["seqid", "myTSS", "type", "strand", "start", "end", "initalScore", "avgScore", "derivScore"])

    ax.boxplot(bw.values(chromo, 1, bw.chroms(chromo)), showfliers= False)
    ax.axhline(y =  noiseLVL, color= "r", linewidth = 1)

    plt.savefig("Boxplot.svg", format = 'svg', dpi=300)
    plt.figure()

    plt.plot(bw.values(chromo, 1 , bw.chroms(chromo)), color = '#DD9837')

    print(f"noiseLVL deriv: {noiseLVL}")
    for start, end in inverseGFF:
        x = np.linspace(start, end, 100)
        #plt.plot(x, -10 * np.ones_like(x) )

    for tss, scoredterms in TSSTermPairings.items():
        plt.plot(tss, 0 , marker = 'o', color = 'black')
        for scoredterm in scoredterms:            
            
            if(scoredterm.strand == '+'):
                estGeneExprMed = np.quantile(bw.values(scoredterm.seqid, tss, end), 0.5)
            else:
                estGeneExprMed = np.quantile(bw.values(scoredterm.seqid, scoredterm.start, tss), 0.5)
            
            #if gene is not expressed, we cannot judge its terminators
            if(estGeneExprMed <= noiseLVL):
                scoredterm.derivScore = "NA"
                data.append([f"{scoredterm.seqid}", f"{tss}", scoredterm.type, scoredterm.strand, scoredterm.start, scoredterm.end, scoredterm.initialScore, scoredterm.avgScore, scoredterm.derivScore])
                continue

            fit_window_start = scoredterm.start - SearchWindow
            fit_window_end = scoredterm.end + SearchWindow

            Wendepunkte = []
            rawValues = np.polynomial.polynomial.Polynomial.fit(range(fit_window_start, fit_window_end), bw.values(chromo, fit_window_start, fit_window_end), deg=5, domain = [])
            #print("Range:" , range(fit_window_start, fit_window_end))
            #print(" normal poly: ",(rawValues.coef))
            firstDeriv = rawValues.deriv()
            secondDeriv = firstDeriv.deriv()
            thirdDeriv = secondDeriv.deriv()
            x = np.linspace(fit_window_start, fit_window_end, 100)
            #y = [np.polyval(rawValues.coef[::-1], i) for i in x]
            # plt.plot(x, rawValues(x))
            # plt.plot(x, firstDeriv(x))
            # plt.plot(x, secondDeriv(x))
            # plt.plot(x, thirdDeriv(x))
            plt.plot(scoredterm.start, rawValues(scoredterm.start) ,'o', label= "Term start", color = '#49257B')
            plt.plot(scoredterm.end, rawValues(scoredterm.end) ,'o' , label= "Term End", color = '#49257B')

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
                plt.plot(root, rawValues(root), marker = "o", color ='yellow')
            
            bestScore = -1
            score = -1

            if(len(Wendepunkte) <= 0):
                print(f"Alarm, keine wendepunkte: {tss}" )
                #scoredterm.derivScore = "NA"

            for wp in Wendepunkte:
                wp = round(wp)
                if(scoredterm.type == "terminator"):
                    score = ScoreArea(scoring_window_offset_intrinsic, scoring_window_size, wp, scoredterm, bw, noiseLVL)
                else:
                    score = ScoreArea(scoring_window_offset_rho, scoring_window_size, wp, scoredterm, bw, noiseLVL)

                if(score > bestScore):
                    bestScore = score
            
            scoredterm.derivScore = bestScore

            data.append([f"{scoredterm.seqid}", f"{tss}", scoredterm.type, scoredterm.strand, scoredterm.start, scoredterm.end, scoredterm.initialScore, scoredterm.avgScore, scoredterm.derivScore])



    with open("AllTermScoring.tsv", 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(data)

    plt.savefig("Wiggle.svg", format = 'svg', dpi=300)
    plt.show()

if __name__ == "__main__":
    bw = sys.argv[1]
    annot = sys.argv[2]
    gffRhoterm = sys.argv[3]
    gffNocornac= sys.argv[4]
    masterTable = sys.argv[5]
    DerivScroring(bw, annot, gffRhoterm, gffNocornac, masterTable)


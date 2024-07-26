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


def ScoreArea(WindowOffsetFromEnd, WindowSize, startOfArea, scoredTerm, bwFile, noiseLvL, iqr):
    if(scoredTerm.strand == '-'):
        WindowSize *= -1
        WindowOffsetFromEnd *= -1
        windStart = max(startOfArea + WindowOffsetFromEnd, 1)
        windEnd = max(windStart + WindowSize, 1)
    else:
        windStart = min(startOfArea + WindowOffsetFromEnd, bwFile.chroms(scoredTerm.seqid))
        windEnd = min(windStart + WindowSize, bwFile.chroms(scoredTerm.seqid))
    
    if(windStart > windEnd):
        b = windEnd
        windEnd = windStart
        windStart = b

    postTermExprQ = np.quantile(bwFile.values(scoredTerm.seqid, windStart, windEnd), 1)

    print(iqr)
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

def DerivScroring(forward_bigwig_path, reverse_bigwig_path, annotationPath, gffRhoterm, gffNocornac, MasterTablePath):    
    #PARAMS
    SearchWindow = 75
    scoring_window_offset_rho = 150
    scoring_window_offset_intrinsic = 20

    scoring_window_size = 25    
    fwbw = bigwig.open(forward_bigwig_path)
    rvbw = bigwig.open(reverse_bigwig_path)

    forward_noise, forward_iqr = CalcbackgroundNoise.InverseOfMasterTableNoise(annotationPath, forward_bigwig_path, MasterTablePath)
    reverse_noise, reverse_iqr = CalcbackgroundNoise.InverseOfMasterTableNoise(annotationPath, reverse_bigwig_path, MasterTablePath)

    TSSTermPairings = AvgScoring.AvgScoreTerminators(gffRhoterm=gffRhoterm,
                                                    gffNocornac= gffNocornac,
                                                    forward_bigwig_path=forward_bigwig_path,
                                                    reverse_bigwig_path=reverse_bigwig_path,
                                                    master_table_path=MasterTablePath,
                                                    annotgff=annotationPath)

    
    forward_chromo = re.match('^[^\.]+', TSSTermPairings[list(TSSTermPairings.keys())[0]][0].seqid).group(0)
    for chrom in fwbw.chroms():
        prefix = re.match('^[^\.]+', chrom).group(0)
        if(prefix == forward_chromo):
            forward_chromo = chrom
            break

    reverse_chromo = re.match('^[^\.]+', TSSTermPairings[list(TSSTermPairings.keys())[0]][0].seqid).group(0)
    for chrom in rvbw.chroms():
        prefix = re.match('^[^\.]+', chrom).group(0)
        if(prefix == reverse_chromo):
            reverse_chromo = chrom
            break
    
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

    info_data = []
    info_data.append(["seqid", "myTSS", "type", "strand", "start", "end", "initalScore", "avgScore", 
                      "derivScore", "bgNoise", "derivScoreWindowStart", "derivScoreWindowEnd", "avgScoreWindowStart", "avgScoreWindowEnd"])

    ax.boxplot(fwbw.values(forward_chromo, 1, fwbw.chroms(forward_chromo)), showfliers= False)
    ax.axhline(y = forward_noise, color= "r", linewidth = 1)
    #ax.axhline(y =  reverse_noise, color= "r", linewidth = 1)

    plt.savefig("Boxplot.svg", format = 'svg', dpi=300)
    plt.figure()

    plt.plot(fwbw.values(forward_chromo, 1 , fwbw.chroms(forward_chromo)), color = '#DD9837')
    plt.plot(rvbw.values(reverse_chromo, 1 , rvbw.chroms(reverse_chromo)), color = '#2267c8')

    print(f"noiseLVL deriv: {forward_noise}")

    for tss, scoredterms in TSSTermPairings.items():
        plt.plot(tss, 0 , marker = 'o', color = 'black')
        for scoredterm in scoredterms:            
            
            if(scoredterm.strand == '+'):
                bw = fwbw
                chromo = forward_chromo
                noiseLVL = forward_noise
                iqr = forward_iqr
                estGeneExprMed = np.quantile(bw.values(scoredterm.seqid, tss, scoredterm.end), 0.5)

                fit_window_start = scoredterm.start - SearchWindow
                fit_window_end = scoredterm.end + SearchWindow
            else:
                bw = rvbw
                chromo = reverse_chromo
                noiseLVL = reverse_noise
                iqr = reverse_iqr
                estGeneExprMed = np.quantile(bw.values(scoredterm.seqid, scoredterm.start, tss), 0.5)

                fit_window_start = scoredterm.start - SearchWindow
                fit_window_end = scoredterm.end + SearchWindow
            
            #if gene is not expressed, we cannot judge its terminators
            if(estGeneExprMed <= noiseLVL):
                scoredterm.derivScore = "NA"
                data.append([f"{scoredterm.seqid}", f"{tss}", scoredterm.type, scoredterm.strand, scoredterm.start, scoredterm.end, scoredterm.initialScore, scoredterm.avgScore, scoredterm.derivScore])
                continue          

            Wendepunkte = []
            fitValues = np.polynomial.polynomial.Polynomial.fit(range(fit_window_start, fit_window_end), bw.values(chromo, fit_window_start, fit_window_end), deg=5, domain = [])
            #print("Range:" , range(fit_window_start, fit_window_end))
            #print(" normal poly: ",(rawValues.coef))
            firstDeriv = fitValues.deriv()
            secondDeriv = firstDeriv.deriv()
            thirdDeriv = secondDeriv.deriv()
            x = np.linspace(fit_window_start, fit_window_end, 100)
            #y = [np.polyval(rawValues.coef[::-1], i) for i in x]
            # plt.plot(x, rawValues(x))
            # plt.plot(x, firstDeriv(x))
            # plt.plot(x, secondDeriv(x))
            # plt.plot(x, thirdDeriv(x))
            plt.plot(scoredterm.start, fitValues(scoredterm.start) ,'o', label= "Term start", color = "yellow")
            plt.plot(scoredterm.end, fitValues(scoredterm.end) ,'o' , label= "Term End", color = '#49257B')

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
                #plt.plot(root, fitValues(root), marker = "o", color ='yellow')
            
            bestScore = -1
            score = -1

            if(len(Wendepunkte) <= 0):
                print(f"Alarm, keine wendepunkte: {tss}" )
                scoredterm.derivScore = "NA"

            INFO_wind_start = 0
            INFO_wind_end = 0

            for wp in Wendepunkte:
                wp = round(wp)
                if(scoredterm.type == "terminator"):
                    score, tmp_start, tmp_end = ScoreArea(scoring_window_offset_intrinsic, scoring_window_size, wp, scoredterm, bw, noiseLVL, iqr)
                else:
                    score, tmp_start, tmp_end = ScoreArea(scoring_window_offset_rho, scoring_window_size, wp, scoredterm, bw, noiseLVL, iqr)

                if(score > bestScore):
                    bestScore = score
                    INFO_wind_start = tmp_start
                    INFO_wind_end = tmp_end
            
            plt.xlim(tss, INFO_wind_end + 50)
            if len(Wendepunkte) > 0:
                plt.ylim(-5,  np.quantile(bw.values(chromo, INFO_wind_start, INFO_wind_end), 1) + 10 * 2)
            else:
                if(scoredterm.strand == '+'):
                    plt.ylim(-5,  np.quantile(bw.values(chromo, tss, scoredterm.end), 1) + 10 * 2)
                else:
                    plt.ylim(-5,  np.quantile(bw.values(chromo, scoredterm.start, tss), 1) + 10 * 2)
            plt.plot(x, fitValues(x))
            
            if scoredterm.strand == '+':
                plt.plot(fwbw.values(forward_chromo, 1 , fwbw.chroms(forward_chromo)), color = '#DD9837')
            else:
                plt.plot(rvbw.values(reverse_chromo, 1 , rvbw.chroms(reverse_chromo)), color = '#2267c8')
                
            plt.legend()
            plt.xlabel('Position')
            plt.ylabel('Value')
            plt.title(f'Plot for Terminator {scoredterm.start} {scoredterm.end}')
            plt.savefig(f'plot_TSS_{tss}.png')  # Save the plot as a PNG file
            plt.close()

            scoredterm.derivScore = bestScore

            data.append([f"{scoredterm.seqid}", 
                         f"{tss}", scoredterm.type, 
                         scoredterm.strand, scoredterm.start, 
                         scoredterm.end, scoredterm.initialScore, 
                         scoredterm.avgScore, 
                         scoredterm.derivScore])
            
            info_data.append([f"{scoredterm.seqid}", 
                         f"{tss}", scoredterm.type, 
                         scoredterm.strand, scoredterm.start, 
                         scoredterm.end, scoredterm.initialScore, 
                         scoredterm.avgScore, 
                         scoredterm.derivScore, noiseLVL ,INFO_wind_start, INFO_wind_end, scoredterm.avgScoreStart, scoredterm.avgScoreEnd])

    with open("AllTermScoringPlusInfo.csv", 'w', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerows(info_data)

    with open("AllTermScoring.tsv", 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(data)

    plt.savefig("Wiggle.svg", format = 'svg', dpi=300)
    plt.show()

if __name__ == "__main__":
    fwbw = sys.argv[1]
    rvbw = sys.argv[2]
    annot = sys.argv[3]
    gffRhoterm = sys.argv[4]
    gffNocornac= sys.argv[5]
    masterTable = sys.argv[6]
    DerivScroring(forward_bigwig_path=fwbw, reverse_bigwig_path=rvbw, annotationPath=annot, gffRhoterm=gffRhoterm, gffNocornac=gffNocornac, MasterTablePath=masterTable)


import math
import pyBigWig as bigwig
import numpy as np
import AvgScoring
import CalcbackgroundNoise
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy.polynomial.polynomial as poly
import sys

def DerivScroring(bigWigPath, annotationPath, terminatorGFFs, MasterTablePath):

    # def deriv(bigwigValues, degree):

    #     if(degree == 0):
    #         return bigwigValues
        
    #     return deriv(np.gradient(bigwigValues, 1), degree-1)
            
    def LinearInterWig(x):
        bigWigValues = bw.values(chromo, 1 , bw.chroms(chromo))

        if(x > len(bigWigValues)):
            return bigWigValues[-1]
        elif(x < 1):
            return bigWigValues[0]

        start = math.floor(x)
        end = math.ceil(x)
        interp = x - start
        return (bigWigValues[start] + ((bigWigValues[end] - bigWigValues[start]) * interp))
    

    bw = bigwig.open(bigWigPath)
    noiseLVL = CalcbackgroundNoise.CalcBackgroundNoise(annotationPath, bigWigPath)

    TSSTermPairings = AvgScoring.AvgScoreTerminators(terminatorGFF, bigWigPath, MasterTablePath, annotationPath)

    chromo = "NC_004703.1"

    inverseGFF = CalcbackgroundNoise.GetInverseOfGFF(annotationPath, bw.chroms(chromo))


    SearchWindow = 75
    Wendepunkte = []
    drawIdx = 0
    getIdx = 0            

    #testf = np.polynomial.Polynomial.fit(range(1,200,1), bw.values(chromo, 1, 200), deg = 10)

    # x = range(1,200,1)
    # plt.plot(x, testf(x))
    # roots = poly.polyroots(testf.coef)
    # real_roots =roots.real[abs(roots.imag) < 1e-5]
    # real_roots = list(filter(lambda r: r > 0, real_roots))
    # print(real_roots)
    #plt.figure()
    fig, ax = plt.subplots()

    for chromo in bw.chroms():
        ax.boxplot(bw.values(chromo, 1, bw.chroms(chromo)), showfliers= False)
        ax.axhline(y =  noiseLVL, color= "r", linewidth = 1)

        plt.savefig("Boxplot.svg", format = 'svg', dpi=300)
        plt.figure()

        plt.plot(bw.values(chromo, 1 , bw.chroms(chromo)), color = '#DD9837')

        for start, end in inverseGFF:
            x = np.linspace(start, end, 100)
            #plt.plot(x, -10 * np.ones_like(x) )

        for tss, scoredterms in TSSTermPairings.items():
            plt.plot(tss, 0 , marker = 'o', color = 'black')
            for scoredterm in scoredterms:

                rawValues = np.polynomial.polynomial.Polynomial.fit(range(scoredterm.start - SearchWindow, scoredterm.end + SearchWindow +1), bw.values(chromo, scoredterm.start - SearchWindow, scoredterm.end + SearchWindow +1), deg=5, domain = [])
                #print("Range:" , range(scoredterm.start - SearchWindow, scoredterm.end + SearchWindow +1))
                #print(" normal poly: ",(rawValues.coef))
                firstDeriv = rawValues.deriv()
                secondDeriv = firstDeriv.deriv()
                thirdDeriv = secondDeriv.deriv()
                x = np.linspace(scoredterm.start - SearchWindow, scoredterm.end + SearchWindow +1, 100)
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
                secondDeriv_real_roots = list(filter(lambda r: r >= scoredterm.start - SearchWindow and r <= scoredterm.end + SearchWindow +1, secondDeriv_real_roots))
                #print("2nd Deriv roots: ", secondDeriv_real_roots)

                #print(poly.polyroots(secondDeriv.coef[::-1]))
                for root in secondDeriv_real_roots:
                    if abs(thirdDeriv(root) > 1e-5):
                        #print("drew", root)
                        #print(secondDeriv(root))
                        Wendepunkte.append(root)
                        plt.plot(root, rawValues(root), marker = "o", color ='yellow')





    plt.savefig("Wiggle.svg", format = 'svg', dpi=300)
    plt.show()

if __name__ == "__main__":
    bw = sys.argv[1]
    annot = sys.argv[2]
    termGFF = sys.argv[3]
    masterTable = sys.argv[4]
    DerivScroring(bw, annot, termGFF, masterTable)


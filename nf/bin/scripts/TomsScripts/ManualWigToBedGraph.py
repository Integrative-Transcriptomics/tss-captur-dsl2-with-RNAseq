import wiggelen as wigg
import wiggelen.intervals as inter
import wiggelen.merge
import FastaReader
import sys
import MergeWigs

def WigToBedGraph(Wigglepath, chromSizesFilePath, outputPath):

    #bedfile = open("/testenv/TestWigs/ManualBED.bed", "x")
    #correctWalker = wigg.fill(wigg.walk(open(path), filler = 0))
    #correctWalker = wigg.fill(wigg.walk(open(Wigglepath)), regions = None, filler = 0, only_edges= False)
    correctWalker = MergeWigs.merge_wigs(Wigglepath)
    #wigg.write(inter.coverage(correctWalker), name='My example')
    BedGraphFile = open(f'{outputPath}/autoBG.bedGraph', 'w')

    sizes = open(chromSizesFilePath, 'r')

    lastVal = 0
    lastReg = ""
    lastPos = 0
    first = True
    start = 0
    totalsum = 0

    sizeDict = {}

    for line in sizes:
        parts = line.split()
        sizeDict[parts[0]] = int(parts[1])

    trueCurrRegLength = 0
    for reg, pos, val in correctWalker:
        if(first):
            lastPos = pos
            lastReg = reg
            lastVal = val
            start = pos
            first = False
            if(pos != 1):
                BedGraphFile.write(f'{lastReg}\t1\t{pos}\t0\n')
            #trueCurrRegLength = FastaReader.get_fasta_seq_length(f"/home/GitRepos/Bachelorarbeit/.venv/dataset/{reg}.fna")
            trueCurrRegLength = sizeDict[reg]

        #If we are entering a new region
        if(reg != lastReg):
            #value might not have changed but still need to fill the current tracked intervall
            BedGraphFile.write(f'{lastReg}\t{start}\t{lastPos}\t{lastVal}\n')
            
            #If true chromosome size > lastPos, there are 0 at end of region which were cut off from .wig
            if(trueCurrRegLength > lastPos):
                BedGraphFile.write(f'{lastReg}\t{lastPos}\t{trueCurrRegLength}\t0\n')
            
            #if next chromosome does not start with position 1, fill missing poses with 0
            if(pos != 1):
                BedGraphFile.write(f'{reg}\t1\t{pos}\t0\n')

            end = pos
            lastVal = val
            lastReg = reg
            start = pos
            #trueCurrRegLength = FastaReader.get_fasta_seq_length(f"/home/GitRepos/Bachelorarbeit/.venv/dataset/{reg}.fna")
            trueCurrRegLength = sizeDict[reg]
            
        #if we still are in the same region
        elif(val != lastVal):        
            BedGraphFile.write(f'{lastReg}\t{start}\t{pos}\t{lastVal}\n')

            end = pos
            lastVal = val
            lastReg = reg
            start = pos

        lastPos = pos
        
        totalsum += val

    #close final intervall of final chromosome
    BedGraphFile.write(f'{reg}\t{start}\t{pos}\t{lastVal}\n')

    #fill potentially end-cut off 0s from chromosome
    if(trueCurrRegLength > lastPos):
        BedGraphFile.write(f'{lastReg}\t{lastPos}\t{trueCurrRegLength}\t0\n')

    #print(pow(totalsum, 1))
    print(BedGraphFile)
    return BedGraphFile
    # BedGraphFile.close()


if __name__ == '__main__':
    WigToBedGraph(sys.argv[1], sys.argv[2], sys.argv[3])
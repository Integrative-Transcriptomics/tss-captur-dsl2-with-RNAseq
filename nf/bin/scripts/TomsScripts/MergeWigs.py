import glob
import wiggelen
import wiggelen.merge
import sys

def merge_wigs(wigsPath):

    walkers = []
    for wigPath in glob.glob(f'{wigsPath}/*.wig'):
        walkers.append(wiggelen.fill(wiggelen.walk(open(wigPath)), regions = None, filler = 0, only_edges= False))
    
    #merged = walkers[0]
    merged = wiggelen.merge.merge(*walkers, merger = wiggelen.merge.mergers['mean'])
    #for i in range(1, len(walkers)):
       # merged = wiggelen.merge.merge(merged, walkers[i], merger= wiggelen.merge.mergers['sum'])

    return merged
    # mergedWigFile = open("MergedWigs.wig", "w")
    # for line in wiggelen.walk(merged):
    #     mergedWigFile.write(line)


if __name__ == "__main__":
    wigsPath = sys.argv[1]
    merge_wigs(wigsPath=wigsPath)
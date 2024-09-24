"""Merger for .wig files

This script merges .wig files located at the specified directory. 

This script requires `wiggelen` to be installed. 
"""


import glob
import wiggelen
import wiggelen.merge
import sys

def merge_wigs(wigsPath):

    walkers = []
    for wigPath in glob.glob(f'{wigsPath}/*.wig'):
        walkers.append(wiggelen.fill(wiggelen.walk(open(wigPath)), regions = None, filler = 0, only_edges= False))
    
    merged = wiggelen.merge.merge(*walkers, merger = wiggelen.merge.mergers['mean'])

    return merged


if __name__ == "__main__":
    wigsPath = sys.argv[1]
    merge_wigs(wigsPath=wigsPath)
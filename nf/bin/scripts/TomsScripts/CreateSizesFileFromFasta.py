"""Helper script for getting size of a genome in a fasta file
"""


import glob
import sys
import FastaReader

def CreateSizesFile(fastaPath, outputSizesPath):
    chromosomeName = 'No chromosome found'

    extensions = ["fa", "fna", "fasta", "frn", "faa", "ffn"]

    files = []
    for ex in extensions:
        files.extend(glob.glob(f'{fastaPath}/*.{ex}'))

    f = open(files[0], 'r')

    for line in f:
        if line.startswith('>'):
            chromosomeName = line.split()[0][1:]
    
    length = FastaReader.get_fasta_seq_length(fastaPath=fastaPath)

    sizesFile = open(f'{outputSizesPath}/autoSizeFile.sizes', 'w')

    sizesFile.write(f'{chromosomeName} {length}\n')

    print(sizesFile)
    return sizesFile

if __name__ == "__main__":
    fastaPath = sys.argv[1]
    outpath = sys.argv[2]
    CreateSizesFile(fastaPath=fastaPath, outputSizesPath=outpath)
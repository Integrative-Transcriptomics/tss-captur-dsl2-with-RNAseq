import sys
import FastaReader

def CreateSizesFile(fastaPath):
    chromosomeName = 'No chromosome found'

    f = open(fastaPath, 'r')
    for line in f:
        if line.startswith('>'):
            chromosomeName = line.split()[0][1:]
    
    length = FastaReader.get_fasta_seq_length(fastaPath=fastaPath)

    sizesFile = open(f'{chromosomeName}.sizes', 'w')

    sizesFile.write(f'{chromosomeName} {length}\n')

    print(sizesFile)
    return sizesFile

if __name__ == "__main__":
    fastaPath = sys.argv[1]
    CreateSizesFile(fastaPath=fastaPath)
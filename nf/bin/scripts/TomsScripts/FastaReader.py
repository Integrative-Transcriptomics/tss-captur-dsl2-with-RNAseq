"""Helper script for getting the full length of a fasta file

Used for determining the true length of genomes, contained within .wig files.
"""


import glob

def get_fasta_seq_length(fastaPath):
    sequence_length = 0
    extensions = ["fa", "fna", "fasta", "frn", "faa", "ffn"]

    files = []
    for ex in extensions:
        files.extend(glob.glob(f'{fastaPath}/*.{ex}'))

    f = open(files[0], 'r')
    
    for line in f:
        if line.startswith('>'):
            continue  # Skip header lines
        sequence_length += len(line.strip())  # Add length of non-empty sequence lines

    return sequence_length
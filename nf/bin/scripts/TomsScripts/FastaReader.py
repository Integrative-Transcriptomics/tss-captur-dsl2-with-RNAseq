
def get_fasta_seq_length(fastaPath):
    sequence_length = 0
    f = open(fastaPath, 'r')
    for line in f:
        if line.startswith('>'):
            continue  # Skip header lines
        sequence_length += len(line.strip())  # Add length of non-empty sequence lines

    return sequence_length
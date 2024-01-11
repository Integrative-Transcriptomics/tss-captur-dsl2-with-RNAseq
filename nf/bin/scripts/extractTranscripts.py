"""Extraction of transcripts classified as RNAs

This script extracts all transcripts classified as RNA with the newly
computed coordinates. It accepts an .tsv file and produces a FASTA file.
This FASTA file can be then used as an input for RNAFold

This script requires `pandas` and `biopython` be installed within the Python
environment you are running this script in.

"""

from Bio import SeqIO, SeqRecord
import argparse
import re
import os.path
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rnas")
    parser.add_argument("--genome_path")
    args = parser.parse_args()
    genome_name = re.split(
        "/", re.split("_allocated_terminators", args.rnas)[0])[-1]

    filelist = os.listdir(args.genome_path)
    fasta_ext = "(fa|fna|fasta|faa|frn|ffn)"
    regex = re.compile("%s.*\.%s$" % (genome_name, fasta_ext))
    match = list(filter(regex.match, filelist))[0]
    file_path = os.path.join(args.genome_path, match)
    with open(file_path, "r") as handle:
        record = list(SeqIO.parse(handle, "fasta"))[0]
    terminator_df = pd.read_csv(args.rnas, index_col=None, sep="\t")
    # terminator_df = terminator_df.loc[terminator_df['type'] == "RNA"]
    terminator_df = terminator_df.loc[terminator_df['type'].isin(["RNA", "COD"])]

    terminator_df.iloc[:, 12:] = terminator_df.iloc[:, 12:].astype(int)
    transcript_list = []
    for index_tss, row in terminator_df.iterrows():
        t_id, start, end, strand = row[["transcript_id", "start_with_terminator",
                                        "end_with_terminator", "strand"]].values

        transcript = record.seq[start: end]
        if strand == "-":
            transcript = transcript.reverse_complement()
        seq = SeqRecord.SeqRecord(
            transcript, id=t_id, name=t_id, description="")
        transcript_list.append(seq)
    with open("%s_rna_transcripts.fasta" % (genome_name), "w") as output_handle:
        SeqIO.write(transcript_list, output_handle, "fasta")

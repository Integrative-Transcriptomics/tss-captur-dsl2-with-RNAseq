"""Transformation of RhoTermPredict output to GFF

This script modifies the output of RhoTermPredict to a GFF format. 
This script expects the .tsv-file of RhoTermPredict and the Genome's NCBI Identifier as an input.

This script requires that `pandas` be installed within the Python
environment you are running this script in.

"""

import pandas as pd
import argparse


def table_to_gff(df, genome_name):
    gff = df.reset_index().rename(
        columns={"Start RUT": "start", "End RUT": "end", "Region": "attributes"})
    gff["seqname"] = genome_name
    gff["frame"] = "."
    gff["source"] = "RhoTermPredict"
    gff["feature"] = "RhoTerminator"
    gff = gff[["seqname", "source", "feature", "start",
               "end", "score", "Strand", "frame", "attributes"]]
    return gff


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rhoterm")
    parser.add_argument("--genome")
    args = parser.parse_args()
    rhoterm_df = pd.read_csv(args.rhoterm, sep="\t")
    rhoterm_gff = table_to_gff(rhoterm_df, args.genome)
    rhoterm_gff.to_csv("%s_rhoterm.gff" % args.genome, sep="\t",
                       header=False, index=False)

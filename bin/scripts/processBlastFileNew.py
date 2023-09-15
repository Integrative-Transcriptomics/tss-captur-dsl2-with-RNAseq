"""BLAST result evaluation

This script evaluates the BLAST-hits and returns the best pairwise alignment
for each transcript. 

For the choice of the functions, please refer to the original thesis. 

This script requires ´pandas´, ´numpy´ and ´scipy´ to be installed. 

"""

import pandas as pd
import numpy as np
import argparse
import os
import re
import time
from scipy.stats import norm

pd.options.mode.chained_assignment = None  # default='warn'


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--threshold", type=int,
                        help="Threshold for allowing a mapped TSS", default=10)
    parser.add_argument(
        "blast_results",  nargs='+')
    args = parser.parse_args()

    filtered_files = [f for f in args.blast_results if is_non_zero_file(f)]
    name_files = [os.path.splitext(f)[0] for f in filtered_files]
    df_filtered_files = [pd.read_csv(
        f, index_col=None, sep="\t", header=None) for f in filtered_files]

    colnames = ["qseqid", "qstart", "qend", "qseq", "sseqid",
                "sstart", "send", "sseq", "evalue", "bitscore", "pident", "frames", "qcovhsp"]
    df_queries = pd.concat(df_filtered_files, keys=name_files)
    df_queries.columns = colnames

    mapped_queries = df_queries["qseqid"].unique()

    max_norm_pdf = norm.pdf(0, 0, 8)

    # This commented code was the first attempt to run the 
    # label via similarity. 
    # Due to time restrictions, this was not included in TSS-captur
    # It works as a reminder of this function
    #
    #
    # df_high_similarity = df_queries[(df_queries.pident >=
    #                                  99) & (df_queries.qcovhsp >= 99)]

    # df_high_similarity["transcript_id"] = df_high_similarity["qseqid"].map(
    #     lambda x: re.findall("\|((?:orphan_|antisense_)\d+)\|", x)[0])
    # df_high_similarity["file_name"] = df_high_similarity.index.map(
    #     lambda x: x[0])
    # df_high_similarity["genome_name"] = df_high_similarity.file_name.map(
    #     lambda x: re.sub("_(?:orphan|antisense)_(?:plus|minus)", "", x))
    # genomes = df_high_similarity["genome_name"].unique()
    # df_high_similarity["qstrand"] = df_high_similarity["file_name"].map(
    #     lambda x: "-" if len(re.findall("minus", x)) > 0 else "+")
    # df_high_similarity["sstrand"] = df_high_similarity.qstrand.where(
    #     df_high_similarity.frames == "1/1", np.where(df_high_similarity.qstrand == "-", "+", "-"))
    # df_high_similarity = df_high_similarity[["transcript_id",
    #                                          "sseqid", "sstart", "send", "qstrand", "sstrand", "pident", "qcovhsp", "genome_name"]]
    # for g in genomes:
    #     df_high_similarity[df_high_similarity.genome_name == g].drop(
    #         "genome_name", axis=1).to_csv("%s_high_similar.tsv" % g, sep="\t", index=False)

    df_queries["eval_qstart"] = df_queries["qstart"].map(
        lambda y: 1 - ((min(max(0, y-args.threshold), args.threshold+20))/(args.threshold+20)))
    df_queries["eval_pident"] = df_queries["pident"].map(
        lambda x: (norm.pdf(max(-25, (x-75)), 0, 8)/max_norm_pdf))
    df_queries["eval_qcovhsp"] = df_queries["qcovhsp"].map(
        lambda z: min(1, (z - 50) / 25))

    df_queries["scored"] = 0.45*df_queries["eval_qstart"]+0.35 * \
        df_queries["eval_pident"]+0.20*df_queries["eval_qcovhsp"]

    df_queries_select_best = df_queries.sort_values(
        "scored", ascending=False, inplace=False).drop_duplicates("qseqid", inplace=False)

    separate_df = {idx: gp.xs(idx, level=0, axis=0)
                   for idx, gp in df_queries_select_best.groupby(level=0, axis=0)}

    for name, gp_df in separate_df.items():
        gp_df.to_csv("%s.tsv" % name, sep="\t", header=None, index=False)

"""CNIT converter to a GFF

This script modifies the output of CNIT to be similar to a GFF file. 
For this, it requires the CNIT-output. 

This script requires `pandas` and `numpy` to be installed. 
"""


import re
import pandas as pd
import numpy as np
import argparse


def create_columns(row):
    row_id_all = row["Transcript ID"]
    # To have the same annotation for both predictors.
    row["index"] = "COD" if row["index"] == "coding" else "RNA"
    # Separate the important information
    row["transcript_id"] = re.findall(
        "\|((?:orphan_|antisense_|internal_)\d+)\|", row_id_all)[0]
    row["position"] = int(re.findall("\|Start:(\d+)\|", row_id_all)[0])
    row["strand"] = re.findall("\|Strand:(\+|\-)", row_id_all)[0]
    # compute the posible coordinates depending on the strand
    if row["strand"] == "+":
        row["PredictionStartCorrected"] = row["position"] + row["start"]
        row["PredictionEndCorrected"] = row["position"] + row["end"]
        row["GenomeStart"] = row["position"]
        row["GenomeEnd"] = row["position"] + row[6]
    else:
        row["PredictionStartCorrected"] = row["position"] - row["end"]
        row["PredictionEndCorrected"] = row["position"] - row["start"]
        row["GenomeEnd"] = row["position"]
        row["GenomeStart"] = row["position"] - row[6]

    # Better ordering
    new_order = [7, 8, 9, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13]
    row = row[new_order]
    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnit_list", nargs="+",)
    args = parser.parse_args()
    genomes = np.unique(
        [re.split("/", re.split("_antisense|_orphan|_internal", x)[0])[-1]for x in args.cnit_list])
    print(args.cnit_list)
    for g in genomes:
        # Return one file per genome, i.e. join every file into one
        joint_genome = pd.DataFrame()
        filt_genome = filter(lambda x: g in x, args.cnit_list)
        for cnit_genome in filt_genome:
            result = pd.read_csv(cnit_genome, sep="\t")
            result = result.apply(create_columns, axis=1)
            joint_genome = pd.concat(
                [joint_genome, result])
        joint_genome.position = joint_genome.position.astype(
            int)
        joint_genome.sort_values("position", inplace=True)
        joint_genome = joint_genome.rename(columns={"start": "PredictionStart",
                                                    "end": "PredictionEnd", "index": "winner"})
        joint_genome.to_csv("%s_evaluated_cnit.tsv" % g, sep="\t",
                            index=False)

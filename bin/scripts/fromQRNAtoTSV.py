from classes.qrnaParser import QRNAparser

import re
import pandas as pd
import numpy as np
import argparse


def correct_gaps(df, seqs):
    merged = df.merge(seqs, on=["position", "strand"])

    def correct_helper(row):
        prediction_start, prediction_end = row["PredictionStart"], row["PredictionEnd"]
        q = row["aligned_query"]
        correct_start, correct_end = q.count(
            "-", 0, prediction_start), q.count("-", 0, prediction_end)
        strand = row["strand"]
        if row["strand"] == "-":
            row["PredictionStartCorrected"] = row["GenomeEnd"] - \
                prediction_end + correct_end
            row["PredictionEndCorrected"] = row["GenomeEnd"] - \
                prediction_start + correct_start
        else:
            row["PredictionStartCorrected"] = row["GenomeStart"] + \
                prediction_start - correct_start
            row["PredictionEndCorrected"] = row["GenomeStart"] + \
                prediction_end - correct_end
        row = row.drop(["aligned_query", "start", "end"])
        new_order = [17, 0, 2, 18, 19, 3, 4, 5,
                     6, 20, 21, 7, 1] + list(range(8, 17))
        return row[new_order]
    merged = merged.apply(correct_helper, axis=1)
    return merged


def process_blast_results(df):

    def create_columns(row):
        row_id_all = row[0]
        row["transcript_id"] = re.findall(
            "\|((?:orphan_|antisense_)\d+)\|", row_id_all)[0]
        row["position"] = int(re.findall("\|Start:(\d+)\|", row_id_all)[0])
        row["strand"] = re.findall("\|Strand:(\+|\-)", row_id_all)[0]
        row["start_query"] = row[1]
        row["end_query"] = row[2]
        row["aligned_query"] = row[3]
        return row[["transcript_id", "position", "strand", "start_query", "end_query", "aligned_query"]]

    df = df.apply(create_columns, axis=1)
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--qrna_list", nargs="+",)
    parser.add_argument("--filteredBlast_list", nargs="+",)
    args = parser.parse_args()
    genomes = np.unique(
        [re.split("/", re.split("_antisense|_orphan", x)[0])[-1]for x in args.qrna_list])
    for g in genomes:
        # Per genome, return one file for nocornac
        qrna_joint_genome = pd.DataFrame()
        filt_genome = filter(lambda x: g in x, args.qrna_list)
        filt_blast_results = filter(lambda x: g in x, args.filteredBlast_list)
        read_blast_results = pd.concat(
            [pd.read_csv(name, sep="\t", header=None) for name in filt_blast_results])

        processed_blast_results = process_blast_results(
            read_blast_results).sort_values("position")
        for qrna_genome in filt_genome:
            qrna_result = QRNAparser(qrna_genome).to_dataframe()
            qrna_joint_genome = pd.concat(
                [qrna_joint_genome, qrna_result])
        qrna_joint_genome.position = qrna_joint_genome.position.astype(
            int)
        qrna_joint_genome.sort_values("position", inplace=True)
        qrna_joint_genome = correct_gaps(
            qrna_joint_genome, processed_blast_results)
        qrna_joint_genome.to_csv("%s_evaluated_qrna.tsv" % g, sep="\t",
                                 index=False)

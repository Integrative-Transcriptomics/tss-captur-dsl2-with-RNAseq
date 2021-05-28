import re
import pandas as pd
import numpy as np
import argparse


def table_to_gff(df, genome_name):
    gff = df.reset_index().drop("score", axis=1).rename(
        columns={"z_score": "score", "winner": "feature", "transcript_id": "attributes"})
    gff["seqname"] = genome_name
    gff["frame"] = "."
    gff["TranscriptStart"] = gff.position.where(
        gff.strand == "+", gff.PredictionStart)
    gff["TranscriptEnd"] = gff.position.where(
        gff.strand == "-", gff.PredictionEnd)
    gff_utr = gff.copy()
    gff_utr["feature"] = "5' UTR"
    gff_utr["TranscriptStart"] = gff_utr.position.where(
        gff_utr.strand == "+", gff_utr.PredictionEnd)
    gff_utr["TranscriptEnd"] = gff_utr.position.where(
        gff_utr.strand == "-", gff_utr.PredictionStart)
    gff_concat = pd.concat([gff, gff_utr], ignore_index=True)[["seqname", "source", "feature", "TranscriptStart",
                                                               "TranscriptEnd", "score", "strand", "frame", "attributes"]]
    return gff_concat


def table_to_crd(gff, t=0):
    crd = gff.reset_index()
    crd_utr = crd[crd["feature"] == "5' UTR"]

    crd = crd[crd["feature"] != "5' UTR"]
    crd_utr["end_utr"] = crd_utr.TranscriptEnd.where(
        crd_utr.strand == "+", crd_utr.TranscriptStart)
    crd_utr = crd_utr[["attributes", "end_utr"]]
    crd = crd[["attributes", "TranscriptStart",
               "TranscriptEnd", "strand", "feature"]]
    crd_out = crd.merge(crd_utr,
                        on="attributes")
    return crd_out


def normalize_score(df):

    agg_df = df.groupby("winner").agg(mean_score=(
        "score", "mean"), std_score=("score", "std")).to_dict()
    df = df.apply(lambda row: normalize_helper(row, agg_df), axis=1)
    return df


def normalize_helper(row, dict):
    row_winner = row["winner"]
    col_mean = abs(dict["mean_score"][row_winner])
    col_std = dict["std_score"][row_winner]
    row["z_score"] = (abs(row.score) - col_mean) / col_std
    return row


def classification_helper(row):
    output_row = row[["transcript_id", "position", "strand", "alignedTo"]]
    if not pd.isnull(row.winner_qrna):
        same_output = row["winner_qrna"] == row["winner_cnit"]
        output_row["same_classification?"] = same_output
        is_cnit_better = row["z_score_cnit"] > row["z_score_qrna"]

        if row["winner_qrna"] == "COD" and not same_output:
            take_class = "qrna"
        elif row["winner_qrna"] == "OTH" or (row["winner_cnit"] == "COD" and not same_output) or is_cnit_better:
            output_row["alignedTo"] = ""
            take_class = "cnit"
            output_row["same_classification?"] = "OTH Found"
        else:
            take_class = "qrna"
    else:
        take_class = "cnit"
        output_row["same_classification?"] = "No QRNA hit"
    output_row["source"] = take_class.upper()
    similar_output = ["score", "z_score", "winner", "GenomeStart", "PredictionStartCorrected",
                      "GenomeEnd", "PredictionEndCorrected", ]
    for label in similar_output:
        output_row["%s" % label] = row["%s_%s" % (label, take_class)]

    return output_row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnit", nargs="+")
    parser.add_argument("--qrna", nargs="+")
    args = parser.parse_args()
    args.cnit.sort()
    args.qrna.sort()
    compare_list = zip(args.cnit, args.qrna)
    for cnit, qrna in compare_list:
        cnit = re.sub(r",|\[|\]", "", cnit)
        qrna = re.sub(r",|\[|\]", "", qrna)
        genome_name = re.split(
            "/", re.split("_evaluated_cnit", cnit)[0])[-1]
        cnit_df = pd.read_csv(cnit, sep="\t", usecols=[
                              0, 1, 2, 4, 5, 9, 10, 11, 12])
        cnit_df = normalize_score(cnit_df)
        qrna_df = pd.read_csv(qrna, sep="\t", usecols=[
                              0, 1, 2, 5, 6, 9, 10, 11, 12, 19, 20, 21])
        qrna_df["score"] = np.where(
            qrna_df.winner == "RNA", qrna_df.sig_rna, (np.where(
                qrna_df.winner == "COD", qrna_df.sig_cod, qrna_df.sig_oth)))
        qrna_df = normalize_score(qrna_df)
        merged_df = cnit_df.merge(
            qrna_df,  how="left", on=["transcript_id", "position", "strand"], suffixes=["_cnit", "_qrna"])
        summarize_df = merged_df.apply(classification_helper, axis=1).rename(columns={
            "PredictionStartCorrected": "PredictionStart", "PredictionEndCorrected": "PredictionEnd"})

        summarize_df["score"] = summarize_df.score.round(3)
        summarize_df["z_score"] = summarize_df.z_score.round(3)
        tmp = summarize_df.alignedTo.str.extract(
            r".*\|.*\|.*\|(?P<alignedTo>.*)\|.*\|(?P<alignedCoords>.*)", expand=True)
        summarize_df["alignedTo"] = tmp["alignedTo"] + \
            " (" + tmp["alignedCoords"] + ")"
        summarize_df.to_csv("%s_final_classification.tsv" %
                            genome_name, sep="\t", index=False)
        gff = table_to_gff(summarize_df, genome_name)
        gff.to_csv("%s.gff" % genome_name, sep="\t",
                   header=False, index=False)
        crd = table_to_crd(gff)
        crd.to_csv("%s.crd" % genome_name, sep="\t",
                   header=False, index=False)

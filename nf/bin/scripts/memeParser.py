"""Parsing of MEME's output

This script returns a tsv-file with the results computed by MEME on 
the promoter regions of the TSS sites.

This tool expects the XML-output of MEME and the Genome's name as an input. 

BioPython is required for this script. 
"""


from Bio.motifs import meme
import argparse
import re
import csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--meme")
    parser.add_argument("--genome")
    args = parser.parse_args()
    with open(args.meme) as f:
        record = meme.read(f)
    with open('%s_motifs_summary.tsv' % args.genome, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ["Motif ID", "Motif Name", "Evalue", "Length", "Occurrences"])
    with open('%s_seqs_meme_results.tsv' % args.genome, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(["transcript_id", "position", "strand",
                             "Motif ID", "Motif Name", "Motif Start", "Motif Length", "p-value"])
    for index, motif in enumerate(record):
        with open('%s_motifs_summary.tsv' % args.genome, 'a') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            tsv_writer.writerow(
                ["motif_%i" % index, motif.name, motif.evalue, motif.length, motif.num_occurrences])
        for ins in motif.instances:
            sequence_id = re.findall(
                ".*((?:orphan_|antisense_|internal_).*)\|Start:(.*)\|Strand:((?:\+|-))", ins.sequence_name)[0]
            with open('%s_seqs_meme_results.tsv' % args.genome, 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow(
                    list(sequence_id) + ["motif_%i" % index, ins.motif_name, ins.start, ins.length, ins.pvalue])

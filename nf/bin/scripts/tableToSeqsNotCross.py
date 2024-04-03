from classes.tableParserNotCrossOtherTSS import GenomeWrapper

import pandas as pd
import argparse
import os


def parse_table(separator, args):
    table_path = args.table_path
    has_diff_conditions = args.conditions
    with open(table_path, newline='') as csvfile:
        parsedtable = pd.read_csv(csvfile, sep=separator)
    is_detected_enriched = parsedtable["enriched"] == 1
    parsedtable = parsedtable[is_detected_enriched]
    wrapper_list = []
    if has_diff_conditions:
        wrapper = GenomeWrapper(args.allow_tss_intersection,
                                args.allow_annotated_tss_intersection, args.allow_gene_intersection)
        wrapper.process_table(parsedtable, has_diff_conditions)
        wrapper_list.append(wrapper)
    else:
        diffGenomes = pd.unique(parsedtable["Genome"])
        for genome_name in diffGenomes:
            wrapper = GenomeWrapper(args.allow_tss_intersection,
                                    args.allow_annotated_tss_intersection, args.allow_gene_intersection)
            wrapper.add_name(genome_name)
            wrapper.process_table(parsedtable)
            wrapper_list.append(wrapper)

    return wrapper_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("table_path")
    parser.add_argument("genome_path")
    parser.add_argument("gff_path")
    parser.add_argument("--conditions", action="store_true",
                        help="Wether the data comes from a comparison of different conditions")
    parser.add_argument("--allow-gene-intersection", "-G", action="store_true",
                        help="Allow to go over annotated genes for the extraction of transcripts")
    parser.add_argument("--allow-annotated-tss-intersection", "-AT", action="store_true",
                        help="Allow to go over TSSs that correspond to annotated genes for the extraction of transcripts")
    parser.add_argument("--allow-tss-intersection", "-T", action="store_true",
                        help="Allow to go over other TSSs that are being analyzed")

    args = parser.parse_args()

    # print(args)
    # print(args.conditions)
    # print(args.genome_path)
    # print(args.gff_path)

    # if args.conditions and (not os.path.isfile(args.genome_path) or not os.path.isfile(args.gff_path)):
    #     parser.error(
    #         "If conditions are been tested, you need to upload only one genome or gff file")
    # # If strains, then check that the paths are directories
    # print(args.conditions)
    # print(os.path.isdir(str(args.genome_path)))
    # print(os.path.isdir(str(args.gff_path)))
    # if (not args.conditions) and ((not os.path.isdir(str(args.genome_path))) or (not os.path.isdir(str(args.gff_path)))):
    #     parser.error(
    #         "If strains are been tested, you need to upload only a directory for genomes or gff files")

    genome_wrapper_list = parse_table("\t", args)

    records = []
    genomes_names = []
    for wrapper in genome_wrapper_list:
        wrapper.add_genome(args.genome_path)
        wrapper.add_gff(args.gff_path)
        genomes_names.append(wrapper.genome_name)
        records = records + wrapper.create_queries()
        with open("tss_avoided_in_%s.tsv" % wrapper.genome_name, "w") as f:
            list_all = wrapper.positions_as_vectors
            list_created = wrapper.created_queries
            not_created = set(list_all) - set(list_created)
            for tss, strand in not_created:
                f.write("%i\t%s\n" % (tss, strand))
    with open('genomes_text.txt', 'w') as f:
        for item in genomes_names:
            f.write("%s\n" % item)

    with open("tss_analyzed.tsv", "w") as file_writer:
        for id_record, length_query in records:
            file_writer.write("%s\t%i\n" % (id_record, length_query))

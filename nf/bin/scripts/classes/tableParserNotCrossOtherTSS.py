import pandas as pd
import math
import os
import re
import random
from pathlib import Path
from Bio import SeqIO, SeqRecord
pd.options.mode.chained_assignment = None  # default='warn'


class GenomeWrapper(object):
    def __init__(self, tss, annotated_tss, genes):
        self.allow_gene_overlap = genes
        self.allow_overlap_other_tss = annotated_tss
        self.allow_overlap_other_same_category = tss
        self.created_queries = []

    def add_name(self, name):
        self.genome_name = name

    def add_genome(self, genome_path):
        is_conditional = self.type == "conditions"
        self.genome_path = self.get_path(genome_path, is_conditional)
        with open(self.genome_path, "r") as handle:
            record = list(SeqIO.parse(handle, "fasta"))[0]
            self.genome_id = record.id
            if is_conditional:
                self.genome_name = record.name.split("|")[-2]

    def add_gff(self, gff_path):
        is_conditional = self.type == "conditions"
        self.gff_path = self.get_path(gff_path, is_conditional, True)
        with open(self.gff_path, newline='') as csvfile:
            parsedGFF = pd.read_csv(
                csvfile, sep="\t", comment="#", header=None)
            is_gene = parsedGFF[2] == "gene"
            parsedGFF = parsedGFF[is_gene]
            parsedGFF["gene_length"] = parsedGFF[4] - parsedGFF[3]

        self.mean_gene_length = math.ceil(parsedGFF["gene_length"].mean())

    def get_max_length_per_tss(self, table):
        '''
        Sets the max length for each TSS to avoid the overlap with other annotated TSSs of any kind
        '''
        data_indexed = index_positions(table)

        table.loc[:, "geneLengthMax"] = table.apply(lambda x: self.get_gene_length_helper(x, data_indexed),
                                                    axis=1)
        return table

    def get_gene_length_helper(self, row, data_indexed):
        data_entry = data_indexed[row.Strand][row.Pos]
        shift1 = data_entry["Shift1"]
        return shift1 if ((shift1 >= 30) or (pd.isnull(shift1))) else data_entry["Shift2"]

    def process_table(self, table, has_dif_conditions=False):
        if has_dif_conditions:
            self.type = "conditions"
            table = table.drop_duplicates(["Pos", "Strand"])
        else:
            self.type = "strains"
            table = table[table["Genome"] == self.genome_name]
        important_columns = ["SuperPos", "SuperStrand", "Pos", "Strand",
                             "Genome", "Locus_tag", "Sequence -50 nt upstream + TSS (51nt)", "geneLengthMax"]
        table = self.get_max_length_per_tss(table)
        # Separate all orphans
        is_orphan = table["Locus_tag"] == "orphan"
        if(is_orphan.any()):
            all_orphans = table[is_orphan].reset_index()

            all_orphans = all_orphans[important_columns]
            all_orphans["type"] = "orphan"
            self.orphans = separate_table(all_orphans)
        else:
            self.orphans = [pd.DataFrame(), pd.DataFrame()]

        # Separate all antisense
        is_antisense = table["Antisense"] == 1
        if(is_antisense.any()):

            not_table = table[~ (is_antisense | is_orphan)]
            index_not_table = (
                not_table[["Pos", "Strand"]].to_records(index=False))
            all_antisense = table[is_antisense]

            # Filtering out repetitive TSS
            # A TSS can be associated with many genes, hence an AS might have been already labelled as P/S/I
            # This need to be removed

            filter_repeated = [not x in index_not_table for x in all_antisense[[
                "Pos", "Strand"]].to_records(index=False)]

            all_antisense = all_antisense[filter_repeated]

            # An AS TSS can also be categorized wrt. many genes, hence the duplicates need to be removed
            all_antisense = all_antisense.drop_duplicates(
                ["Pos", "Strand"]).reset_index()
            # print(len(all_antisense))

            all_antisense = all_antisense[important_columns]
            all_antisense["type"] = "antisense"

            self.antisense = separate_table(all_antisense)
        else:
            self.antisense = [pd.DataFrame(), pd.DataFrame()]
        #self.list_positions()

        #seperate all internals
        is_internal = table["Internal"] == 1
        if(is_internal.any()):
            not_table = table[~ (is_antisense | is_orphan | is_internal)]

            index_not_table =(
                not_table[["Pos", "Strand"]].to_records(index=False))
            
            all_internals = table[is_internal]

                        # Filtering out repetitive TSS
            # A TSS can be associated with many genes, hence an AS might have been already labelled as P/S/I
            # This need to be removed

            filter_repeated = [not x in index_not_table for x in all_internals[[
                "Pos", "Strand"]].to_records(index=False)]
                
            all_internals = all_internals[filter_repeated]

            # An AS TSS can also be categorized wrt. many genes, hence the duplicates need to be removed
            all_internals = all_internals.drop_duplicates(
                ["Pos", "Strand"]).reset_index()
            
            all_internals = all_internals[important_columns]
            all_internals["type"] = "internal"

            self.internals = separate_table(all_antisense)
        else:
            self.internals = [pd.DataFrame(), pd.DataFrame()]
        self.list_positions()


    def list_positions(self):
        """
        Creates a list of all TSS-positions that were labelled as sole AS or orphans
        """
        positions_as_vectors = []
        for x in self.antisense + self.orphans + self.internals:
            if len(x) > 0:
                vectors = [(v["Pos"], v["Strand"])
                           for v in x.to_dict("records")]
                positions_as_vectors = positions_as_vectors + vectors
        self.positions_as_vectors = positions_as_vectors

    def create_queries(self):

        orphan_records, orphan_promoters = self.helper_create_queries(
            self.orphans, "orphan")

        antisense_records, antisense_promoters = self.helper_create_queries(
            self.antisense, "antisense")

        internal_records, internal_promoters = self.helper_create_queries(self.internals, "internal")

        promoter_regions = antisense_promoters + orphan_promoters + internal_promoters

        if (len(promoter_regions) > 0):  # Create file only if there is a sequence to save in
            if (len(promoter_regions) > 1000):
                random.seed(123)
                promoter_regions = random.sample(promoter_regions, 1000)
            with open("%s_%s.fasta" % (self.genome_name, "promoter_regions"), "w") as output_handle:
                SeqIO.write(promoter_regions, output_handle, "fasta")
        records = orphan_records + antisense_records + internal_records

        ids_of_records = [[item.id, len(item.seq)]
                          for sublist in records for item in sublist]

        filenames = ["orphan_plus", "orphan_minus",
                     "antisense_plus", "antisense_minus", "internal_plus", "internal_minus"]
        for name, seqs in zip(filenames, records):
            if (len(seqs) > 0):  # Create file only if there is a sequence to save in
                with open("%s_%s_queries.fasta" % (self.genome_name, name), "w") as output_handle:
                    SeqIO.write(seqs, output_handle, "fasta")
        return ids_of_records

    def get_gene_positions(self):
        with open(self.gff_path, newline='') as csvfile:
            parsedGFF = pd.read_csv(
                csvfile, sep="\t", comment="#", header=None)
        is_gene = parsedGFF[2] == "gene"
        parsedGFF = parsedGFF[is_gene]
        genes_in_plus = list(parsedGFF[parsedGFF[6] == "+"][3])
        genes_in_minus = list(parsedGFF[parsedGFF[6] == "-"][4])

        return genes_in_plus, genes_in_minus

    def get_genome(self):
        with open(self.genome_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                genome = record
        return genome

    def extract_queries(self, table, found_genes, tss_type, strand, mean_length):
        genome = self.get_genome()
        genome_length = len(genome.seq)
        list_records = []
        list_promoters = []
        index_found_genes = 0
        # To add the column to distance with the previous or next TSS
        table = shift_table(table)
        for index_tss, row in table.iterrows():

            tss_start = row["Pos"]
            previous_tss = row["prevPos"]
            next_tss = row["nextPos"]
            is_searching = True
            while is_searching:
                if index_found_genes >= len(found_genes)-1:
                    index_found_genes = len(found_genes)-1
                    is_searching = False
                elif tss_start < found_genes[index_found_genes]:
                    if strand == "-":
                        index_found_genes -= 1  # Go one backwards to find the closest gene
                        index_found_genes = max(index_found_genes, 0)
                    is_searching = False
                else:
                    index_found_genes += 1

            dist_tss_gene = abs(found_genes[index_found_genes] - tss_start)
            # Compute the distance to the next TSS
            if strand == "+":
                # for the plus strand we extend in the forward direction
                dist_tss_next_tss = next_tss - \
                    tss_start if not math.isnan(
                        next_tss) else mean_length  # if last TSS then take the mean length
            else:
                dist_tss_next_tss = tss_start - \
                    previous_tss if not math.isnan(  # Here we extend in the other direction
                        previous_tss) else mean_length
            dist_tss_next_tss = int(dist_tss_next_tss)
            dist_all_tss = int(mean_length if pd.isnull(
                row["geneLengthMax"]) else row["geneLengthMax"])
            distances_to_account = [mean_length]
            if not self.allow_gene_overlap:
                distances_to_account.append(dist_tss_gene)
            if not self.allow_overlap_other_tss:
                distances_to_account.append(dist_all_tss)
            if not self.allow_overlap_other_same_category:
                distances_to_account.append(dist_tss_next_tss)

            # Limit length to 1000 to avoid problems with QRNA
            length = min(int(min(distances_to_account)), 1000)

            if strand == "+":
                start_query = tss_start - 1  # Correct for python starting at 0
                end_query = tss_start - 1 + length
                promoter_start = max(start_query - 50, 0)
                promoter_end = start_query + 1  # Correct for python not including final query

            else:
                start_query = max(tss_start - length, 0)
                # Not correcting since the position it takes at the end is tss_start -1
                end_query = tss_start
                promoter_start = tss_start - 1  # Correct for position in python
                promoter_end = min(promoter_start + 51, genome_length)
            # For debugging reasons, avoid overtaking the lenght of the genome
            end_query = min(end_query, genome_length)

            query = genome.seq[start_query:end_query]
            promoter_query = genome.seq[promoter_start:promoter_end]

            if strand == "-":
                query = query.reverse_complement()
                promoter_query = promoter_query.reverse_complement()
            if len(query) > 0:
                record = self.produce_record(
                    query, tss_type, index_tss, tss_start, strand)

                list_records.append(record)
            if len(promoter_query) > 0:
                promoter_record = self.produce_record(
                    promoter_query, tss_type, index_tss, tss_start, strand, True)
                list_promoters.append(promoter_record)
            self.created_queries.append((tss_start, strand))
        return list_records, list_promoters

    def produce_record(self, query, tss_type, index, tss_pos, strand, is_promoter=False):
        tmp = ""
        if not is_promoter:
            tmp = "%s|" % self.genome_id
        id = "%s%s_%i|Start:%i|Strand:%s" % (tmp,
                                             tss_type, index, tss_pos, strand)
        id = id.replace("||", "|")
        name = "%s|%s_%i"
        return SeqRecord.SeqRecord(query, id=id, name=name, description="")

    def helper_create_queries(self, tables, tss_type,):
        plus_table, minus_table = tables
        mean_length = self.mean_gene_length
        genes_in_plus, genes_in_minus = self.get_gene_positions()
        plus_queries, plus_promoters = self.wrapper_table(
            plus_table, "+", mean_length, genes_in_plus, tss_type)
        minus_queries, minus_promoters = self.wrapper_table(
            minus_table, "-", mean_length, genes_in_minus, tss_type)
        return [plus_queries, minus_queries], plus_promoters + minus_promoters

    def get_path(self, path, is_conditional, is_gff=False):
        if is_conditional:
            return path
        else:
            filelist = os.listdir(path)
            fasta_ext = "(fa|fna|fasta|faa|frn|ffn)"
            gff_ext = "(gff|gff3)"
            print("isgff: ", is_gff, path, filelist)            
            regex = re.compile(
                r".*\.(%s)$" % (gff_ext if is_gff else fasta_ext))
            match = list(filter(regex.match, filelist))[0]
            return os.path.join(path, match)

    def wrapper_table(self, table, strand, mean_length, gene_list, tss_type):
        result_queries = []
        result_promoters = []
        if len(table) > 0:
            table_mod = compare_distances_tss(table, strand)
            result_queries, result_promoters = self.extract_queries(
                table_mod, gene_list, tss_type, strand, mean_length)
        return result_queries, result_promoters


def index_table(table, filter_mask, strand):
    subtable = table[filter_mask].drop_duplicates("Pos")
    index_direction = -1 if strand == '+' else 1
    def shifting(x): return abs(subtable["Pos"].shift(
        index_direction * x) - subtable["Pos"])
    subtable["Shift1"] = shifting(1)
    subtable["Shift2"] = shifting(2)
    indexed = subtable.set_index("Pos").to_dict('index')
    return indexed


def index_positions(table):
    sorted_table = table[["Pos", "Strand"]].sort_values("Pos")
    plus_entries = sorted_table["Strand"] == "+"
    indexed_plus = index_table(sorted_table, plus_entries, "+")
    indexed_minus = index_table(sorted_table, ~plus_entries, "-")
    data_indexed = {"+": indexed_plus, "-": indexed_minus}
    return data_indexed


def separate_table(table):
    table_plus = table["Strand"] == "+"
    return [table[table_plus].sort_values("Pos"), table[~table_plus].sort_values("Pos")]


def compare_distances_tss(table, strand):
    '''
    Filters the table in such a way that every uncl. TSS has a distance of 30 nt to its unclassified previous TSS. 
    Otherwise it is removed to leave a larger path to the previous TSS. 
    '''
    mod_table = shift_table(table)
    if strand == "+":
        mod_table = mod_table[(
            (mod_table.Pos - mod_table.prevPos) >= 30) | mod_table.prevPos.isnull()]

    else:
        mod_table = mod_table[((mod_table.nextPos -
                                mod_table.Pos) >= 30) | mod_table.nextPos.isnull()]
    return mod_table.drop(["prevPos", "nextPos"], axis=1)


def shift_table(table):
    mod_table = table
    mod_table["prevPos"] = table["Pos"].shift()
    mod_table["nextPos"] = table["Pos"].shift(-1)
    return mod_table

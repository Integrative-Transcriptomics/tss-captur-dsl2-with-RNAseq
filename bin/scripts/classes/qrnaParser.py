import re
import pandas as pd


class QRNAparser(object):
    def __init__(self, path):
        classed_sqs = []

        with open(path, 'r', encoding='utf-8') as qrna:
            for l in qrna:
                logoddpattern = re.compile("logodd\w*\s=\s*(-?[0-9]*.[0-9]*)")
                sigpattern = re.compile("sigm\w*\s=\s*(-?[0-9]*.[0-9]*)")
                normalpattern = re.compile(
                    "\s[{OTH}|{COD}|{RNA}]\w*\s=\s*(-?[0-9]*.[0-9]*)")
                winner_coor = "%s ends \*\(.\) =\s*\((\w*)\.\..*\.\.(\w*)\)*\s"
                if l[0] == ">":
                    removed_ws = l[1:].split()[0]
                    if "Start:" in l:
                        splitted_info = removed_ws.split("|")

                        genome = "|".join(splitted_info[0:4])
                        [transcript_ID, tss_pos, strand,
                            blasted_positions, scores] = splitted_info[-5:]
                        tss_pos = tss_pos.split(":")[1]
                        strand = strand.split(":")[1]
                        [start, end] = blasted_positions.split("-")

                        current_seq = QRNAelement(
                            genome, transcript_ID, tss_pos, strand)
                        current_seq.add_positions(start, end)
                    else:
                        current_seq.add_alignedTo(removed_ws)
                elif re.match("(winner)\s=", l):
                    current_seq.add_winner(l.split("=")[1].strip())
                elif len(re.findall(winner_coor % "RNA", l)) > 0:
                    found_coor = re.findall(winner_coor % "RNA", l)
                    current_seq.add_pred_rna(*found_coor)
                elif len(re.findall(winner_coor % "COD", l)) > 0:
                    found_coor = re.findall(winner_coor % "COD", l)
                    current_seq.add_pred_cod(*found_coor)
                elif len(normalpattern.findall(l)) > 0:
                    found = normalpattern.findall(l)

                    current_seq.add_normal(*found)
                elif len(logoddpattern.findall(l)) > 0:
                    found = logoddpattern.findall(l)
                    current_seq.add_log(*found)
                elif len(sigpattern.findall(l)) > 0:
                    found = sigpattern.findall(l)
                    current_seq.add_sig(*found)
                    classed_sqs.append(current_seq)
        self.classed_sqs = classed_sqs

    def to_dataframe(self):
        list_export = [x.asdict() for x in self.classed_sqs]
        return (pd.DataFrame(list_export).set_index(["genome", "ID"]))


class QRNAelement(object):

    def __init__(self, genome, ID, position, strand):
        self.genome = genome
        self.ID = ID
        self.position = position
        self.strand = strand

    def add_positions(self, start, end):
        self.start = start
        self.end = end

    def add_scores(self, evalue, bitscore):
        self.evalue = evalue
        self.bitscore = bitscore

    def add_pred_rna(self, coor):
        self.rna_start = coor[0]
        self.rna_end = coor[1]

    def add_pred_cod(self, coor):
        self.cod_start = coor[0]
        self.cod_end = coor[1]

    def add_winner(self, winner):
        self.winner = winner

    def add_alignedTo(self, alignedTo):
        self.alignedTo = alignedTo

    def add_sig(self, oth, cod, rna):
        self.sig_oth = oth
        self.sig_cod = cod
        self.sig_rna = rna

    def add_normal(self, oth, cod, rna):
        self.normal_oth = oth
        self.normal_cod = cod
        self.normal_rna = rna

    def add_log(self, oth, cod, rna):
        self.log_oth = oth
        self.log_cod = cod
        self.log_rna = rna

    def asdict(self):
        prediction_coord = prediction_wrapper(self, self.winner)
        gStart = translate_coord(
            self.position, self.start, self.end, self.strand)
        gEnd = translate_coord(self.position,  self.end,
                               self.start, self.strand)
        return {
            "genome": self.genome,
            "ID": self.ID,
            "position": self.position,
            "alignedTo": self.alignedTo,
            # "evalue": self.evalue,
            # "bitscore": self.bitscore,
            "strand": self.strand,
            "start": self.start,
            # -1 for correcting since blast takes 0 as 1
            "GenomeStart": gStart,
            "GenomeEnd": gEnd,
            "PredictionStart": int(prediction_coord[0]),
            "PredictionEnd": int(prediction_coord[1]),
            "end": self.end,
            "winner": self.winner,
            "normal_oth": self.normal_oth,
            "normal_cod": self.normal_cod,
            "normal_rna": self.normal_rna,
            "log_oth": self.log_oth,
            "log_cod": self.log_cod,
            "log_rna": self.log_rna,
            "sig_cod": self.sig_cod,
            "sig_rna": self.sig_rna,
            "sig_oth": self.sig_oth}


def prediction_wrapper(obj, winner):
    if (winner == "RNA"):
        return [obj.rna_start, obj.rna_end]
    elif (winner == "COD"):
        return [obj.cod_start, obj.cod_end]
    else:
        return [obj.start, obj.end]


def translate_coord(pos, plus_corr, minus_corr, strand):
    return int(pos) + (int(plus_corr)-1 if strand == "+" else -1 * (int(minus_corr))+1)

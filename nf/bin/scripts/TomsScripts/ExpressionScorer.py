import csv
import re
import CalcbackgroundNoise
import DerivScoring
import AvgScoring
import pyBigWig as bigwig
import gff3_parser
import pandas as pd
import numpy as np

def FindFirstUpstreamTSS(position, minimum_gene_length, tssList, strand):
    lastPos = -1
    if strand == '+':
        for tss in tssList:
            if(tss >= position):
                return lastPos
            
            if(tss <= position - minimum_gene_length):
                lastPos = tss
        return lastPos
    else:
        for tss in tssList:
            if(tss > position + minimum_gene_length):
                return tss
    return -1

def score_genome_region(forward_bigwig_path, 
                      reverse_bigwig_path, annotation_path, 
                      gffRhoterm, gffNocornac, master_table_path):
    

    MINIMUM_GENE_LENGTH = 15



    forward_noise, forward_iqr = CalcbackgroundNoise.InverseOfMasterTableNoise(annotation_path, forward_bigwig_path, master_table_path)
    reverse_noise, reverse_iqr = CalcbackgroundNoise.InverseOfMasterTableNoise(annotation_path, reverse_bigwig_path, master_table_path)

    fwbw = bigwig.open(forward_bigwig_path)
    rvbw = bigwig.open(reverse_bigwig_path)

    rhotermdata = gff3_parser.parse_gff3(gffRhoterm, verbose = False, parse_attributes = False)
    nocornacdata = gff3_parser.parse_gff3(gffNocornac, verbose = False, parse_attributes = False)

    all_term_data = pd.concat([rhotermdata, nocornacdata], ignore_index = True)

    del rhotermdata
    del nocornacdata

    forward_chromosome = re.match(r'^[^\.]+', rhotermdata.loc[1, "Seqid"]).group(0)

    for chrom in fwbw.chroms():
        prefix = re.match(r'^[^\.]+', chrom).group(0)
        if(prefix == forward_chromosome):
            forward_chromosome = chrom
            break

    reverse_chromosome = re.match(r'^[^\.]+', rhotermdata.loc[1, "Seqid"]).group(0)

    for chrom in rvbw.chroms():
        prefix = re.match(r'^[^\.]+', chrom).group(0)
        if(prefix == reverse_chromosome):
            reverse_chromosome = chrom
            break
    

    all_terms_intervalls = zip(all_term_data.loc[:, "Start"].astype(int), all_term_data.loc[:, "End"].astype(int), all_term_data.loc[:, "Score"],
                                all_term_data.loc[:, "Type"], all_term_data.loc[: ,"Strand"])
    

    MasterTable = pd.read_csv(master_table_path, sep="\t")

    positive_tss = sorted(MasterTable[MasterTable['SuperStrand'] == '+']['SuperPos'].tolist())
    negative_TSS = sorted(MasterTable[MasterTable["SuperStrand"] == '-']['SuperPos'].tolist())

    del MasterTable
    with open("AllTermScoringPlusInfo.csv", 'w', newline ='') as info, open("AllTermScoring.tsv", 'w', newline = "") as allTermScoring:
        infoWriter = csv.writer(info, delimiter=",")
        termWriter = csv.writer(allTermScoring, delimiter='\t')

        termWriter.writerow(["seqid", "myTSS", "type", "strand", "start", "end", "initalScore", "avgScore", "dropScore", "derivScore"])
        infoWriter.writerow(["seqid", "myTSS", "type", "strand", "start", "end", "initalScore", "avgScore",
                      "derivScore", "bgNoise", "derivScoreWindowStart", "derivScoreWindowEnd", "avgScoreWindowStart", "avgScoreWindowEnd"])

        for start, end, score, type, strand in all_terms_intervalls:
            if(type != 'terminator' and type != 'RhoTerminator'):
                continue
            
            #Get TSS
            if(strand == '+'):
                bw = fwbw
                chromosome = forward_chromosome
                noiseLvL = forward_noise
                iqr = forward_iqr
                myTss = FindFirstUpstreamTSS(start, MINIMUM_GENE_LENGTH, positive_tss, '+')
                #No TSS found, no gene here, we assume
                if(myTss <= -1):
                    continue
                estGeneExprMed = np.quantile(bw.values(chromosome, myTss, end), 0.5)
            else:
                bw = rvbw
                chromosome = reverse_chromosome
                noiseLvL = reverse_noise
                iqr = reverse_iqr
                myTss = FindFirstUpstreamTSS(end, MINIMUM_GENE_LENGTH, negative_TSS, '-')
                #print(f"Start: {start} End: {end} tss: {myTss}")
                #No TSS found, no gene here, we assume
                if(myTss <= -1):
                    continue
                estGeneExprMed = np.quantile(bw.values(chromosome, start, myTss), 0.5)

            #If gene is not expressed, don't score at all
            if(estGeneExprMed <= noiseLvL):
                termWriter.writerow([chromosome, myTss, type, strand, start, end, score, "NA", "NA", "NA"])

            #Get avg score
            #get deriv score
            #get drop score
            #write to info and non info tsv


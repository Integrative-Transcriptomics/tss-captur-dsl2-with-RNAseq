"""Calculator of Background noise along bigWig file

This script calculates the background noise of a biugwig file by averaging the expression level across 
all regions present in the provided annotation (gff) file and NOT present in the MasterTable. Also 
return the IQR of the expression levels in these regions.

This script requires `gff3_parser` `numpy` `pandas` to be installed. 
"""

import glob
import gff3_parser
import numpy as np
import pandas as pd

def InverseOfMasterTableNoise(annot_path, bw, masterTable_path):
    chunksize = 10**6
    print(masterTable_path)
    master_table_chunks = pd.read_csv(masterTable_path, delimiter='\t', chunksize=chunksize)

    allMasterTableLocusTags = {}
    for chunk in master_table_chunks:
        locus_tags = list(map(str.strip, chunk['Locus_tag'].astype(str).tolist()))
        super_strands = list(map(str.strip, chunk['SuperStrand'].astype(str).tolist()))
        for locus_tag, strand in zip(locus_tags, super_strands):
            allMasterTableLocusTags[locus_tag] = strand

    files = glob.glob(f"{annot_path}/*.gff")
    files.extend(list(glob.glob(f"{annot_path}/*.gff3")))
    annot = gff3_parser.parse_gff3(files[0], verbose = False, parse_attributes = True) 
    chromosome = annot.loc[1, "Seqid"]

    accepted_types = ['gene']
    unexpressedCDS = []

    allnames = []
    for start, end, feature, strand, name in zip(
            annot.loc[:, "Start"].astype(int),
            annot.loc[:, "End"].astype(int),
            annot.loc[:, "Type"],
            annot.loc[:, "Strand"],
            annot.loc[:, "Name"].astype(str)):
        if (feature in accepted_types):
            allnames.append(name)
            if allMasterTableLocusTags.get(name, None) == strand:
                unexpressedCDS.append((start, end))

    fused_unexpressed_CDS = []

    fused_unexpressed_CDS = unexpressedCDS

    all_expr_vals = []

    iqr = np.quantile(bw.values(chromosome, 1, bw.chroms(chromosome)), 0.75) - np.quantile(bw.values(chromosome, 1, bw.chroms(chromosome)), 0.25)
    
    for start, end in fused_unexpressed_CDS:
        vals = bw.values(chromosome, start, end)
        all_expr_vals.extend(vals)

    all_expr_vals = np.array(all_expr_vals)
    
    return np.quantile(all_expr_vals, 0.5), iqr
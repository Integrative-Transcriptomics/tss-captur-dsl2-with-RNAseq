import glob
import itertools
import gff3_parser
import pyBigWig as bigwig
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def CalcBackgroundNoise(annotationPath, bigWigPath):

    bw = bigwig.open(bigWigPath)


    files = glob.glob(f"{annotationPath}/*.gff")
    print(f"annnot path bg: {annotationPath}")
    print(f"anoot globs bg: {files}")
    annot = gff3_parser.parse_gff3(files[0], verbose = False, parse_attributes = False)


    chromosome = annot.loc[1, "Seqid"]


    inverse = GetInverseOfGFF(annotationPath, bw.chroms(chromosome))

    #print(np.quantile(valueCopy, 1))
    #print(nonTranscr)

    quantiles = []
    for intervall in inverse:
        quantiles.append(np.quantile(bw.values(chromosome, intervall[0], intervall[1]), 0.5))
        
    return(np.quantile(quantiles, 0.5))



def GetInverseOfGFF(annotationPath, chromSize):


    files = glob.glob(f"{annotationPath}/*.gff")
    print(f"annnot path inverse: {annotationPath}")
    print(f"anoot globs inverse: {files}")
    annot = gff3_parser.parse_gff3(files[0], verbose = False, parse_attributes = False)  

    acceptedtypes = ["CDS", "gene"]

    #valueCopy = bw.values(chromosome, 1, bw.chroms(chromosome))
    transcr = []
    #print(bw.chroms(chromosome))

    #print(len(valueCopy))
    #might need to filter for certain regions here as well
    for start, end, feature, strand in zip(annot.loc[:, "Start"].astype(int).tolist(), annot.loc[:, "End"].astype(int).tolist(), annot.loc[:, "Type"],annot.loc[:, "Strand"]):
        if(strand != "+" or (not feature in acceptedtypes)):
            continue
        #valueCopy[start:end] = itertools.repeat("out", (end - start))
        transcr.append((start,end))
    #valueCopy = [i for i in valueCopy if i!= "out"]

    start = 0

    # Sort intervals by start position
    transcr.sort(key=lambda x: x[0])

    fused_transcr = []
    current_start, current_end = transcr[0]

    for start, end in transcr[1:]:
        # If there's no gap, extend the current interval
        if start <= current_end:
            current_end = max(current_end, end)
        # If there's a gap, add the current interval and start a new one
        else:
            fused_transcr.append((current_start, current_end))
            current_start, current_end = start, end

    # Add the last interval
    fused_transcr.append((current_start, current_end))

    inverse = []

    ##if first transcribed intervall starts after 1, pre intervall region is not trranscribed
    if(1 < fused_transcr[0][0]):
        inverse.append((1, fused_transcr[0][0]))

    for en in enumerate(fused_transcr):
        if(en[0] < (len(fused_transcr) -1)):
            #add intervall from end of current intervall till start of next one
            inverse.append((en[1][1], fused_transcr[en[0] + 1][0]))

    #if end of last transcribed intervall is before end of chromosome, add last untranscribed intervall till end
    if(chromSize > fused_transcr[-1][1]):
        inverse.append((fused_transcr[-1][1], chromSize))
    
    return inverse

def InverseOfMasterTableNoise(annot_path, bigwig_path, masterTable_path):
    chunksize = 10**6
    print(masterTable_path)
    master_table_chunks = pd.read_csv(masterTable_path, delimiter='\t', chunksize=chunksize)

    allMasterTableLocusTags = {}
    for chunk in master_table_chunks:
        locus_tags = list(map(str.strip, chunk['Locus_tag'].astype(str).tolist()))
        super_strands = list(map(str.strip, chunk['SuperStrand'].astype(str).tolist()))
        for locus_tag, strand in zip(locus_tags, super_strands):
            allMasterTableLocusTags[locus_tag] = strand

    # for locustag in allMasterTableLocusTags:
    #     print(locustag)

    print(len(allMasterTableLocusTags))

    files = glob.glob(f"{annot_path}/*.gff")
    files = files.extend(glob.glob(f"{annot_path}/*.gff3"))
    print(f"annnot path inverse: {annot_path}")
    print(f"anoot globs inverse: {files}")
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
        if (feature in accepted_types): #and strand == '+':
            allnames.append(name)
            #print(name)
            if allMasterTableLocusTags.get(name, None) == strand:
                unexpressedCDS.append((start, end))

    print(len(allnames))
    # print(len(allMasterTableLocusTags.intersection(allnames)))
    #print(len(allMasterTableLocusTags - set(allnames)))

        # if feature not in accepted_types or name in allMasterTableLocusTags:
        #     continue
        # else:
        #     unexpressedCDS.append((start,end))

    current_start, current_end = unexpressedCDS[0]
    fused_unexpressed_CDS = []

    # for start, end in unexpressedCDS[1:]:
    #     if start <= current_end:
    #         current_end = max(current_end, end)
    #     else:
    #         fused_unexpressed_CDS.append((current_start, current_end))
    #         current_start, current_end = start, end

    # Add the last interval
    #fused_unexpressed_CDS.append((current_start, current_end))

    fused_unexpressed_CDS = unexpressedCDS

    print(len(unexpressedCDS))

    all_expr_vals = []

    with bigwig.open(bigwig_path) as bw:
        IQR = np.quantile(bw.values(chromosome, 1, bw.chroms(chromosome)), 0.75) - np.quantile(bw.values(chromosome, 1, bw.chroms(chromosome)), 0.25)
        for start, end in fused_unexpressed_CDS:
            vals = bw.values(chromosome, start, end)
            #print(vals)
            all_expr_vals.extend(vals)

    all_expr_vals = np.array(all_expr_vals)
    uniques, counts = np.unique(all_expr_vals, return_counts=True)
    print(len(all_expr_vals))

    plt.bar(uniques, counts, width=1)
    plt.figure()
    fig, ax = plt.subplots()
    ax.boxplot(all_expr_vals, showfliers=False)
    ax.axhline(y = np.quantile(all_expr_vals, 0.5), color= "r", linewidth = 1)
    plt.show()

    
    return np.quantile(all_expr_vals, 0.5), IQR
/*
    Extracts transcripts from the MasterTable of TSSpredator
    Returns queries for classification, promoter regions, a .txt file with genome names, and a .tsv file with ignored and analyzed positions
*/
process MASTERTOFASTA {
    label params.dataPrepLabel
    publishDir "$outputPath/Queries", pattern: "*.tsv", mode: params.pubDirMode
    publishDir "$outputPath/Queries", pattern: "*.fasta" , mode: params.pubDirMode
    publishDir "$outputPath", pattern: "genomes_text.txt" , mode: params.pubDirMode

    input: 
    val masterTable
    val genomePath
    val gffPath
    val outputPath

    output: 
    path "genomes_text.txt", emit: fileText 
    path "*_queries.fasta", emit: queries
    path "*_promoter_regions.fasta", emit: promoters
    path "tss_analyzed.tsv", emit: summaryTranscripts
    path "*.tsv", emit: analyzedTSS

    script:
    """
        python3 $params.pyScriptTableToQuery $masterTable $genomePath $gffPath
    """
}

/*
    Extracts the corresponding TaxIDs from the Taxonomy tree of NCBI using the NCBI-IDs of the given genomes (as a single call)
*/
process GETGENOMESLCA {
    label params.dataPrepLabel
    tag "$ids"
    
    input: 
    val ids

    output: 
    stdout emit: taxID

    script:
    """
        efetch -db nuccore -id $ids -format docsum | xtract -pattern DocumentSummary -element TaxId 
    """
}

/*
    Extracts the LCA from the species in question using the given TaxIDs
*/
process GETLCAID {
    label params.dataPrepLabel

    input: 
    val ids

    output: 
    path "common_species_id.txt", emit: commonSpeciesID

    script: 
    """
        python3 $params.pyScriptCommonSpecies $ids > common_species_id.txt
    """
}

/*
    Creates a TaxIDlist for BLAST with all corresponding species under the extracted LCA
*/
process GETBLASTID {
    label params.dataPrepLabel

    input: 
    val taxIDList

    output: 
    path "taxidlist.taxid", emit: taxIDListFile

    script:
    """
        get_species_taxids.sh -t $taxIDList > taxidlist.taxid
    """
}

/*
    Runs BLAST against the nt-database using the extracted queries and the restricted TaxIDs
*/
process BLASTFASTA {
    label params.dataPrepLabel

    input:
    each query
    path taxIDList
    env BLASTDB
    val outputPath

    output:
    path "*.blastn", emit: blastFiles

    script: 
    """
        blastn -query $query -db nt_prok -out ${query.baseName - "_queries"}.blastn -outfmt "6 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore pident frames qcovhsp" -task dc-megablast -taxidlist $taxIDList
    """
}

/*
    Evaluates the BLAST-hits and extracts the corresponding hits for the pairwise alignment (for QRNA)
*/
process EVALBLAST {
    label params.dataPrepLabel
    publishDir "$outputPath/BlastFiltered", mode: params.pubDirMode

    input: 
    path blastFiles
    val outputPath

    output:
    path "*.tsv", emit: blastFiltered

    // path "evaluation*" into evaluation_ch
    // TODO: Still add the evaluation table, since there might be transcripts without any hit. 
    script:
    """
        python3 $params.pyEvaluateBlast $blastFiles -t 10
    """
}
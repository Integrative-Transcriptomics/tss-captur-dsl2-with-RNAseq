/*
    Extracts transcripts from the MasterTable of TSSpredator
    Returns queries for classification, promoter regions, a .txt file with genome names, and a .tsv file with ignored and analyzed positions
*/
process MASTERTOFASTA {
    label "DataPreparation"
    publishDir "$output_path/queries", pattern: "*.tsv", mode: 'copy'
    publishDir "$output_path/queries", pattern: "*.fasta" , mode: 'copy'
    publishDir "$output_path", pattern: "genomes_text.txt" , mode: 'copy'

    input: 
    val table
    val genomes
    val gffs
    val output_path

    output: 
    path "genomes_text.txt", emit: file_text 
    path "*_queries.fasta", emit: queries
    path "*_promoter_regions.fasta", emit: promoters
    path "tss_analyzed.tsv", emit: summary_transcripts
    path "*.tsv", emit: analyzed_tss

    script:
    """
        python3 $params.pyScriptTableToQuery $table $genomes $gffs
    """
}

/*
    Extracts the corresponding TaxIDs from the Taxonomy tree of NCBI using the NCBI-IDs of the given genomes (as a single call)
*/
process GETGENOMESLCA {
    label "DataPreparation"
    tag "$ids"

    input: 
    val ids

    output: 
    stdout emit: tax_id

    script:
    """
        efetch -db nuccore -id $ids -format docsum | xtract -pattern DocumentSummary -element TaxId 
    """
}

/*
    Extracts the LCA from the species in question using the given TaxIDs
*/
process GETLCAID {
    label "DataPreparation"

    input: 
    val ids

    output: 
    path "common_species_id.txt", emit: common_species_id

    script: 
    """
        python3 $params.pyScriptCommonSpecies $ids > common_species_id.txt
    """
}

/*
    Creates a TaxIDlist for BLAST with all corresponding species under the extracted LCA
*/
process GETBLASTID {
    label "DataPreparation"

    input: 
    val taxidlist

    output: 
    path "taxidlist.taxid", emit: taxidlist_file

    script:
    """
        get_species_taxids.sh -t $taxidlist > taxidlist.taxid
    """
}

/*
    Runs BLAST against the nt-database using the extracted queries and the restricted TaxIDs
*/
process BLASTFASTA {
    label "DataPreparation"

	input:
    each query
    path taxidlist
    env BLASTDB
    val output_path

    output:
    path "*.blastn", emit: blasted_files

    script: 
    """
        blastn -query $query -db nt -out ${query.baseName - "_queries"}.blastn -outfmt "6 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore pident frames qcovhsp" -task dc-megablast -taxidlist $taxidlist
    """
}

/*
    Evaluates the BLAST-hits and extracts the corresponding hits for the pairwise alignment (for QRNA)
*/
process EVALBLAST {
    label "DataPreparation"
    publishDir "$output_path/filtered_blast", mode: 'copy'

    input: 
    path blasted_files
    val output_path

    output:
    path "*.tsv", emit: filtered_queries_ch

    // path "evaluation*" into evaluation_ch
    // TODO: Still add the evaluation table, since there might be transcripts without any hit. 
    script:
    """
        python3 $params.pyEvaluateBlast $blasted_files -t 10
    """
}
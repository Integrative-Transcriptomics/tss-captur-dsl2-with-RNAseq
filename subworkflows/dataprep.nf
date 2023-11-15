include { MASTERTOFASTA, GETGENOMESLCA, GETLCAID, GETBLASTID, BLASTFASTA, EVALBLAST } from '../modules/dataprep'

workflow DATAPREPARATION { 
    take:
        table_ch
        genomes_path
        gff_path
        blastDB
        output_path

    main:
        MASTERTOFASTA(table_ch, genomes_path, gff_path, output_path)
        GETGENOMESLCA(MASTERTOFASTA.out.file_text.map{ it.getText().split("\n") }.flatten().collect().map{ it.join(",") })
        GETLCAID(GETGENOMESLCA.out.tax_id.map{ it.split("\n").join(" ") })
        GETBLASTID(GETLCAID.out.common_species_id.readLines().map{ it[0].findAll(/\d{1,10}/)[0] })
        BLASTFASTA(MASTERTOFASTA.out.queries, GETBLASTID.out.taxidlist_file, blastDB, output_path)
        EVALBLAST(BLASTFASTA.out.blasted_files, output_path)

    emit:
        queries = MASTERTOFASTA.out.queries
        promoters = MASTERTOFASTA.out.promoters
        summary_transcripts = MASTERTOFASTA.out.summary_transcripts
        filtered_queries = EVALBLAST.out.filtered_queries_ch
}
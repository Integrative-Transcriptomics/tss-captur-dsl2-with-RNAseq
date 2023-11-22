include { MASTERTOFASTA; GETGENOMESLCA; GETLCAID; GETBLASTID; BLASTFASTA; EVALBLAST } from '../modules/dataprep'

workflow DATAPREPARATION { 
    take:
        masterTable
        genomesPath
        gffPath
        outputPath

    main:
        MASTERTOFASTA(masterTable, genomesPath, gffPath, outputPath)
        GETGENOMESLCA(MASTERTOFASTA.out.fileText.map{ it.getText().split("\n") }.flatten().collect().map{ it.join(",") })
        GETLCAID(GETGENOMESLCA.out.taxID.map{ it.split("\n").join(" ") })
        GETBLASTID(GETLCAID.out.commonSpeciesID.readLines().map{ it[0].findAll(/\d{1,10}/)[0] })
        BLASTFASTA(MASTERTOFASTA.out.queries, GETBLASTID.out.taxIDListFile, params.blastDB, outputPath)
        EVALBLAST(BLASTFASTA.out.blastFiles, outputPath)

    emit:
        queries = MASTERTOFASTA.out.queries
        promoters = MASTERTOFASTA.out.promoters
        summaryTranscripts = MASTERTOFASTA.out.summaryTranscripts
        blastFiltered = EVALBLAST.out.blastFiltered
}
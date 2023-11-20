#!/usr/bin/env nextflow
import java.nio.file.Paths

// Input Manager
table_ch =      Channel.value("${projectDir}/${params.inputTable}")
genomes_path =  Channel.value("${projectDir}/${params.inputGenomes}")
genomes_ext =   Channel.fromPath(["fa", "fna", "fasta", "frn", "faa", "ffn"].collect{ "${params.inputGenomes}*.${it}" })
gff_path =      Channel.value("${projectDir}/${params.inputGFFs}")
output_path =   Channel.value("${projectDir}/${params.outputDir}")

/*
// TODO: Help message
log.info 
"""\
            TSS-CAPTUR
 ===================================
 MasterTable : ${table_ch}
 Genomes path: ${genomes_path}
 GFFs path   : ${gff_path}
 Output path : ${output_path}
 """
 */

include { DATAPREPARATION } from './subworkflows/dataprep'
include { MEME } from './modules/meme'
include { CLASSIFICATION } from './subworkflows/classification'
include { TERMINATORPREDICTION } from './subworkflows/terminatorpred'
include { RNAFOLD } from './modules/rnafold'
include { CREATEREPORT } from './modules/report'

workflow {
    DATAPREPARATION(table_ch, genomes_path, gff_path, output_path)
    MEME(DATAPREPARATION.out.promoters, output_path)
    CLASSIFICATION(DATAPREPARATION.out.queries, DATAPREPARATION.out.filtered_queries.collect(), output_path)
    TERMINATORPREDICTION(genomes_ext,
                        projectDir, 
                        CLASSIFICATION.out.crd_files, 
                        DATAPREPARATION.out.summary_transcripts, 
                        output_path)
    RNAFOLD(TERMINATORPREDICTION.out.allocation, genomes_path, output_path)
    CREATEREPORT(RNAFOLD.out.output_figures.collect(), MEME.out.motifResult.collect(), output_path)
}
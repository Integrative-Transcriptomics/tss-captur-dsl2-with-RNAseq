#!/usr/bin/env nextflow

// Info message
log.info """\

            TSS-CAPTUR
 ===================================
 MasterTable : ${params.masterTable}
 Genomes path: ${params.genomesPath}
 GFFs path   : ${params.gffPath}
 Forward/reverse Wiggle paths: "forward: ${params.forwardWigglesPath}  reverse: ${params.reverseWigglesPath}"
 Output path : ${params.outputPath}
 """
 
include { DATAPREPARATION } from './subworkflows/dataprep'
include { MEME } from './modules/meme'
include { CLASSIFICATION } from './subworkflows/classification'
include { TERMINATORPREDICTION } from './subworkflows/terminatorpred'
include { RNAFOLD } from './modules/rnafold'
include { CREATEREPORT } from './modules/report'
include { CLEANWORKDIR } from './modules/cleaner'
include { TERMINATORALLOC } from './subworkflows/terminatoralloc.nf'
include { WIGGLEANALYSIS } from './subworkflows/WiggleAnalysis'
 
workflow {
    genomesExt = Channel.fromPath(["fa", "fna", "fasta", "frn", "faa", "ffn"].collect{ "${params.genomesPath}*.${it}" })

    DATAPREPARATION(params.masterTable, params.genomesPath, params.gffPath, params.outputPath)
    MEME(DATAPREPARATION.out.promoters, params.outputPath)
    CLASSIFICATION(DATAPREPARATION.out.queries, DATAPREPARATION.out.blastFiltered.collect(), params.outputPath)
    TERMINATORPREDICTION(genomesExt,
                        projectDir,  
                        params.outputPath)
    
    if(!(params.forwardWigglesPath == "Gar nix" || params.forwardWigglesPath == "" || params.reverseWigglesPath == "Gar nix" || params.reverseWigglesPath == ""))
    {
         WIGGLEANALYSIS(params.forwardWigglesPath, params.reverseWigglesPath, params.gffPath, params.genomesPath, projectDir, TERMINATORPREDICTION.out.gffNocornac, TERMINATORPREDICTION.out.gffRhoterm, params.masterTable)
         TERMINATORALLOC(TERMINATORPREDICTION.out.gffNocornac.collect(),
                    TERMINATORPREDICTION.out.gffRhoterm.collect(),
                    CLASSIFICATION.out.crdFiles, 
                    DATAPREPARATION.out.summaryTranscripts,
                    WIGGLEANALYSIS.out.wiggleTerms,
                    params.outputPath)
    }
    else
    {
         TERMINATORALLOC(TERMINATORPREDICTION.out.gffNocornac.collect(),
                    TERMINATORPREDICTION.out.gffRhoterm.collect(),
                    CLASSIFICATION.out.crdFiles, 
                    DATAPREPARATION.out.summaryTranscripts,
                    params.nofile,
                    params.outputPath)
    }
    
    RNAFOLD(TERMINATORALLOC.out.allocation, params.genomesPath, params.outputPath)

    CREATEREPORT(RNAFOLD.out.outputFigures.collect(), MEME.out.motifResult.collect(), params.outputPath) | collect | CLEANWORKDIR
}

workflow.onComplete = {
    
}

// Writes the error message to the output path
workflow.onError = {
    new File(params.errorPath).text = workflow.errorReport
}
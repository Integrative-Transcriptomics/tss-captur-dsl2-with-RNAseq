include { FINDTERMINATORS } from '../modules/terminatorpred'

workflow TERMINATORALLOC {
    take:
        nocornacgffs
        rhotermgffs
        crdFiles
        summaryTranscripts 
        wigglePredictedTerminators
        outputPath
    
    main:
        FINDTERMINATORS(nocornacgffs, 
                        rhotermgffs, 
                        crdFiles.collect(), 
                        summaryTranscripts,
                        wigglePredictedTerminators, 
                        outputPath)
    
    emit:
        allocation = FINDTERMINATORS.out.terminatorsAllocation
}
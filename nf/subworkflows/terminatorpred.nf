include { RHOTERMPREDICT; NOCORNAC; FINDTERMINATORS } from '../modules/terminatorpred'

workflow TERMINATORPREDICTION {
    take:
        genomesExt
        projectDir
        crdFiles
        summaryTranscripts 
        outputPath
    
    main:
        RHOTERMPREDICT(genomesExt, outputPath)
        NOCORNAC(params.nocornacConfig, projectDir, genomesExt, outputPath)
        FINDTERMINATORS(NOCORNAC.out.gffNocornac.collect(), 
                        RHOTERMPREDICT.out.gffRhoterm.collect(), 
                        crdFiles.collect(), 
                        summaryTranscripts, 
                        outputPath)
    
    emit:
        allocation = FINDTERMINATORS.out.terminatorsAllocation
}
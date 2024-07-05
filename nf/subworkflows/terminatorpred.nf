include { RHOTERMPREDICT; NOCORNAC; FINDTERMINATORS } from '../modules/terminatorpred'

workflow TERMINATORPREDICTION {
    take:
        genomesExt
        projectDir
        outputPath
    
    main:
        RHOTERMPREDICT(genomesExt, outputPath)
        NOCORNAC(params.nocornacConfig, projectDir, genomesExt, outputPath)
    
    emit:
        gffRhoterm = RHOTERMPREDICT.out.gffRhoterm
        gffNocornac = NOCORNAC.out.gffNocornac
}
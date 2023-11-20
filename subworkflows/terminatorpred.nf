include { RHOTERMPREDICT; NOCORNAC; FINDTERMINATORS } from '../modules/terminatorpred'

workflow TERMINATORPREDICTION {
    take:
        genomes_ext
        nocornacConfig
        projectDir
        crd_files
        summary_transcripts 
        output_path
    
    main:
        RHOTERMPREDICT(genomes_ext, output_path)
        NOCORNAC(nocornacConfig, projectDir, genomes_ext, output_path)
        FINDTERMINATORS(NOCORNAC.out.nocornac_gffs.collect(), 
                        RHOTERMPREDICT.out.rhoterm_gffs.collect(), 
                        crd_files.collect(), 
                        summary_transcripts, 
                        output_path)
    
    emit:
        allocation = FINDTERMINATORS.out.terminators_allocation
}
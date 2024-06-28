include { WIGGLESCORETERMINATORS; SIZESFROMFASTA; WIGTOBEDGRAPH; BEDGRAPHTOBIGWIG;} from '../modules/wiggleanalysis'

workflow WIGGLEANALYSIS {
    take:
        wigglePath
        annotationPath
        fastaPath
        projectDir
        terminatorsAllocation
        MasterTable

    main:
        //build sizes file from fasta
        //convert wigs to big wigs (get length of chromosome from fasta)
        SIZESFROMFASTA(fastaPath)
        WIGTOBEDGRAPH(wigglePath, SIZESFROMFASTA.out.sizesFile)
        BEDGRAPHTOBIGWIG(WIGTOBEDGRAPH.out.bgFile, SIZESFROMFASTA.out.sizesFile)

        //score bigwigs using terminators from rhotermpred and trasntermhp, output scoring to file
        WIGGLESCORETERMINATORS(BEDGRAPHTOBIGWIG.out.bwFile, annotationPath, terminatorsAllocation, MasterTable)

    emit:
        //emit output file with terminators and scoring here
        wiggleTerms = WIGGLESCORETERMINATORS.out.RescoredTerminators
}
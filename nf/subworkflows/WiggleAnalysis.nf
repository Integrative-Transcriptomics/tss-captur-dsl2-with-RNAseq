include { 
    WIGGLESCORETERMINATORS;
    SIZESFROMFASTA;
    WIGTOBEDGRAPH as ForwardWGTB; 
    WIGTOBEDGRAPH as ReverseWGTB;
    BEDGRAPHTOBIGWIG as ForwardBGTBW;
    BEDGRAPHTOBIGWIG as ReverseBGTBW;} from '../modules/wiggleanalysis'

workflow WIGGLEANALYSIS {
    take:
        forwardWigglespath
        reverseWigglesPath
        annotationPath
        fastaPath
        projectDir
        gffNocornac
        gffRhoterm
        MasterTable

    main:
        //build sizes file from fasta
        //convert wigs to big wigs (get length of chromosome from fasta)
        SIZESFROMFASTA(fastaPath)

        ForwardWGTB(forwardWigglespath, SIZESFROMFASTA.out.sizesFile, "forward")
        ForwardBGTBW(ForwardWGTB.out.bgFile, SIZESFROMFASTA.out.sizesFile, "forward")

        ReverseWGTB(reverseWigglesPath, SIZESFROMFASTA.out.sizesFile, "reverse")
        ReverseBGTBW(ReverseWGTB.out.bgFile, SIZESFROMFASTA.out.sizesFile, "reverse")

        //score bigwigs using terminators from rhotermpred and trasntermhp, output scoring to file
        WIGGLESCORETERMINATORS(ForwardBGTBW.out.bwFile, ReverseBGTBW.out.bwFile, annotationPath, gffNocornac, gffRhoterm, MasterTable)

    emit:
        //emit output file with terminators and scoring here
        wiggleTerms = WIGGLESCORETERMINATORS.out.RescoredTerminators
}
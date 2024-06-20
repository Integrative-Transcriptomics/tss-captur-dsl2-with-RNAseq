
workflow WIGGLEANALYSIS{
    take:
        wigglePath
        annotationPath
        fastaPath
        projectDir
        terminatorsAllocation
        MasterTable

    main:
        //build sizes file from fasta
        //convert wigs to big wigs (get length of chromosome from gff)
        WIGGLESTOBWS(wigglePath, fastaPath, projectDir)

        //score bigwigs using terminators from rhotermpred and trasntermhp, output scoring to file
        WIGGLESCORETERMINATORS(WIGGLESTOBWS.out.bwFile, annotationPath, terminatorsAllocation,  MasterTable)

    emit:
        //emit output file with terminators and scoring here
}

process WIGGLESCORETERMINATORS
{
    input:
    path bigWigFile
    path annotationPath
    path TerminatorAlloc
    path MasterTable

    output
    path ScoredTermGFF, emit RescoredTerminators

    script:
    """
    python $params.pyDerivScoring ${bigWigFile} ${annoationPath} ${TerminatorAlloc} ${MasterTable}
    """

}

process WIGGLESTOBWS
{
    input:
    path wigglePath
    path fastaPath
    path bigWigOutputPath

    output:
    path bwPath, emit: bwFile

    script:
    """
        sizes_file = $(python $params.pySizesFileFromFasta $fastaPath)
        bg_file = $(python $params.pyManualWigToBedGraph $wigglePath $sizes_file CapturNextflow)
        $params.bedGraphToBigWig $bg_file $sizes_file $bigWigOutputPath
    """
}
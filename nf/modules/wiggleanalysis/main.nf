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
process WIGGLESCORETERMINATORS
{
    input:
    path bigWigFile
    path annotationPath
    path TerminatorAlloc
    path MasterTable

    output:
    path "test", emit: RescoredTerminators

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
        python3 $params.pySizesFileFromFasta $fastaPath $bigWigOutputPath
        python3 $params.pyManualWigToBedGraph $wigglePath $bigWigOutputPath/autoSizeFile.sizes CapturNextflow
        $params.bedGraphToBigWig CapturNextflow.bedGraph $bigWigOutputPath/autoSizeFile.sizes $bigWigOutputPath
    """
}
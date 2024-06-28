process WIGGLESCORETERMINATORS
{
    container 'testdocker'

    input:
    path bigWigFile
    path annotationPath
    path TerminatorAlloc
    path MasterTable

    output:
    path "test", emit: RescoredTerminators

    script:
    """
    python $params.pyDerivScoring $bigWigFile $annotationPath $TerminatorAlloc $MasterTable
    """
}

// process WIGGLESTOBWS
// {
//     input:
//     path wigglePath
//     path fastaPath
//     path bigWigOutputPath

//     output:
//     path bwPath, emit: bwFile



//     script:
//     """
//         python3 $params.pySizesFileFromFasta $fastaPath $bigWigOutputPath
//         python3 $params.pyManualWigToBedGraph $wigglePath $bigWigOutputPath/autoSizeFile.sizes CapturNextflow
//         $params.bedGraphToBigWig CapturNextflow.bedGraph $bigWigOutputPath/autoSizeFile.sizes $bigWigOutputPath
//     """
// }

process SIZESFROMFASTA
{
    container 'testdocker'

    input:
    path fastaPath

    output:
    path "*.sizes", emit: sizesFile

    script:
    """
        python3 $params.pySizesFileFromFasta $fastaPath \$PWD
    """

}

process WIGTOBEDGRAPH
{
    container 'testdocker'
    
    input:
    path wigglePath
    path sizesFile

    output: 
    path "*.bedGraph", emit: bgFile

    script:
    """
        python3 $params.pyManualWigToBedGraph $wigglePath $sizesFile \$PWD
    """

}

process BEDGRAPHTOBIGWIG
{
    container 'testdocker'

    input:
    file bgFile
    file sizesFile

    output: 
    path "*.bw", emit: bwFile

    script:
    """
        $params.bedGraphToBigWig $bgFile $sizesFile \$PWD/autoBW.bw
    """
}
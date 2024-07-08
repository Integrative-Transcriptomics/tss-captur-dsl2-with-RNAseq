process WIGGLESCORETERMINATORS
{
    container 'testdocker'

    input:
    path forwardBigWig
    path reverseBigWig
    path annotationPath
    path gffNocornac
    path gffRhoterm
    path MasterTable

    output:
    path "*.tsv", emit: RescoredTerminators

    script:
    """
    python $params.pyDerivScoring $forwardBigWig $reverseBigWig $annotationPath $gffNocornac $gffRhoterm $MasterTable
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
    val fileName

    output: 
    path "${fileName}.bedGraph", emit: bgFile

    script:
    """
        python3 $params.pyManualWigToBedGraph $wigglePath $sizesFile $fileName
    """

}

process BEDGRAPHTOBIGWIG
{
    container 'testdocker'

    input:
    file bgFile
    file sizesFile
    val fileName

    output: 
    path "${fileName}.bw", emit: bwFile

    script:
    """
        $params.bedGraphToBigWig $bgFile $sizesFile ${fileName}.bw
    """
}
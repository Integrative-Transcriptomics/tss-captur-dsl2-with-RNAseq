/*
    Runs a Motif Analysis on the promoter region using MEME
    Returns the output of MEME and a .tsv file for data parsing
*/
process MEME {
    label params.memeLabel
    publishDir "$outputPath/MotifAnalysis/${promoter.baseName - "_promoter_regions"}", mode: params.pubDirMode

    input: 
    each file(promoter)
    val outputPath

    output:
    path "*", emit: motifResult

    script:
    """
        meme $promoter -dna -nmotifs $params.motifNumber -minw 5 -maxw 20 
        python3 $params.pyMemeParser --meme meme_out/meme.xml --genome ${promoter.baseName - "_promoter_regions"}
    """
}
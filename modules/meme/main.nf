/*
    Runs a Motif Analysis on the promoter region using MEME
    Returns the output of MEME and a .tsv file for data parsing
*/
process MEME {
    label "MotifAnalysis"
    publishDir "$output_path/MotifAnalysis/${promoter.baseName - "_promoter_regions"}", mode: 'copy'

    input: 
    each file(promoter)
    val output_path

    output:
    path "*", emit: motifResult

    script:
    """
    meme $promoter -dna -nmotifs $params.motifNumber -minw 5 -maxw 20 
    python3 $pyMemeParser --meme meme_out/meme.xml --genome ${promoter.baseName - "_promoter_regions"}
    """
}
/*
    Runs a MFE analysis on the transcripts classified as RNAs
    Returns a .jpg file to visualize the secondary structure of the transcripts
*/
process RNAFOLD {
    label params.rnaFoldLabel
    publishDir "$outputPath/SecondaryStructure/${terminator.baseName - "_allocated_terminators"}/Visualizations", pattern: "*.jpg", mode: params.pubDirMode
    publishDir "$outputPath/SecondaryStructure/${terminator.baseName - "_allocated_terminators"}/", pattern: "*.tsv", mode: params.pubDirMode

    input:
    each terminator
    val genomes
    val outputPath

    output:
    path "*.jpg" 
    path "*.tsv" 
    path "*.fasta"
    stdout emit: outputFigures

    script:
    """
        python3 $params.pyExtractRNATranscripts --rnas $terminator --genome_path $genomes
        $params.rnafold --noLP -i *.fasta > rnaFold.out
        gawk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""), \$0);N++;next;} {match(\$0, /(.*)\\s\\((.*)\\)/, ary); if (length(ary)>0) {printf("%s\\t%s",ary[1],ary[2]);} else { printf("%s\\t", \$0)}} END {printf("\\n");}' rnaFold.out > ${terminator.baseName - "_allocated_terminators"}.tsv
        gmt psconvert *.ps
    """
}
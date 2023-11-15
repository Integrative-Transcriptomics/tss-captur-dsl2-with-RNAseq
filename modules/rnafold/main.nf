/*
    Runs a MFE analysis on the transcripts classified as RNAs
    Returns a .jpg file to visualize the secondary structure of the transcripts
*/
process RNAFOLD {
    label "SecondaryStruct"
    publishDir "$output_path/SecondaryStructure/${terminator.baseName - "_allocated_terminators"}/Visualizations", pattern: "*.jpg", mode: 'copy'
    publishDir "$output_path/SecondaryStructure/${terminator.baseName - "_allocated_terminators"}/", pattern: "*.tsv", mode: 'copy'
    // TODO: publishDir "$output_path/SecondaryStructure/${terminator.baseName - "_allocated_terminators"}/Transcripts", pattern: "*.fasta", mode: 'copy'

    input:
    each terminator
    val genomes
    val output_path

    output:
    path "*.jpg" 
    path "*.tsv" 
    path "*.fasta"
    stdout emit: output_figures

    """
    python3 $pyExtractRNATranscripts --rnas $terminator --genome_path $genomes
    $rnafold --noLP -i *.fasta > rnaFold.out
    gawk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""), \$0);N++;next;} {match(\$0, /(.*)\\s\\((.*)\\)/, ary); if (length(ary)>0) {printf("%s\\t%s",ary[1],ary[2]);} else { printf("%s\\t", \$0)}} END {printf("\\n");}' rnaFold.out > ${terminator.baseName - "_allocated_terminators"}.tsv
    gmt psconvert *.ps
    """
}
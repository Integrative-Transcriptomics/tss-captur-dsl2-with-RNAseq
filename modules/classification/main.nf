/*
    Runs CNIT for the classification
*/
process CNIT {
    label "Classification"

    input: 
    each val(query)

    output:
    path "${query.baseName}.index", emit: cnit_results

    shell:
    """
    python2 $cnit -f $query -o ${query.baseName} -m "pl"
    mv ${query.baseName}/${query.baseName}.index .
    """
}

/*
    Evaluates output of CNIT and converts it to a .tsv file 
*/
process EVALCNIT {
    label "Classification"
    publishDir "$output_path/Classification/CNIT/Evaluation"

    input: 
    path cnit_index
    val output_path

    output:
    path "*.tsv", emit: cnit_eval

    shell:
    """
    python3 $pyFromCNITtoTSV --cnit_list $cnit_index
    """
}

/*
    Runs QRNA for the classification
*/
process QRNA {
    label "Classification"
    cache "lenient"
    // TODO: publishDir "$output_path/Classification/QRNA/Results", mode: 'copy' // Put in one folder QRNA/result

    input: 
    each blasted_file
    env QRNADB
    val output_path

    output:
    path "${blasted_file.baseName}.qrna", emit: qrna_normal

    """
    gawk 'BEGIN { OFS = "\\n"} {print ">" \$1 "|" \$2 "-" \$3 "|" \$9"," \$10"," \$11 ", Frames:" \$12", Score:" \$17, \$4, ">" \$5 "|" \$6 "-" \$7, \$8 "\\n"}' $blasted_file > ${blasted_file.baseName}.fa
    $eqrna --ones ${blasted_file.baseName}.fa > ${blasted_file.baseName}.qrna
    """
}

/*
    Transcribes the output of QRNA to a .tsv file and corrects the position of each UTR region using the BLAST files
*/
process EVALQRNA{
    label "Classification"
    publishDir "$output_path/Classification/QRNA/Evaluation", pattern: "*.tsv", mode: 'copy'

    input:
    path path_qrna_files
    path path_filtered_blast_files
    val output_path

    output:
    path "*.tsv", emit: qrna_eval

    """
        python3 $pyFromQRNAtoTSV --qrna_list $path_qrna_files --filteredBlast_list $path_filtered_blast_files
    """
}

/*
    Compares the output of CNIT and QRNA to call the final classification on each transcript
*/
process COMPARECNITQRNA {
    label "Classification"
    publishDir "$output_path/Classification",  mode: 'copy'

    input:
    path qrna_eval_file
    path cnit_eval_file
    val output_path

    output:
    path "*final_classification.tsv", emit: final_eval_classification
    path "*.crd", emit: crd_files
    path "*.gff"

    shell:     
    """
        python3 $pyDecideClass --qrna $qrna_eval_file --cnit $cnit_eval_file
    """
}
/*
    Runs CNIT for the classification
*/
process CNIT {
    label params.classificationLabel

    input: 
    each query

    output:
    path "${query.baseName}.index", emit: cnitResults

    script:
    """
        python2 $params.cnit -f $query -o ${query.baseName} -m "pl"
        mv ${query.baseName}/${query.baseName}.index .
    """
}

/*
    Evaluates output of CNIT and converts it to a .tsv file 
*/
process EVALCNIT {
    label params.classificationLabel
    publishDir "$outputPath/Classification/CNIT/Evaluation"

    input: 
    path cnitIndex
    val outputPath

    output:
    path "*.tsv", emit: cnitEval

    script:
    """
        python3 $params.pyFromCNITtoTSV --cnit_list $cnitIndex
    """
}

/*
    Runs QRNA for the classification
*/
process QRNA {
    label params.classificationLabel
    cache "lenient"

    input: 
    each blastFilteredFile
    env QRNADB
    val outputPath

    output:
    path "${blastFilteredFile.baseName}.qrna", emit: qrnaNormal

    script:
    """
        gawk 'BEGIN { OFS = "\\n"} {print ">" \$1 "|" \$2 "-" \$3 "|" \$9"," \$10"," \$11 ", Frames:" \$12", Score:" \$17, \$4, ">" \$5 "|" \$6 "-" \$7, \$8 "\\n"}' $blastFilteredFile > ${blastFilteredFile.baseName}.fa
        $params.eqrna --ones ${blastFilteredFile.baseName}.fa > ${blastFilteredFile.baseName}.qrna
    """
}

/*
    Transcribes the output of QRNA to a .tsv file and corrects the position of each UTR region using the BLAST files
*/
process EVALQRNA{
    label params.classificationLabel
    publishDir "$outputPath/Classification/QRNA/Evaluation", pattern: "*.tsv", mode: params.pubDirMode

    input:
    path qrnaFiles
    path blastFiltered
    val outputPath

    output:
    path "*.tsv", emit: qrnaEval
    
    script:
    """
        python3 $params.pyFromQRNAtoTSV --qrna_list $qrnaFiles --filteredBlast_list $blastFiltered
    """
}

/*
    Compares the output of CNIT and QRNA to call the final classification on each transcript
*/
process COMPARECNITQRNA {
    label params.classificationLabel
    publishDir "$outputPath/Classification",  mode: params.pubDirMode

    input:
    path qrnaEvalFile
    path cnitEvalFile
    val outputPath

    output:
    path "*final_classification.tsv", emit: finalEvalClassification
    path "*.crd", emit: crdFiles
    path "*.gff"

    script:     
    """
        python3 $params.pyDecideClass --qrna $qrnaEvalFile --cnit $cnitEvalFile
    """
}
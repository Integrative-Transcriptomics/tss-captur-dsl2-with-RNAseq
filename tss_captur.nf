#!/usr/bin/env nextflow
import java.nio.file.Paths

params.inputTable =  "/input/CampySubset/MasterTableSubset.tsv" 
params.inputGenomes = "/input/CampySubset/fasta/" 
params.inputGFFs = "/input/CampySubset/gff" 
params.motifNumber = 5 

params.outputDir = "/output/CampySubset_correctedThresholdNocornac_NotCross"
params.conditions = "FALSE"

// Input Manager
nameProjDir = projectDir.toString()
table_ch = Channel.value(Paths.get(nameProjDir, params.inputTable))
genomes_path = Channel.value(Paths.get(nameProjDir, params.inputGenomes))
genomes_ext = ["fa", "fna", "fasta", "frn", "faa", "ffn"].collect{"${params.inputGenomes}/*.${it}"}
gff_path = Channel.value(Paths.get(nameProjDir,params.inputGFFs))
output_path = Channel.value(Paths.get(nameProjDir, params.outputDir))

// Own scripts
pythonScriptTableToQuery=Paths.get(nameProjDir, "/bin/scripts/tableToSeqsNotCross.py")
pyScriptCommonSpecies = Paths.get(nameProjDir,"/bin/scripts/processCommonIDs.py")
pyEvaluateBlast = Paths.get(nameProjDir, "/bin/scripts/processBlastFileNew.py")
pyFromQRNAtoTSV =  Paths.get(nameProjDir, "/bin/scripts/fromQRNAtoTSV.py")
pyFromCNITtoTSV =  Paths.get(nameProjDir, "/bin/scripts/fromCNITtoTSV.py")
pyTransformRhoTermToGFF =  Paths.get(nameProjDir, "/bin/scripts/fromRhoTermToGFF.py")
pyDecideClass =  Paths.get(nameProjDir, "/bin/scripts/fromTablesToCRD.py")
pyAllocateTerminators = Paths.get(nameProjDir, "/bin/scripts/allocateTerminators.py")
pyExtractRNATranscripts = Paths.get(nameProjDir, "/bin/scripts/extractTranscripts.py")
pyMemeParser = Paths.get(nameProjDir, "/bin/scripts/memeParser.py")
pyCreateReport = Paths.get(nameProjDir, "/bin/app.py")


// Programs used
eqrna =  Paths.get(nameProjDir,"/bin/eqrna")
templates =  Paths.get(nameProjDir,"/bin/templates")
staticfiles =  Paths.get(nameProjDir,"/bin/static")
eqrnaLib =  Paths.get(nameProjDir,"/bin/lib")
nocornac =  Paths.get(nameProjDir, "/bin/nocornac.jar")
nocornacConfig =  Paths.get(nameProjDir,"/bin/config.conf")
cnit = Paths.get(nameProjDir, "/bin/CNCI2.py")
pyRhoTermPredict = Paths.get(nameProjDir, "/bin/RhoTermPredict.py")
rnafold = Paths.get(nameProjDir, "/bin/RNAfold")



Channel.fromPath(genomes_ext).into { genomes_nocornac; genomes_rhoterm }


/*
    Extracts transcripts from the MasterTable of TSSpredator
    Returns queries for classification, promoter regions, a .txt file with genome names, and a .tsv file with ignored and analyzed positions
*/
process MASTERTOFASTA {
    label "DataPreparation"
    publishDir "$output_path/queries", pattern: "*.tsv", mode: 'copy'
    publishDir "$output_path/queries", pattern: "*.fasta" , mode: 'copy'
    publishDir "$output_path", pattern: "genomes_text.txt" , mode: 'copy'

    input: 
    val table
    val genomes
    val gffs
    val output_path

    output: 
    path "genomes_text.txt", emit: file_text 
    path "*_queries.fasta", emit: queries
    path "*_promoter_regions.fasta", emit: promoters
    path "tss_analyzed.tsv", emit: summary_transcripts
    path "*.tsv", emit: analyzed_tss

    shell:
    """
        python3 $pythonScriptTableToQuery $table $genomes $gffs
    """
}

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

    shell:
    """
    meme $promoter -dna -nmotifs $params.motifNumber -minw 5 -maxw 20 
    python3 $pyMemeParser --meme meme_out/meme.xml --genome ${promoter.baseName - "_promoter_regions"}
    """
}

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

cnit_results.collect().set{cnit_collected}

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

file_text.map{it.getText().split("\n")}.set{ fixed_genomes } 
genomes_to_process = fixed_genomes.flatten().collect()
genomes_to_process.map{it.join(",")}.set{acc_codes_to_id}

/*
    Extracts the corresponding TaxIDs from the Taxonomy tree of NCBI using the NCBI-IDs of the given genomes (as a single call)
*/
process GETGENOMESLCA {
    label "DataPreparation"
    tag "$ids"

    input: 
    val ids

    output: 
    stdout, emit: tax_id

    shell:
    """
        efetch -db nuccore -id $ids -format docsum | xtract -pattern DocumentSummary -element TaxId 
    """
}

tax_id.map{it.split("\n").join(" ")}.set{ tax_id_test } 

/*
    Extracts the LCA from the species in question using the given TaxIDs
*/
process GETLCAID {
    label "DataPreparation"

    input: 
    val ids

    output: 
    path "common_species_id.txt", emit: common_species_id

    shell: 
    """
        python3 $pyScriptCommonSpecies $ids > common_species_id.txt
    """
}

// Adding the RegEx avoids problem when the Taxa Id is downloaded for the first time. 
common_species_id.readLines().map{it[0].findAll(/\d{1,10}/)[0]}.set{cut_common_species_id}

/*
    Creates a TaxIDlist for BLAST with all corresponding species under the extracted LCA
*/
process GETBLASTID {
    label "DataPreparation"

    input: 
    val taxidlist

    output: 
    path "taxidlist.taxid", emit: taxidlist_file

    shell:
    """
        get_species_taxids.sh -t $taxidlist > taxidlist.taxid
    """
}

/*
    Runs BLAST against the nt-database using the extracted queries and the restricted TaxIDs
*/
process BLASTFASTA {
    label "DataPreparation"

	input:
    each val(query)
    path taxidlist
    env BLASTDB
    val output_path

    output:
    path "*.blastn", emit: blasted_files

    shell: 
    """
    blastn -query $query -db nt -out ${query.baseName - "_queries"}.blastn  -outfmt "6 qseqid qstart qend qseq sseqid sstart send  sseq evalue bitscore pident frames qcovhsp" -task dc-megablast  -taxidlist $taxidlist
    """
}

/*
    Evaluates the BLAST-hits and extracts the corresponding hits for the pairwise alignment (for QRNA)
*/
process EVALBLAST {
    label "DataPreparation"
    publishDir "$output_path/filtered_blast", mode: 'copy'

    input: 
    path blasted_files
    val output_path

    output:
    path "*.tsv", emit: filtered_queries_ch

    // TODO: Still add the evaluation table, since there might be transcripts without any hit. 
    """
        python3 $pyEvaluateBlast $blasted_files -t 10
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

qrna_normal.collect().set{qrna_evaluate}
filtered_queries_for_qrna_check.collect().set{blasted_files_qrna_check}

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

cnit_eval.collect().set{cnit_eval_collected}
qrna_eval.collect().set{qrna_eval_collected}

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

/*
    Identifies possible Rho-Dependent terminators for the found genomes
*/
process RHOTERMPREDICT {
    label "TerminatorPrediciton"
    publishDir "$output_path/Terminators/RhoTermPredict",  mode: 'copy'

    input:
    path genome
    val output_path

    output:
    path "*.gff", emit: rhoterm_gffs
    
    shell:
    """
        python3 $pyRhoTermPredict $genome --output ${genome.baseName}_rhoterm
        python3 $pyTransformRhoTermToGFF --rhoterm *${genome.baseName}_rhoterm.tsv --genome ${genome.baseName}
    """
}

/*
    Identifies Rho-Independent terminators for the found genomes using TransTermHP
*/
process NOCORNAC {
    label "TerminatorPrediciton"
    publishDir "$output_path/Terminators/nocoRNAc",  mode: 'copy'

    input:
    val config_file
    val proj_dir
    path genome
    val output_path

    output: 
    path "*.gff", emit: nocornac_gffs

    """
    sed 's@dataPath = data@dataPath =  $proj_dir/nocornac/data@g' $config_file > config_temp.conf
    sed 's@transtermPath = progs/@transtermPath =  $proj_dir/bin/progs/@g' config_temp.conf > config.conf
    java -Xmx1G -jar $nocornac  -genomeFastaFile $genome -gffOutFile ${genome.baseName}_nocornac.gff -terminators
    """
}

nocornac_gffs.collect().set{nocornac_collected}
rhoterm_gffs.collect().set{rhoterm_collected}
crd_files.collect().set{crd_files_collected}

/*
    Finds the possible terminators for each transcript
*/
process FINDTERMINATORS {
    label "TerminatorPrediciton"
    publishDir "$output_path/Terminators",   mode: 'copy'

    input:
    path nocornac
    path rhoterm
    path crd
    path tsv
    val output_path

    output:
    path "*.tsv", emit: terminators_allocation 

    """
    python3 $pyAllocateTerminators --nocornac $nocornac --rhoterm $rhoterm --crd $crd --tssAnalyzed $tsv
    """
}

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
    stdout, emit: output_figures

    """
    python3 $pyExtractRNATranscripts --rnas $terminator --genome_path $genomes
    $rnafold --noLP -i *.fasta > rnaFold.out
    gawk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""), \$0);N++;next;} {match(\$0, /(.*)\\s\\((.*)\\)/, ary); if (length(ary)>0) {printf("%s\\t%s",ary[1],ary[2]);} else { printf("%s\\t", \$0)}} END {printf("\\n");}' rnaFold.out > ${terminator.baseName - "_allocated_terminators"}.tsv
    gmt psconvert *.ps
    """
}

output_figures.collect().set{all_figures}
motifResult.collect().set{resultingMotifs}

/*
    Generates a HTML report using Flask by copying static files and templates needed for the report creation script
    Returns the report in the specified output path
*/
process CREATEREPORT {
    label "Report"
    publishDir "$output_path/", pattern: "*.tsv", mode: 'copy'

    input:
    val templates
    val staticfiles
    val all_figures
    val motif_done
    val output_path
	
	output:
	path "*.tsv" 

    """
    cp -r $staticfiles static
    cp -r $templates templates
    python3 $pyCreateReport --path $output_path
    """
}
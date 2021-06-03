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


/**
    Extract the transcripts from the MasterTable of TSSpredator. 
    Return the given queries for classification, the promoter regions, as
    well as a txt file with the name of the genomes and a tsv file with ignored and analyzed positions
*/
process fromMasterToFasta{
    publishDir "$output_path/queries", pattern: "*.tsv", mode: 'copy'
    publishDir "$output_path/queries", pattern: "*.fasta" , mode: 'copy'
    publishDir "$output_path", pattern: "genomes_text.txt" , mode: 'copy'

    input: 
    val table from table_ch
    val genomes from genomes_path
    val gffs from gff_path
    val output_path from output_path

    output: 

    path "genomes_text.txt" into file_text 
    path "*_queries.fasta" into queries, queries_CNIT
    path "*_promoter_regions.fasta" into promoters
    path "tss_analyzed.tsv" into summary_transcripts
    path "*.tsv" into analyzed_tss

    shell:


    if (params.conditions == "TRUE")  
    """
        python3 $pythonScriptTableToQuery $table $genomes $gffs  --conditions
    """

    else 
   """
        python3 $pythonScriptTableToQuery $table $genomes $gffs
    """

}

/**
    Uses meme to run a Motif Analysis on the promoter region
    It returns the output of meme but also a TSV file for own parsing of the data. 
*/
process motifAnalysis{
    publishDir "$output_path/MotifAnalysis/${promoter.baseName - "_promoter_regions"}", mode: 'copy'

    input: 
    each promoter from promoters
    val output_path from output_path


    output:
    path "*" into motifResult


    """
    meme $promoter -dna -nmotifs $params.motifNumber -minw 5 -maxw 20 
    python3 $pyMemeParser --meme meme_out/meme.xml --genome ${promoter.baseName - "_promoter_regions"}
    """


}

process runCNIT{
    // publishDir "$output_path/Classification/CNIT/Results", mode: 'copy'

    input: 
    each query from queries
    val output_path from output_path


    output:
    path "${query.baseName}.index" into cnit_results

    """
    python2 $cnit -f $query -o ${query.baseName} -m "pl"
    mv ${query.baseName}/${query.baseName}.index .
    """
}

cnit_results.collect().set{cnit_collected}


/*
    Evaluate output for cnit and converts it to a tsv file. 
*/
process evaluateCNIT {
    publishDir "$output_path/Classification/CNIT/Evaluation"

    input: 
    path cnit_index from cnit_collected
    val output_path from output_path

    output:
    path "*.tsv" into cnit_eval

    """
    python3 $pyFromCNITtoTSV --cnit_list $cnit_index
    """
}

file_text.map{it.getText().split("\n")}.set{ fixed_genomes } 
genomes_to_process = fixed_genomes.flatten().collect()
genomes_to_process.map{it.join(",")}.set{acc_codes_to_id}

/**
    Using the NCBI-IDs of the given genomes, 
    it extracts the corresponding TaxIDs from the Taxonomy tree of NCBI
**/
process getLCAofGenomes{ // makes one call for all genomes
    tag "$ids"
    input: 
    val ids from acc_codes_to_id

    output: 
    stdout tax_id

    shell:
    """
        efetch -db nuccore -id $ids -format docsum | xtract -pattern DocumentSummary -element TaxId 

    """
    }

tax_id.map{it.split("\n").join(" ")}.set{ tax_id_test } 

/**
    From the given TaxIds, it extracts the LCA from the species in question
**/
process getIDofLCA {
    input: 
    val ids from tax_id_test

    output: 
    path "common_species_id.txt" into common_species_id

    shell: 

    """
        python3 $pyScriptCommonSpecies $ids > common_species_id.txt
    """
}


common_species_id.readLines().map{it[0]}.set{cut_common_species_id}

/**
    Creates a TaxIDlist for blast with all corresponding species under the extracted LCA
**/
process getBlastIDs {
    input: 
    val taxidlist from cut_common_species_id

    output: 
    path "taxidlist.taxid" into taxidlist_file

    """
        get_species_taxids -t $taxidlist > taxidlist.taxid
    """
}

/**
    Blast against the nt-database using the extracted queries and the restricted TaxIDs
*/

process blastFasta{
    // publishDir "$output_path/blasted_queries", mode: 'copy'

	input:
    val output_path from output_path
    each query from queries
    path taxidlist from taxidlist_file
    env BLASTDB from "/tmp"

    output:
    path "*.blastn" into blasted_files
    shell: 

    """
    blastn -query $query -db nt -out ${query.baseName - "_queries"}.blastn  -outfmt "6 qseqid qstart qend qseq sseqid sstart send  sseq evalue bitscore pident frames qcovhsp" -task dc-megablast  -taxidlist $taxidlist

    """
}


/**
    Runs the evaluations for the BLAST-hits and extracts the corresponding hits for the pairwise alignment
**/
process evaluateBlast {
    publishDir "$output_path/filtered_blast", mode: 'copy'


    input: 
    path blasted_files from blasted_files
    val output_path from output_path



    output:
    path "*.tsv" into filtered_queries_ch, filtered_queries_for_qrna_check
    // path "evaluation*" into evaluation_ch
    // TODO: Still add the evaluation table, since there might be transcripts without any hit. 
    """
        python3 $pyEvaluateBlast $blasted_files -t 10
    """


}

/**
    runs QRNA for the classification
*/
process runQRNAnormal{
    cache "lenient"
    // publishDir "$output_path/Classification/QRNA/Results", mode: 'copy' // Put in one folder QRNA/result

    input: 
    each blasted_file from filtered_queries_ch
    val output_path from output_path
    env QRNADB from eqrnaLib


    output:
    path "${blasted_file.baseName}.qrna" into qrna_normal

    """
    gawk 'BEGIN { OFS = "\\n"} {print ">" \$1 "|" \$2 "-" \$3 "|" \$9"," \$10"," \$11 ", Frames:" \$12", Score:" \$17, \$4, ">" \$5 "|" \$6 "-" \$7, \$8 "\\n"}' $blasted_file > ${blasted_file.baseName}.fa
    $eqrna --ones ${blasted_file.baseName}.fa > ${blasted_file.baseName}.qrna
    """
}
//awk 'BEGIN { OFS = "\n"} {print ">" $1 "|" $2 "-" $3 "|" $9"," $10"," $11 ", Frames:" $12", Score:" $17, $4, ">" $5 "|" $6 "-" $7, $8 "\n"}' 


qrna_normal.collect().set{qrna_evaluate}
filtered_queries_for_qrna_check.collect().set{blasted_files_qrna_check}

/**
    Transcribes the output of QRNA to a TSV file and corrects the position of each UTR region using the blast files.  
*/
process evaluateQRNA{
    publishDir "$output_path/Classification/QRNA/Evaluation", pattern: "*.tsv", mode: 'copy'


    input:
    path path_qrna_files from qrna_evaluate
    path path_filtered_blast_files from blasted_files_qrna_check
    val output_path from output_path

    output:
    path "*.tsv" into qrna_eval

    """
        python3 $pyFromQRNAtoTSV --qrna_list $path_qrna_files --filteredBlast_list $path_filtered_blast_files
    """

}


cnit_eval.collect().set{cnit_eval_collected}
qrna_eval.collect().set{qrna_eval_collected}

/**
    Compares the output of CNIT and QRNA to call the final classification on each transcript
**/
process compareCNITandQRNA {
    publishDir "$output_path/Classification",  mode: 'copy'

    input:
    path qrna_eval_file from qrna_eval_collected
    path cnit_eval_file from cnit_eval_collected
    val output_path from output_path


    output:
    path "*final_classification.tsv" into final_eval_classification
    path "*.gff" into gff_classification
    path "*.crd" into crd_files

    shell:     
    """
        python3 $pyDecideClass --qrna $qrna_eval_file --cnit $cnit_eval_file
    """

}
/**
    Identifies possible Rho-Dependent terminators for the found genomes
*/
process rhotermpredict {
    publishDir "$output_path/Terminators/RhoTermPredict",  mode: 'copy'

    input:
    path genome from genomes_rhoterm
    val output_path from output_path

    output:
    path "*.gff" into rhoterm_gffs
    
    shell:
    """
        python3 $pyRhoTermPredict $genome --output ${genome.baseName}_rhoterm
        python3 $pyTransformRhoTermToGFF --rhoterm *${genome.baseName}_rhoterm.tsv --genome ${genome.baseName}
    """


}

/**
    Identifies Rho-Independent terminators for the found genomes 
*/
process nocornac_gffs {
    publishDir "$output_path/Terminators/nocoRNAc",  mode: 'copy'

    input:
    val config_file from nocornacConfig
    val proj_dir from projectDir
    val output_path from output_path
    path genome from genomes_nocornac


    output: 
    path "*.gff" into nocornac_gffs

    """
    sed 's@dataPath = data@dataPath =  $proj_dir/nocornac/data@g' $config_file > config_temp.conf
    sed 's@transtermPath = progs/@transtermPath =  $proj_dir/bin/progs/@g' config_temp.conf > config.conf
    java -Xmx1G -jar $nocornac  -genomeFastaFile $genome -gffOutFile ${genome.baseName}_nocornac.gff -terminators
    """
}

nocornac_gffs.collect().set{nocornac_collected}
rhoterm_gffs.collect().set{rhoterm_collected}
crd_files.collect().set{crd_files_collected}


// 
// Finding the possible terminators for each transcript

process findingTerminators {
    publishDir "$output_path/Terminators",   mode: 'copy'

    input:
    path nocornac from nocornac_collected
    path rhoterm from rhoterm_collected
    path crd from crd_files_collected
    path tsv from summary_transcripts
    val output_path from output_path


    output:
    path "*.tsv" into terminators_allocation 

    """
    python3 $pyAllocateTerminators --nocornac $nocornac --rhoterm $rhoterm --crd $crd --tssAnalyzed $tsv
    """
}

/**
Runs a MFE analysis on the transcripts classified as RNAs.
Returns also a JPG file to visualize the Sec. Structure of the transcripts
*/
process runRNAfold{
    // publishDir "$output_path/SecondaryStructure/${terminator.baseName - "_allocated_terminators"}/Transcripts", pattern: "*.fasta", mode: 'copy'
    publishDir "$output_path/SecondaryStructure/${terminator.baseName - "_allocated_terminators"}/Visualizations", pattern: "*.jpg", mode: 'copy'
    publishDir "$output_path/SecondaryStructure/${terminator.baseName - "_allocated_terminators"}/", pattern: "*.tsv", mode: 'copy'

    input:
    each terminator from terminators_allocation 
    val output_path from output_path
    val genomes from genomes_path

    output:
    path "*.jpg" 
    path "*.tsv" 
    path "*.fasta"
    stdout into output_figures

    """
    python3 $pyExtractRNATranscripts --rnas $terminator --genome_path $genomes
    $rnafold --noLP -i *.fasta > rnaFold.out
    gawk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""), \$0);N++;next;} {match(\$0, /(.*)\\s\\((.*)\\)/, ary); if (length(ary)>0) {printf("%s\\t%s",ary[1],ary[2]);} else { printf("%s\\t", \$0)}} END {printf("\\n");}' rnaFold.out > ${terminator.baseName - "_allocated_terminators"}.tsv
    gmt psconvert *.ps
    
    """

}

output_figures.collect().set{all_figures}
motifResult.collect().set{resultingMotifs}

process createReport {

    publishDir "$output_path/", pattern: "*.tsv", mode: 'copy'
    input:
    val output_path from output_path
    val templates from templates
    val staticfiles from staticfiles
    val all_figures from all_figures
    val motif_done from resultingMotifs
	
	output:
	path "*.tsv" 


    """
    cp -r $staticfiles static
    cp -r $templates templates
    python3 $pyCreateReport --path $output_path
    """

}
//    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""), $0);N++;next;} {match($0, /(.*)\s\((.*)\)/, ary); if (length(ary)>0) {printf("%s\t%s",ary[1],ary[2]);} else { printf("%s\t", $0)}} END {printf("\n");}' rnaFold.out >tabbedRNAfold.out


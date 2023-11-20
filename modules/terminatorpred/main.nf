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
    
    script:
    """
        python3 $pyRhoTermPredict $genome --output ${genome.baseName}_rhoterm
        python3 $pyTransformRhoTermToGFF --rhoterm *${genome.baseName}_rhoterm.tsv --genome ${genome.baseName}
    """
}

/*
    Identifies Rho-Independent terminators for the found genomes
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

    script:
    """
    sed 's@dataPath = data@dataPath =  $proj_dir/nocornac/data@g' $config_file > config_temp.conf
    sed 's@transtermPath = progs/@transtermPath =  $proj_dir/bin/progs/@g' config_temp.conf > config.conf
    java -Xmx1G -jar $nocornac  -genomeFastaFile $genome -gffOutFile ${genome.baseName}_nocornac.gff -terminators
    """
}

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

    script:
    """
    python3 $pyAllocateTerminators --nocornac $nocornac --rhoterm $rhoterm --crd $crd --tssAnalyzed $tsv
    """
}
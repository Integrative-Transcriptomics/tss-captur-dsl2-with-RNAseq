/*
    Identifies possible Rho-Dependent terminators for the found genomes
*/
process RHOTERMPREDICT {
    label params.terminatorPredLabel
    publishDir "$outputPath/Terminators/RhoTermPredict",  mode: params.pubDirMode

    input:
    path genome
    val outputPath

    output:
    path "*.gff", emit: gffRhoterm
    
    script:
    """
        python3 $params.rhoTermPredict $genome --output ${genome.baseName}_rhoterm
        python3 $params.pyTransformRhoTermToGFF --rhoterm *${genome.baseName}_rhoterm.tsv --genome ${genome.baseName}
    """
}

/*
    Identifies Rho-Independent terminators for the found genomes
*/
process NOCORNAC {
    label params.terminatorPredLabel
    publishDir "$outputPath/Terminators/nocoRNAc",  mode: params.pubDirMode

    input:
    val config
    val projDir
    path genome
    val outputPath

    output: 
    path "*.gff", emit: gffNocornac

    script:
    """
        sed 's@dataPath = data@dataPath = $projDir/nocornac/data@g' $config > config_temp.conf
        sed 's@transtermPath = progs/@transtermPath = $projDir/bin/progs/@g' config_temp.conf > config.conf
        java -Xmx1G -jar $params.nocornac -genomeFastaFile $genome -gffOutFile ${genome.baseName}_nocornac.gff -terminators
    """
}

/*
    Finds the possible terminators for each transcript
*/
process FINDTERMINATORS {
    label params.terminatorPredLabel
    publishDir "$outputPath/Terminators",   mode: params.pubDirMode

    input:
    path nocornac
    path rhoterm
    path crd
    path tsv
    val outputPath

    output:
    path "*.tsv", emit: terminatorsAllocation 

    script:
    """
        python3 $params.pyAllocateTerminators --nocornac $nocornac --rhoterm $rhoterm --crd $crd --tssAnalyzed $tsv
    """
}
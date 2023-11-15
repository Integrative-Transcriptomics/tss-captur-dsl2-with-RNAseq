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
eqrna = Paths.get(nameProjDir, "/bin/eqrna")
templates = Paths.get(nameProjDir, "/bin/templates")
staticfiles = Paths.get(nameProjDir, "/bin/static")
eqrnaLib = Paths.get(nameProjDir, "/bin/lib")
nocornac = Paths.get(nameProjDir, "/bin/nocornac.jar")
nocornacConfig = Paths.get(nameProjDir, "/bin/config.conf")
cnit = Paths.get(nameProjDir, "/bin/CNCI2.py")
pyRhoTermPredict = Paths.get(nameProjDir, "/bin/RhoTermPredict.py")
rnafold = Paths.get(nameProjDir, "/bin/RNAfold")
blastDB = "/tmp"

// TODO: Imports

workflow {
    // TODO: Pipeline logic
}
./nextflow run tss_captur.nf --inputTable input/CampySubset/MasterTableSubset.tsv --inputGenomes input/CampySubset/fasta/ --inputGFFs input/CampySubset/gff/ --outputDir output/InterfaceTesting/ --blastdb /ceph/ibmi/it/references/library/nt_new  -with-docker mwittep/tsscaptur  -with-dag flowchart.mmd -resume
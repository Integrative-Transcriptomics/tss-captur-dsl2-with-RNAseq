/*
    Deletes all files in the Nextflow work directory
*/
process CLEANWORKDIR {
    input:
    val report

    script:
    """
    cd ../../..
    rm -rf work/*
    """
}
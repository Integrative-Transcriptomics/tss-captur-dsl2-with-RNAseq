/*
    Deletes all files in the Nextflow work directory older than x days
*/
process CLEANWORKDIR {
    input:
    val _

    script:
    """
    find ../../ -type d -mtime +${params.cleanOlder} -exec rm -rf {} \;
    """
}
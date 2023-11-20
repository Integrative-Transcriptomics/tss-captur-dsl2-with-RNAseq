/*
    Generates a HTML report using Flask by copying static files and templates needed for the report creation script
    Returns the report in the specified output path
*/
process CREATEREPORT {
    label "Report"
    publishDir "$output_path/", pattern: "*.tsv", mode: 'copy'

    input:
    val all_figures
    val motif_done
    val output_path
	
	output:
	path "*.tsv" 

    script:
    """
        cp -r $params.staticfiles static
        cp -r $params.templates templates
        python3 $params.pyCreateReport --path $output_path
    """
}
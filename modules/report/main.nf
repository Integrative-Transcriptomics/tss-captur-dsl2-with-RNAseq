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
include { CNIT, EVALCNIT, QRNA, EVALQRNA, COMPARECNITQRNA } from '../modules/classification'

workflow CLASSIFICATION {
    take:
        queries
        filtered_queries
        eqrnaLib
        output_path

    main:
        CNIT(queries)
        EVALCNIT(CNIT.out.cnit_results.collect(), output_path)
        QRNA(filtered_queries, eqrnaLib, output_path)
        EVALQRNA(QRNA.out.qrna_normal.collect(),filtered_queries, output_path)
        COMPARECNITQRNA(EVALQRNA.out.qrna_eval.collect(), EVALCNIT.out.cnit_eval.collect(), output_path)

    emit:
        crd_files = COMPARECNITQRNA.out.crd_files
}
include { CNIT; EVALCNIT; QRNA; EVALQRNA; COMPARECNITQRNA } from '../modules/classification'

workflow CLASSIFICATION {
    take:
        queries
        blastFiltered
        outputPath

    main:
        CNIT(queries)
        EVALCNIT(CNIT.out.cnitResults.collect(), outputPath)
        QRNA(blastFiltered, params.eqrnaLib, outputPath)
        EVALQRNA(QRNA.out.qrnaNormal.collect(), blastFiltered, outputPath)
        COMPARECNITQRNA(EVALQRNA.out.qrnaEval.collect(), EVALCNIT.out.cnitEval.collect(), outputPath)

    emit:
        crdFiles = COMPARECNITQRNA.out.crdFiles
}
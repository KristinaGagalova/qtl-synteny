#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QTL_SYNTENY_FLOW } from './workflows/qtl-synteny'

workflow {
    QTL_SYNTENY_FLOW()
}

nextflow.enable.dsl = 2

include { PREPARE_EVIDENCE } from '../subworkflows/local/prepare_evidence'
include { SYNTENY_VIEWS    } from '../subworkflows/local/synteny_views'

workflow QTL_SYNTENY_FLOW {

    main:
    required = [
        source_genome  : params.source_genome,
        target_genome  : params.target_genome,
        source_gff     : params.source_gff,
        target_gff     : params.target_gff,
        source_proteins: params.source_proteins,
        qtl_bed        : params.qtl_bed
    ]
    required.each { k, v -> if (!v) error "Missing required parameter: --${k}" }
    if (!params.skip_rbh && !params.target_proteins)
        error "--target_proteins is required unless --skip_rbh"

    ch_versions = Channel.empty()
    no_file     = file("${projectDir}/assets/NO_FILE")

    meta_src = [id: params.source_name ?: file(params.source_genome).simpleName]
    meta_tgt = [id: params.target_name ?: file(params.target_genome).simpleName]

    ch_source_fa  = Channel.value([meta_src, file(params.source_genome, checkIfExists: true)])
    ch_target_fa  = Channel.value([meta_tgt, file(params.target_genome, checkIfExists: true)])
    ch_source_gff = Channel.value([meta_src, file(params.source_gff,    checkIfExists: true)])
    ch_target_gff = Channel.value([meta_tgt, file(params.target_gff,    checkIfExists: true)])
    ch_source_pep = Channel.value([meta_src, file(params.source_proteins, checkIfExists: true)])
    ch_target_pep = params.target_proteins
        ? Channel.value([meta_tgt, file(params.target_proteins, checkIfExists: true)])
        : Channel.value([meta_tgt, no_file])
    ch_qtl = Channel.value(file(params.qtl_bed, checkIfExists: true))

    ch_source_expr = params.source_expression
        ? Channel.value(file(params.source_expression, checkIfExists: true)) : Channel.value(no_file)
    ch_target_expr = params.target_expression
        ? Channel.value(file(params.target_expression, checkIfExists: true)) : Channel.value(no_file)

    PREPARE_EVIDENCE(ch_source_gff, ch_target_gff, ch_source_pep, ch_target_pep, ch_target_fa)
    ch_versions = ch_versions.mix(PREPARE_EVIDENCE.out.versions)

    SYNTENY_VIEWS(
        ch_qtl,
        ch_source_fa,
        ch_target_fa,
        PREPARE_EVIDENCE.out.target_index,
        PREPARE_EVIDENCE.out.source_genes,
        PREPARE_EVIDENCE.out.target_genes,
        PREPARE_EVIDENCE.out.miniprot,
        PREPARE_EVIDENCE.out.rbh,
        ch_source_expr,
        ch_target_expr
    )
    ch_versions = ch_versions.mix(SYNTENY_VIEWS.out.versions)

    emit:
    views    = SYNTENY_VIEWS.out.views
    index    = SYNTENY_VIEWS.out.index
    regions  = SYNTENY_VIEWS.out.regions
    versions = ch_versions
}

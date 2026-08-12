/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PREPARE_EVIDENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Everything that is genome-wide and therefore built ONCE, then reused by
    every QTL:

      * gene tables for both annotations
      * minimap2 index of the target             -> STEP 1 (DNA synteny)
      * miniprot placement of source proteins    -> STEP 3 (homology, no
        target annotation needed)
      * MMseqs2 easy-rbh between the proteomes   -> STEP 4 (orthology, needs
        both annotations)

    Steps 3 and 4 are independent of step 1 here; they are only intersected
    with the syntenic intervals later, per QTL, in ANNOTATE_REGIONS.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GFF_TO_GENE_BED as GENES_SOURCE } from '../../modules/local/gff_to_gene_bed'
include { GFF_TO_GENE_BED as GENES_TARGET } from '../../modules/local/gff_to_gene_bed'
include { MINIMAP2_INDEX                  } from '../../modules/local/minimap2_index'
include { MINIPROT_INDEX                  } from '../../modules/local/miniprot_index'
include { MINIPROT_ALIGN                  } from '../../modules/local/miniprot_align'
include { MMSEQS_RBH                      } from '../../modules/local/mmseqs_rbh'

workflow PREPARE_EVIDENCE {

    take:
    ch_source_gff
    ch_target_gff
    ch_source_pep
    ch_target_pep
    ch_target_fa

    main:
    ch_versions = Channel.empty()
    no_file     = file("${projectDir}/assets/NO_FILE")

    GENES_SOURCE(ch_source_gff)
    GENES_TARGET(ch_target_gff)
    ch_versions = ch_versions.mix(GENES_SOURCE.out.versions)
    ch_versions = ch_versions.mix(GENES_TARGET.out.versions)

    // STEP 1 support: one target index for every per-QTL alignment
    MINIMAP2_INDEX(ch_target_fa)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    // STEP 3: source proteins spliced onto the target genome
    if (!params.skip_miniprot) {
        MINIPROT_INDEX(ch_target_fa)
        MINIPROT_ALIGN(ch_source_pep, MINIPROT_INDEX.out.index)
        ch_versions = ch_versions.mix(MINIPROT_INDEX.out.versions)
        ch_versions = ch_versions.mix(MINIPROT_ALIGN.out.versions)
        ch_miniprot = MINIPROT_ALIGN.out.gff.map { meta, f -> f }
    } else {
        ch_miniprot = Channel.value(no_file)
    }

    // STEP 4: reciprocal best hits between the two annotations
    if (!params.skip_rbh) {
        MMSEQS_RBH(ch_source_pep, ch_target_pep)
        ch_versions = ch_versions.mix(MMSEQS_RBH.out.versions)
        ch_rbh = MMSEQS_RBH.out.rbh.map { meta, f -> f }
    } else {
        ch_rbh = Channel.value(no_file)
    }

    emit:
    source_genes = GENES_SOURCE.out.genes.map { meta, f -> f }
    target_genes = GENES_TARGET.out.genes.map { meta, f -> f }
    target_index = MINIMAP2_INDEX.out.index
    miniprot     = ch_miniprot
    rbh          = ch_rbh
    versions     = ch_versions
}

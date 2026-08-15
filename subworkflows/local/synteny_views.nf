/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SYNTENY_VIEWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Scattered per QTL:

      STEP 1  slice the QTL interval out of the source assembly and align it
              against the WHOLE target with minimap2
      STEP 2  rank the target sequences it hit by aligned bp, keep top N -
              these are the syntenic regions, defined by DNA alone
      STEP 3  keep the miniprot placements of this QTL's source proteins that
              fall INSIDE those regions
      STEP 4  keep the reciprocal best hits whose target gene falls INSIDE
              those regions
      then    build the interactive view. Genes displayed (and used for the
              expression panel) are all source genes in the QTL plus all
              target genes in the regions - including those with no protein
              link at all.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACT_SOURCE_REGION } from '../../modules/local/extract_source_region'
include { MINIMAP2_REGION       } from '../../modules/local/minimap2_region'
include { RANK_DNA_REGIONS      } from '../../modules/local/rank_dna_regions'
include { ANNOTATE_REGIONS      } from '../../modules/local/annotate_regions'
include { BUILD_SYNTENY_HTML    } from '../../modules/local/build_synteny_html'
include { BUILD_INDEX_HTML      } from '../../modules/local/build_index_html'
include { EXTRACT_SELECTED_SCAFFOLDS } from '../../modules/local/extract_selected_scaffolds'

workflow SYNTENY_VIEWS {

    take:
    ch_qtl
    ch_source_fa
    ch_target_fa
    ch_target_index
    ch_source_genes
    ch_target_genes
    ch_miniprot
    ch_rbh
    ch_source_expr
    ch_target_expr
    ch_source_ann
    ch_target_ann

    main:
    ch_versions = Channel.empty()

    ch_qtl_rows = ch_qtl
        .splitCsv(sep: '\t', strip: true)
        .filter { row -> row.size() >= 3 && !row[0].startsWith('#') &&
                         !row[0].startsWith('track') && !row[0].startsWith('browser') }
        .map { row ->
            qid  = row.size() > 3 && row[3] ? row[3] : "${row[0]}:${row[1]}-${row[2]}"
            safe = qid.replaceAll(/[^A-Za-z0-9._-]/, '_')
            [id: safe, qtl_id: qid, chrom: row[0],
             start: row[1] as long, end: row[2] as long]
        }

    // STEP 1
    EXTRACT_SOURCE_REGION(ch_qtl_rows, ch_source_fa)
    MINIMAP2_REGION(EXTRACT_SOURCE_REGION.out.fasta, ch_target_index)
    ch_versions = ch_versions.mix(EXTRACT_SOURCE_REGION.out.versions.first())
    ch_versions = ch_versions.mix(MINIMAP2_REGION.out.versions.first())

    // STEP 2
    RANK_DNA_REGIONS(MINIMAP2_REGION.out.paf)
    ch_versions = ch_versions.mix(RANK_DNA_REGIONS.out.versions.first())

    // STEPS 3 + 4: intersect the genome-wide protein evidence with the
    // regions. The PAF is passed through too, so ANNOTATE_REGIONS can split
    // source coverage into gene-bearing vs gene-less scaffolds.
    ch_ann_in = RANK_DNA_REGIONS.out.regions
        .join(MINIMAP2_REGION.out.paf)
        .combine(ch_source_genes)
        .combine(ch_target_genes)
        .combine(ch_miniprot)
        .combine(ch_rbh)

    ANNOTATE_REGIONS(ch_ann_in)
    ch_versions = ch_versions.mix(ANNOTATE_REGIONS.out.versions.first())

    // interactive view per QTL
    ch_html_in = ANNOTATE_REGIONS.out.regions
        .join(ANNOTATE_REGIONS.out.links)
        .join(ANNOTATE_REGIONS.out.genes)
        .join(MINIMAP2_REGION.out.paf)
        .join(ANNOTATE_REGIONS.out.coverage)

    BUILD_SYNTENY_HTML(ch_html_in, ch_source_expr, ch_target_expr, ch_source_ann, ch_target_ann)
    ch_versions = ch_versions.mix(BUILD_SYNTENY_HTML.out.versions.first())

    // FASTA of the selected target sequences, for downstream work
    EXTRACT_SELECTED_SCAFFOLDS(ANNOTATE_REGIONS.out.regions, ch_target_fa)
    ch_versions = ch_versions.mix(EXTRACT_SELECTED_SCAFFOLDS.out.versions.first())

    BUILD_INDEX_HTML(
        ANNOTATE_REGIONS.out.regions.map { meta, f -> f }.collectFile(
            name: 'all_regions.tsv', keepHeader: true, skip: 1),
        ch_qtl,
        BUILD_SYNTENY_HTML.out.html.map { meta, f -> f }.collect()
    )
    ch_versions = ch_versions.mix(BUILD_INDEX_HTML.out.versions)

    emit:
    regions  = ANNOTATE_REGIONS.out.regions
    links    = ANNOTATE_REGIONS.out.links
    views     = BUILD_SYNTENY_HTML.out.html
    scaffolds = EXTRACT_SELECTED_SCAFFOLDS.out.full
    windows   = EXTRACT_SELECTED_SCAFFOLDS.out.windows
    index    = BUILD_INDEX_HTML.out.html
    versions = ch_versions
}

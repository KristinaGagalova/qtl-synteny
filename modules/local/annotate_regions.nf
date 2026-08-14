process ANNOTATE_REGIONS {
    tag "$meta.qtl_id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(regions), path(paf), path(src_genes), path(tgt_genes), path(miniprot), path(rbh)

    output:
    tuple val(meta), path("*.regions.tsv")  , emit: regions
    tuple val(meta), path("*.links.tsv")    , emit: links
    tuple val(meta), path("*.genes.tsv")    , emit: genes
    tuple val(meta), path("*.coverage.tsv") , emit: coverage
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def mp_arg  = miniprot.name == 'NO_FILE' ? '' : "--miniprot ${miniprot}"
    def rbh_arg = rbh.name      == 'NO_FILE' ? '' : "--rbh ${rbh}"
    """
    annotate_regions.py \\
        --regions ${regions} \\
        --paf ${paf} \\
        --qtl-id '${meta.qtl_id}' \\
        --qtl-chrom '${meta.chrom}' \\
        --qtl-start ${meta.start} \\
        --qtl-end ${meta.end} \\
        --source-genes ${src_genes} \\
        --target-genes ${tgt_genes} \\
        ${mp_arg} \\
        ${rbh_arg} \\
        --out-links ${meta.id}.links.tsv \\
        --out-genes ${meta.id}.genes.tsv \\
        --out-regions ${meta.id}.regions.tsv \\
        --out-coverage ${meta.id}.coverage.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.regions.tsv ${meta.id}.links.tsv ${meta.id}.genes.tsv ${meta.id}.coverage.tsv versions.yml
    """
}

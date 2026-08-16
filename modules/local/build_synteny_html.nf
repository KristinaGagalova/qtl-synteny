process BUILD_SYNTENY_HTML {
    tag "$meta.qtl_id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(regions), path(links), path(genes), path(paf), path(coverage)
    path source_expr
    path target_expr
    path source_annotation
    path target_annotation

    output:
    tuple val(meta), path("*.html")             , emit: html
    tuple val(meta), path("*.match_report.txt") , emit: match_report
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def se_arg = source_expr.name == 'NO_FILE' ? '' : "--source-expr ${source_expr}"
    def te_arg = target_expr.name == 'NO_FILE' ? '' : "--target-expr ${target_expr}"
    def sa_arg = source_annotation.name == 'NO_FILE' ? '' : "--source-annotation ${source_annotation}"
    def ta_arg = target_annotation.name == 'NO_FILE' ? '' : "--target-annotation ${target_annotation}"
    """
    build_synteny_html.py \\
        --qtl-id '${meta.qtl_id}' \\
        --qtl-chrom '${meta.chrom}' \\
        --qtl-start ${meta.start} \\
        --qtl-end ${meta.end} \\
        --regions ${regions} \\
        --links ${links} \\
        --genes ${genes} \\
        --paf ${paf} \\
        --coverage ${coverage} \\
        ${se_arg} \\
        ${te_arg} \\
        ${sa_arg} \\
        ${ta_arg} \\
        --match-report ${meta.id}.match_report.txt \\
        --out ${meta.id}.html \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.html ${meta.id}.match_report.txt versions.yml
    """
}

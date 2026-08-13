process EXTRACT_SELECTED_SCAFFOLDS {
    tag "$meta.qtl_id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(regions)
    tuple val(meta2), path(target_fa)

    output:
    tuple val(meta), path("*.scaffolds.fa")       , emit: full
    tuple val(meta), path("*.windows.fa")         , emit: windows
    tuple val(meta), path("*.selected_scaffolds.tsv"), emit: list
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    extract_selected_scaffolds.py \\
        --regions ${regions} \\
        --qtl-id '${meta.qtl_id}' \\
        --target-fasta ${target_fa} \\
        --out-full ${meta.id}.scaffolds.fa \\
        --out-window ${meta.id}.windows.fa \\
        --out-list ${meta.id}.selected_scaffolds.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.scaffolds.fa ${meta.id}.windows.fa ${meta.id}.selected_scaffolds.tsv versions.yml
    """
}

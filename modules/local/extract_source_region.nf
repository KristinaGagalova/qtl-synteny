process EXTRACT_SOURCE_REGION {
    tag "$meta.qtl_id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    val meta
    tuple val(meta2), path(source_fa)

    output:
    tuple val(meta), path("*.src.fa"), emit: fasta
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    extract_source_region.py \\
        --source-fasta ${source_fa} \\
        --chrom '${meta.chrom}' \\
        --start ${meta.start} \\
        --end ${meta.end} \\
        --out ${meta.id}.src.fa \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.src.fa versions.yml
    """
}

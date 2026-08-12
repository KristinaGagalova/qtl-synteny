process RANK_DNA_REGIONS {
    tag "$meta.qtl_id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(paf)

    output:
    tuple val(meta), path("*.dna_regions.tsv"), emit: regions
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rank_dna_regions.py \\
        --paf ${paf} \\
        --qtl-id '${meta.qtl_id}' \\
        --src-chrom '${meta.chrom}' \\
        --src-start ${meta.start} \\
        --src-end ${meta.end} \\
        --out ${meta.id}.dna_regions.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.dna_regions.tsv versions.yml
    """
}

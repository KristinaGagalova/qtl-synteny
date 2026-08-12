process MINIMAP2_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::minimap2=2.28"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0' :
        'biocontainers/minimap2:2.28--he4a0461_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '-x asm20'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Built once and reused by every QTL, so no QTL re-indexes the target.
    minimap2 ${args} -t ${task.cpus} -d ${prefix}.mmi ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mmi versions.yml
    """
}

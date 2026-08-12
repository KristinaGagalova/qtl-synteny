process MINIPROT_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::miniprot=0.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/miniprot:0.13--h577a1d6_0' :
        'biocontainers/miniprot:0.13--h577a1d6_0' }"

    input:
    tuple val(meta) , path(pep)    // source proteins
    tuple val(meta2), path(index)  // miniprot index of the target genome

    output:
    tuple val(meta), path("*.miniprot.gff"), emit: gff
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    miniprot \\
        -t ${task.cpus} \\
        --gff \\
        ${args} \\
        ${index} \\
        ${pep} \\
        > ${prefix}.miniprot.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miniprot: \$(miniprot --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.miniprot.gff versions.yml
    """
}

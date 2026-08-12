process MINIMAP2_REGION {
    tag "$meta.qtl_id"
    label 'process_medium'

    conda "bioconda::minimap2=2.28"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0' :
        'biocontainers/minimap2:2.28--he4a0461_0' }"

    input:
    tuple val(meta), path(src_fa)
    tuple val(meta2), path(target_index)

    output:
    tuple val(meta), path("*.region.paf"), emit: paf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-x asm20'
    """
    # STEP 1: the QTL interval (query) against the WHOLE target assembly.
    # This is what discovers the syntenic regions - it does not depend on any
    # gene or ortholog evidence, so it still works when RBH comes back empty.
    if [ -s ${src_fa} ]; then
        minimap2 -c ${args} -t ${task.cpus} ${target_index} ${src_fa} > ${meta.id}.region.paf
    else
        : > ${meta.id}.region.paf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.region.paf versions.yml
    """
}

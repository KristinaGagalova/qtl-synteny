process MMSEQS_RBH {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mmseqs2=15.6f452"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_0' :
        'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_0' }"

    input:
    tuple val(meta) , path(source_pep)
    tuple val(meta2), path(target_pep)

    output:
    tuple val(meta), path("*.rbh.m8") , emit: rbh
    tuple val(meta), path("*.hits.m8"), emit: hits, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // easy-rbh returns reciprocal best hits directly, so the second reverse
    // search is no longer needed. The optional easy-search adds the broader
    // (non-RBH) hit set, which matters on a fragmented target where a real
    // orthologue can lose reciprocity to a tandem duplicate.
    def do_search = params.mmseqs_also_search
        ? "mkdir -p tmp_search && mmseqs easy-search ${source_pep} ${target_pep} ${prefix}.hits.m8 tmp_search --threads ${task.cpus} ${args2}"
        : "true"
    """
    mkdir -p tmp_rbh

    mmseqs easy-rbh \\
        ${source_pep} ${target_pep} ${prefix}.rbh.m8 tmp_rbh \\
        --threads ${task.cpus} ${args}
    ${do_search}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extra  = params.mmseqs_also_search ? "touch ${prefix}.hits.m8" : ""
    """
    touch ${prefix}.rbh.m8
    ${extra}
    touch versions.yml
    """
}

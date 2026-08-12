process GFF_TO_GENE_BED {
    tag "$meta.id"
    label 'process_low'
    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"
    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.genes.tsv"), emit: genes
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gff_to_gene_bed.py --gff ${gff} --out ${prefix}.genes.tsv ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'seqid\tstart\tend\tgene_id\ttranscript_id\tstrand\n' > ${prefix}.genes.tsv
    touch versions.yml
    """
}

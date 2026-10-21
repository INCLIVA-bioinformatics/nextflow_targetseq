// modules/incliva/picard/markduplicates/main.nf
process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/picard:2.18.6--py36_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"),              emit: bam
    tuple val(meta), path("${prefix}.duplicates_metrics"), emit: metrics
    path "versions.yml",                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_dedup"
    """
    picard MarkDuplicates \\
        I=${bam} \\
        O=${prefix}.bam \\
        M=${prefix}.duplicates_metrics \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_dedup"
    """
    touch ${prefix}.bam
    touch ${prefix}.duplicates_metrics
    touch versions.yml
    """
}

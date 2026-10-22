// modules/incliva/picard/collectalignmentsummarymetrics/main.nf
process PICARD_COLLECTALIGNMENTSUMMARYMETRICS {
    tag "$meta.id"
    label 'process_low'

    // Misma versión que rutina.sh (picard v2.18.6)
    container 'quay.io/biocontainers/picard:2.18.6--py36_0'

    input:
    tuple val(meta), path(bam)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("${prefix}.alignment_metrics"), emit: metrics
    path "versions.yml",                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    picard CollectAlignmentSummaryMetrics \\
        I=${bam} \\
        O=${prefix}.alignment_metrics \\
        R=${fasta} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectAlignmentSummaryMetrics --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment_metrics
    touch versions.yml
    """
}

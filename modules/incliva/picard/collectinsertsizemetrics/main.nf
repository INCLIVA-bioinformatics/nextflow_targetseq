// modules/incliva/picard/collectinsertsizemetrics/main.nf
process PICARD_COLLECTINSERTSIZEMETRICS {
    tag "$meta.id"
    label 'process_low'

    // Misma versión que rutina.sh (picard v2.18.6)
    container 'quay.io/biocontainers/picard:2.18.6--py36_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.insert_size_metrics"), emit: metrics
    path "versions.yml",                                   emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // El contenedor no incluye R: picard escribe las métricas y después falla al dibujar el
    // histograma (H=). rutina.sh sólo utiliza las métricas, así que se acepta si el fichero existe.
    """
    picard CollectInsertSizeMetrics \\
        I=${bam} \\
        O=${prefix}.insert_size_metrics \\
        H=${prefix}.insert_size_histogram.pdf \\
        ${args} \\
        || test -s ${prefix}.insert_size_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectInsertSizeMetrics --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.insert_size_metrics
    touch versions.yml
    """
}

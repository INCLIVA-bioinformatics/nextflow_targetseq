// modules/incliva/bamclipper/main.nf
// Soft-clipping de las secuencias de primers (paneles de amplicones), igual que rutina.sh.
// Usa el bamclipper.sh del contenedor (v1.0.0): mismo resultado que la v1.1.1 de rutina.sh.
process BAMCLIPPER {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/bamclipper:1.0.0--hdfd78af_2'

    input:
    tuple val(meta), path(bam), path(bai)
    path bedpe

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml",                    emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_primerclipped"
    // bamclipper.sh v1.0.0 no expone su versión: se toma del tag del contenedor
    def VERSION = '1.0.0'
    """
    bamclipper.sh \\
        -b ${bam} \\
        -p ${bedpe} \\
        -n ${task.cpus} \\
        ${args}

    mv ${bam.baseName}.primerclipped.bam ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamclipper: ${VERSION}
        samtools: \$(samtools --version-only)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_primerclipped"
    """
    touch ${prefix}.bam
    touch versions.yml
    """
}

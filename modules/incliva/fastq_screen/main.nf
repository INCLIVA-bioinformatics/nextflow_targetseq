// modules/incliva/fastq_screen/main.nf
process FASTQ_SCREEN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/fastq-screen:0.15.3'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_screen.txt"),  emit: txt
    tuple val(meta), path("*_screen.html"), emit: html
    path "versions.yml",                    emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    fastq_screen \\
        --outdir . \\
        --threads ${task.cpus} \\
        ${args} \\
        ${reads[0]} ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq_screen: \$(fastq_screen --version 2>&1 | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch fake_screen.txt
    touch fake_screen.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq_screen: \$(fastq_screen --version 2>&1 | sed 's/^.*version //')
    END_VERSIONS
    """
}
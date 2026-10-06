// modules/incliva/prinseq/main.nf
process PRINSEQ {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    //container 'quay.io/biocontainers/prinseq:0.20.4'
    container 'quay.io/biocontainers/prinseq-plus-plus:1.2.4--h077b44d_8'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_good_{1,2}.fastq"), emit: good_reads
    tuple val(meta), path("*_bad_{1,2}.fastq"),  emit: bad_reads, optional: true
    tuple val(meta), path("*_singletons_{1,2}.fastq"), emit: singletons, optional: true
    tuple val(meta), path("*.log"),              emit: log
    path "versions.yml",                         emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    prinseq++ \\
        -fastq ${reads[0]} \\
        -fastq2 ${reads[1]} \\
        ${args} \\
        -out_good ${prefix}_good_1.fastq \\
        -out_good2 ${prefix}_good_2.fastq \\
        -out_bad ${prefix}_bad_1.fastq \\
        -out_bad2 ${prefix}_bad_2.fastq \\
        -out_single ${prefix}_singletons_1.fastq \\
        -out_single2 ${prefix}_singletons_2.fastq \\
        2> ${prefix}_prinseq.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prinseq: \$(prinseq++ --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_good_1.fastq
    touch ${prefix}_good_2.fastq
    touch ${prefix}_prinseq.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prinseq: \$(prinseq++ --version 2>&1 | head -1)
    END_VERSIONS
    """
}
// modules/incliva/cutadapt/main.nf
process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/cutadapt:4.9--py311haab0aaa_3'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log"),      emit: log
    path "versions.yml",                 emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cutadapt \\
        -j ${task.cpus} \\
        ${args} \\
        -o ${prefix}_R1.fastq.gz \\
        -p ${prefix}_R2.fastq.gz \\
        ${reads[0]} ${reads[1]} \\
        > ${prefix}_cutadapt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1.fastq.gz
    touch ${prefix}_R2.fastq.gz
    touch ${prefix}_cutadapt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
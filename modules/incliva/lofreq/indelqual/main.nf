// modules/incliva/lofreq/indelqual/main.nf
process LOFREQ_INDELQUAL {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/lofreq:2.1.5--py311h6385618_16'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml",                 emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    lofreq indelqual \\
        -u 20,20 \\
        -f ${fasta} \\
        -o ${prefix}_indelqual.bam \\
        ${bam}

    mv ${prefix}_indelqual.bam ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq 2>&1 | grep -i version | head -n1 | sed 's/.*version //')
    END_VERSIONS
    """
}
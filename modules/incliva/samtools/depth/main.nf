// modules/incliva/samtools/depth/main.nf
process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.depth.tsv"), emit: tsv
    path "versions.yml",                 emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools depth \\
        -d 9898989898 \\
        -a \\
        -b ${bed} \\
        ${bam} \\
        > ${prefix}.depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools 2>&1 | grep -i version | head -n1 | sed 's/.*version //')
    END_VERSIONS
    """
    
}
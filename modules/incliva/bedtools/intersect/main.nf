// modules/incliva/bedtools/intersect/main.nf
process BEDTOOLS_INTERSECT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    //container 'quay.io/biocontainers/bedtools:2.31.1'
    
    def extra_args = ''

    input:
    tuple val(meta), path(a_file)  // archivo A (BED o BAM)
    tuple val(meta), path(b_file)  // archivo B (BED)
    //val extra_args = ''            // argumentos adicionales opcionales

    output:
    tuple val(meta), path("*.intersected.bed"), emit: bed
    tuple val(meta), path("*.intersected.bam"), emit: bam, optional: true
    path "*.versions", emit: versions

    script:
    def prefix = meta.id ?: a_file.baseName
    def is_bam = a_file.toString().endsWith('.bam')
    def output_file = is_bam ? "${prefix}.intersected.bam" : "${prefix}.intersected.bed"

    """
    bedtools intersect \\
        -a ${a_file} \\
        -b ${b_file} \\
        ${extra_args} \\
        > ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools //')
    END_VERSIONS
    """
}
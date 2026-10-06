// modules/incliva/scripts/get_covered_bed_regions.nf
process GET_COVERED_BED_REGIONS {
    tag "$meta.id"
    label 'process_light'

    conda 'bioconda::bedtools'
    container 'quay.io/biocontainers/bedtools:2.31.1'

    input:
    tuple val(meta), path(depth_tsv)
    val cov_threshold

    output:
    tuple val(meta), path("*.covered_regions.bed"), emit: bed
    path "*.versions", emit: versions

    script:
    """
    awk -v cov=${cov_threshold} '\$3 >= cov {print \$1"\\t"\$2-1"\\t"\$2}' \\
        ${depth_tsv} | \\
        bedtools merge -i - > \\
        ${meta.id}.covered_regions.bed
    """
}
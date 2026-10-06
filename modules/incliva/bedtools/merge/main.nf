// modules/incliva/bedtools/merge/main.nf
process BEDTOOLS_MERGE {
    tag "$meta.id"
    label 'process_light'

    conda "${moduleDir}/environment.yml"
    //container 'quay.io/biocontainers/bedtools:2.31.1'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.merged.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    """
    bedtools merge -i ${bed} > ${meta.id}.merged.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools 2>&1 | grep -i version | head -n1 | sed 's/.*version //')
    END_VERSIONS
    """
}
// modules/incliva/bedtools/sort/main.nf
process BEDTOOLS_SORT {
    tag "$meta.id"
    label 'process_light'

    conda "${moduleDir}/environment.yml"
    //container 'quay.io/biocontainers/bedtools:2.31.1'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.sorted.bed"), emit: sorted
    path "versions.yml", emit: versions

    script:
    //    sort -k1,1 -k2,2n ${bed} > ${meta.id}.sorted.bed

    """
    bedtools sort -i ${bed} > ${meta.id}.sorted.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | head -n1 | sed 's/sort (GNU coreutils) //')
    END_VERSIONS
    """
}
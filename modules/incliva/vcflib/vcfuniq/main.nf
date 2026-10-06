// modules/incliva/vcflib/vcfuniq/main.nf
process VCFLIB_VCFUNIQ {
    tag "$meta.id"
    label 'process_light'

    //conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/vcflib:1.0.15--h3fa9d83_1'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.uniq.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    vcfuniq <(zcat ${vcf}) | bgzip -c > ${meta.id}.uniq.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: \$(vcfuniq --version 2>&1 | head -1)
    END_VERSIONS
    """
}
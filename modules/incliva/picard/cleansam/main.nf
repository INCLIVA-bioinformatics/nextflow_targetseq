// modules/incliva/picard/cleansam/main.nf
process PICARD_CLEANSAM {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/picard:2.18.6--py36_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path "versions.yml",                 emit: versions


    script:
    """
    picard CleanSam \\
        I=${bam} \\
        O=${meta.id}.cleansam.tmp.bam\\
        VALIDATION_STRINGENCY=LENIENT
    
    mv ${meta.id}.cleansam.tmp.bam ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard --version)
    END_VERSIONS
    """


}
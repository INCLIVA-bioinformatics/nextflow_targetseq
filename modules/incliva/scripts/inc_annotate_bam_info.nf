// modules/incliva/scripts/inc_annotate_bam_info.nf
process INC_ANNOTATE_BAM_INFO {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::python pysam'
    container 'python:3.11'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta), path(bam)
    val annotation

    output:
    tuple val(meta), path("*.annotated.vcf"), emit: vcf
    path "*.versions", emit: versions

    script:
    script_path = file("${projectDir}/bin/inc_annotate_bam_info.py")

    """
    python ${script_path} \\
        -i ${vcf} \\
        -o ${meta.id}.annotated.vcf \\
        -b ${bam} \\
        -a ${annotation}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
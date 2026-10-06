// modules/incliva/scripts/inc_annotate_format_lofreq.nf
process INC_ANNOTATE_FORMAT_LOFREQ {
    tag "$meta.id"
    label 'process_light'

    conda 'bioconda::python'
    container 'python:3.11'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.formatted.vcf"), emit: vcf
    path "*.versions", emit: versions

    script:
    script_path = file("${projectDir}/bin/inc_annotate_format_lofreq.py")

    """
    python ${script_path} \\
        -i ${vcf} \\
        -o ${meta.id}.formatted.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
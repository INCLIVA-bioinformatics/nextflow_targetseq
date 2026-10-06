// modules/incliva/scripts/inc_select_variants.nf
process INC_SELECT_VARIANTS {
    tag "$meta.id"
    label 'process_light'

    conda 'bioconda::python'
    container 'python:3.11'

    input:
    tuple val(meta), path(vcf), path(tbi)
    val min_af
    val min_depth
    val caller

    output:
    tuple val(meta), path("*.selected.vcf"), emit: vcf
    path "*.versions", emit: versions

    script:
    script_path = file("${projectDir}/bin/inc_select_variants.py")

    """
    python ${script_path} \\
        -i ${vcf} \\
        -o ${meta.id}.selected.vcf \\
        -c ${caller} \\
        --vaf ${min_af} \\
        --depth ${min_depth}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
// modules/incliva/scripts/inc_tag_reliable.nf
process INC_TAG_RELIABLE {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::python bedtools'
    container 'python:3.11'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta), path(primers_bed)
    val tag

    output:
    tuple val(meta), path("*.tagged.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tagged.vcf.gz.tbi"), emit: tbi
    path "*.versions", emit: versions

    script:
    script_path = file("${projectDir}/bin/inc_tag_reliable.py")

    """
    python ${script_path} \\
        -i ${vcf} \\
        -o ${meta.id}.tagged.vcf \\
        -p ${primers_bed} \\
        -t ${tag}

    bgzip -c ${meta.id}.tagged.vcf > ${meta.id}.tagged.vcf.gz
    tabix -p vcf ${meta.id}.tagged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
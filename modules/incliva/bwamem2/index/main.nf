// modules/incliva/bwamem2/index/main.nf
process BWAMEM2_INDEX {
    tag "bwa_mem2_index"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9a/9ac054213e67b3c9308e409b459080bbe438f8fd6c646c351bc42887f35a42e7/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:e1f420694f8e42bd' }"
        
    input:
    tuple val(meta), path(index_files)  // lista de archivos del índice

    output:
    tuple val(meta), path("bwamem2"), emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p bwamem2
    cp ${index_files} bwamem2/
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: pre-generated index copied
    END_VERSIONS
    """
}
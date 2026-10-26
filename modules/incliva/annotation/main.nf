// ==========================================
// Anotación de rutina.sh (función annotate_variants())
// ==========================================
// Mismas herramientas y versiones que rutina.sh (VEP 110, SnpSift 4.3t). Los scripts de rutina están en bin/
// (copiados de scripts_rutina_stats/; sólo inc_samples_contamination cambia: el nombre de la corrida se pasa por
// argumento, porque en rutina se deduce de la ruta /media/scratch1/<corrida>/...).
// Los scripts de python2 se ejecutan con environment_python2.yml y los de python3/R con environment.yml.


// ==========================================
// PROCESS: INC_VEP
// ==========================================
// VEP 110 con los mismos argumentos que rutina.sh (caché merged offline + plugin SpliceAI). vep_cache es el
// directorio que rutina monta en /opt/vep/.vep (homo_sapiens_merged/110_GRCh38 y Plugins/).
process INC_VEP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ensembl-vep=110.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:110.1--pl5321h2a3209d_0' :
        'biocontainers/ensembl-vep:110.1--pl5321h2a3209d_0' }"

    input:
    tuple val(meta), path(vcf)
    path fasta
    path vep_cache
    tuple path(spliceai_snv), path(spliceai_snv_tbi)
    tuple path(spliceai_indel), path(spliceai_indel_tbi)

    output:
    tuple val(meta), path("${prefix}_annotated_vep.vcf"), emit: vcf
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    vep \\
        -i ${vcf} \\
        --use_transcript_ref \\
        --force_overwrite \\
        --offline \\
        --cache \\
        --merged \\
        --dir ${vep_cache} \\
        --dir_plugins ${vep_cache}/Plugins \\
        --fasta ${fasta} \\
        --sift b \\
        --polyphen b \\
        --gene_phenotype \\
        --numbers \\
        --vcf_info_field ANN \\
        --terms SO \\
        --hgvs \\
        --shift_hgvs 1 \\
        --canonical \\
        ${args} \\
        --biotype \\
        --xref_refseq \\
        --max_af \\
        --af_gnomad \\
        --pubmed \\
        --minimal \\
        --vcf \\
        --af \\
        --symbol \\
        --domains \\
        --fork ${task.cpus} \\
        --buffer_size 200 \\
        --plugin SpliceAI,snv=${spliceai_snv},indel=${spliceai_indel},cutoff=${params.spliceai_cutoff} \\
        -o ${prefix}_annotated_vep.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: INC_ANNOTATE_HOTSPOTS
// ==========================================
// bcftools annotate con cancerhotspots (HS en INFO)
process INC_ANNOTATE_HOTSPOTS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data' :
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(vcf)
    tuple path(hotspots), path(hotspots_tbi)

    output:
    tuple val(meta), path("${prefix}_annotated_uniq_h.vcf"), emit: vcf
    path "versions.yml"                                    , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '##INFO=<ID=HS,Number=1,Type=String,Description="Tag variants in cancer hotspots cancerhotspot.org">' > hs_header.txt

    bcftools annotate \\
        -a ${hotspots} \\
        -c CHROM,FROM,TO,HS \\
        -h hs_header.txt \\
        ${vcf} \\
        > ${prefix}_annotated_uniq_h.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: INC_SNPSIFT_ANNOTATE
// ==========================================
// SnpSift annotate (COSMIC, mixed_db, INCLIVA somatic y OncoKB). Falla si SnpSift lanza una excepción (sale con 0). ext.args: opciones de SnpSift (p.ej. -info ...);
// ext.suffix: sufijo del vcf de salida.
process INC_SNPSIFT_ANNOTATE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::snpsift=4.3.1t"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:4.3.1t--hdfd78af_3' :
        'biocontainers/snpsift:4.3.1t--hdfd78af_3' }"

    input:
    tuple val(meta), path(vcf)
    tuple path(database), path(database_index)

    output:
    tuple val(meta), path("${prefix}_${suffix}.vcf"), emit: vcf
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: 'snpsift'
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ''
    """
    SnpSift ${avail_mem} \\
        annotate \\
        ${args} \\
        ${database} \\
        ${vcf} \\
        > ${prefix}_${suffix}.vcf \\
        2> snpsift.log
    cat snpsift.log >&2

    # SnpSift sale con código 0 aunque se caiga a mitad (p.ej. base de datos truncada: "Invalid file pointer"),
    # dejando un vcf incompleto: si hay una excepción de Java en el log, el proceso falla
    if grep -q "Exception" snpsift.log; then
        echo "ERROR: SnpSift annotate ha fallado con ${database} (vcf de salida incompleto)" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: INC_MV_COSMIC_ID
// ==========================================
// inc_mv_cosmic_id_to_annotation_v1.2.py (python2): pasa los IDs de COSMIC de la columna ID a Existing_variation de ANN
process INC_MV_COSMIC_ID {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment_python2.yml"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}_annotated_uniq_h_c_cok.vcf"), emit: vcf
    path "versions.yml"                                          , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python2 ${projectDir}/bin/inc_mv_cosmic_id_to_annotation_v1.2.py \\
        -i ${vcf} \\
        -o ${prefix}_annotated_uniq_h_c_cok.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: INC_ADD_CONSERVATION
// ==========================================
// inc_add_conservation_v1.0.py (python2): Conservation en INFO a partir de conservation.bed del panel
process INC_ADD_CONSERVATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment_python2.yml"

    input:
    tuple val(meta), path(vcf)
    path conservation_bed

    output:
    tuple val(meta), path("${prefix}_annotated_cons.vcf"), emit: vcf
    path "versions.yml"                                  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python2 ${projectDir}/bin/inc_add_conservation_v1.0.py \\
        ${vcf} \\
        ${prefix}_annotated_cons.vcf \\
        ${conservation_bed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: INC_VCF_TO_CSV
// ==========================================
// inc_vcf_to_csv_v1.2.py + inc_aa_3to1_v1.0.py → <corrida>_annotated.csv
// inc_aa_3to1 es python2, así que se ejecuta aparte (INC_AA_3TO1).
process INC_VCF_TO_CSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(vcf)
    path list_of_transcripts

    output:
    tuple val(meta), path("${prefix}_annotated.csv"), emit: csv
    path "versions.yml"                             , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tr_list = list_of_transcripts ? "-r ${list_of_transcripts}" : ''
    def amplicon = params.amplicon ? '--amplicon' : ''
    """
    python3 ${projectDir}/bin/inc_vcf_to_csv_v1.2.py \\
        -f ${vcf} \\
        -o ${prefix}_annotated.csv \\
        -m ${params.af} \\
        -n 1 \\
        ${tr_list} \\
        ${amplicon}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pyvcf3: \$(python3 -c 'import vcf; print(vcf.VERSION)')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: INC_AA_3TO1
// ==========================================
// inc_aa_3to1_v1.0.py (python2): HGVSp con aminoácidos de una letra. El script sobrescribe el csv de entrada,
// así que se trabaja sobre una copia en el directorio de salida.
process INC_AA_3TO1 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment_python2.yml"

    input:
    tuple val(meta), path(csv, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}_annotated.csv"), emit: csv
    path "versions.yml"                             , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp ${csv} ${prefix}_annotated.csv
    python2 ${projectDir}/bin/inc_aa_3to1_v1.0.py \$PWD/${prefix}_annotated.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: INC_SAMPLES_CONTAMINATION
// ==========================================
// inc_samples_contamination_v2.0.R → <corrida>_sample_contamination.csv (variantes compartidas entre muestras,
// a partir de las columnas *_FREQUENCY del csv anotado)
process INC_SAMPLES_CONTAMINATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(csv)

    output:
    path "QC/mapping/*_sample_contamination.csv", emit: csv
    path "versions.yml"                         , emit: versions

    script:
    // El script lee <dir>/variant_calling/*annotated.csv y escribe en <dir>/QC/mapping/
    """
    mkdir -p variant_calling QC/mapping
    ln -s \$PWD/${csv} variant_calling/${csv.name}

    Rscript ${projectDir}/bin/inc_samples_contamination_v2.0.R \$PWD/ ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: INC_SPLIT_ANNOTATED_CSV
// ==========================================
// inc_split_annotated_csv_v1.0.R → un csv por muestra (<muestra>.csv, con '.' en lugar de '-' en el nombre)
process INC_SPLIT_ANNOTATED_CSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(csv)

    output:
    path "split/*.csv"  , emit: csv
    path "versions.yml" , emit: versions

    script:
    """
    mkdir split
    Rscript ${projectDir}/bin/inc_split_annotated_csv_v1.0.R -c ${csv} -o split

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """
}

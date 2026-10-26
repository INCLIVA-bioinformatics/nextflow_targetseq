// ==========================================
// Estadísticas de rutina.sh (función statistics() y final de annotate_variants())
// ==========================================
// Los cálculos reproducen línea a línea los de rutina.sh. Los scripts de rutina están en bin/
// (copiados de scripts_rutina_stats/ con un único cambio: el nombre de la corrida o del panel se
// pasa por argumento, porque en rutina se deduce de la ruta /media/scratch1/<corrida>/...).
//
// Sin contenedor: se ejecutan con conda (igual que GET_COVERED_BED).


// ==========================================
// PROCESS: GLOBAL_MAPPING_STATS
// ==========================================
// Una fila del <corrida>_global_mapping_stats.csv por muestra + ficheros de QC/mapping por muestra:
//   - <sample>_non_covered_regions_<cov>x.bed
//   - <sample>_coverage_plot.pdf          (inc_calculate_IQR_and_coverage_plot_v1.0.R)
//   - <sample>_af_plot.png                (inc_af_v1.0.R)
//   - <sample>_flagstat.txt               (inc_per_read_mapping.py)
//   - <sample>_mapped_by_chromosome.png y <sample>_chrY_chrX.png (inc_chr_mapped_reads_plots.R, si xy_stats)
process GLOBAL_MAPPING_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(raw_r1), path(prinseq_r1), path(sorted_metrics), path(ontarget_sorted_metrics), path(ontarget_processed_metrics), path(insert_metrics), path(dup_metrics), path(bam), path(bai), path(depth), path(non_covered), path(vcf)
    path target_merged
    path inc_db

    output:
    tuple val(meta), path("${meta.id}_mapping_stats.csv")                  , emit: csv
    tuple val(meta), path("${meta.id}_non_covered_regions_${params.cov}x.bed"), emit: non_covered_regions
    tuple val(meta), path("${meta.id}_coverage_plot.pdf")                  , emit: coverage_plot
    tuple val(meta), path("${meta.id}_af_plot.png")                        , emit: af_plot, optional: true
    tuple val(meta), path("${meta.id}_flagstat.txt")                       , emit: flagstat
    tuple val(meta), path("${meta.id}_idxstats")                           , emit: idxstats
    tuple val(meta), path("${meta.id}_{mapped_by_chromosome,chrY_chrX}.png"), emit: xy_plots, optional: true
    path "versions.yml"                                                    , emit: versions

    script:
    def sample = meta.id
    def cov    = params.cov
    def header = "Sample,Rawdata,%LowQReads,%MappedReads,%Paired_reads,%DuplicateReads,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,Nt_100x,Nt_250x,Nt_500x,Mean_insert,Insert_SD,NumVars"
    if (params.xy_stats) {
        header += ",chrY/chrX,chrX_mapped_reads,chrY_mapped_reads"
    }
    """
    samtools idxstats ${bam} > ${sample}_idxstats

    mapping_reads=`head ${sorted_metrics} | tail -1 | awk '{print \$6}'`
    ontarget_mapping_reads=`head ${ontarget_sorted_metrics} | tail -1 | awk '{print \$6}'`
    dedup_reads=`head ${ontarget_processed_metrics} | tail -1 | awk '{print \$6}'`

    raw_data=`zcat -f ${raw_r1} | wc -l | awk '{print \$1/2}'`

    if [ "${params.has_umi}" == "true" ]
    then
        tmp_duplicates=`echo "\${ontarget_mapping_reads} - \${dedup_reads}" | bc`
        duplicates=`echo "scale=4; \${tmp_duplicates}/\${raw_data}" | bc | awk '{print \$1*100}'`
    elif [ "${params.amplicon}" == "true" ]
    then
        duplicates="0"
    else
        duplicates=`head -8 ${dup_metrics} | tail -1 | awk '{print \$10*100}'`
    fi

    # inc_per_read_mapping.py: % de "properly paired" de samtools flagstat sobre el BAM final
    samtools flagstat ${bam} > ${sample}_flagstat.txt
    per_read=`grep "properly paired" ${sample}_flagstat.txt | cut -d'(' -f2 | cut -d'%' -f1`

    all_mapped=`echo "scale=4; \${mapping_reads}/\${raw_data}" | bc | awk '{print \$1*100}'`

    prinseq_reads=`zcat -f ${prinseq_r1} | wc -l | awk '{print \$1/2}'`
    prinseq_per=`echo "scale=4; \${prinseq_reads}/\${raw_data}" | bc | awk '{print 100-(\$1*100)}'`

    bam_target_reads=`head ${ontarget_sorted_metrics} | tail -1 | awk '{print \$6}'`
    usable_reads=`echo "scale=4; \${dedup_reads}/\${raw_data}" | bc | awk '{print \$1*100}'`
    specificity=`echo "scale=4; \${bam_target_reads}/\${mapping_reads}" | bc | awk '{print \$1*100}'`

    all_vars=`zcat ${vcf} | grep -v "#" -c || true`

    mean_insert_tmp=`head -8 ${insert_metrics} | tail -1 | awk '{print \$6}'`
    mean_insert=`echo "scale=2; \${mean_insert_tmp}/1" | bc`
    deviation_tmp=`head -8 ${insert_metrics} | tail -1 | awk '{print \$7}'`
    deviation=`echo "scale=2; \${deviation_tmp}/1" | bc`

    design=`awk '{print \$3-\$2+1}' ${target_merged} | awk 'BEGIN {sum=0} {sum+=\$1} END {print sum}'`

    # ${non_covered} = bedtools subtract -a target.bed -b <sample>_covered_regions (GET_COVERED_BED)
    bedtools intersect -a ${non_covered} -b ${inc_db} -wao | awk '{if (\$5 ~ /chr/) print \$12"\\t"\$1"\\t"\$2"\\t"\$3"\\t"\$4; else print "--\\t"\$1"\\t"\$2"\\t"\$3"\\t"\$4}' | sed 's/.*CDS=//' | sed 's/;CNT.*\\tchr/\\tchr/' | awk '{print \$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$1}' > ${sample}_non_covered_regions_${cov}x.bed

    sample250=`awk '{if (\$3>=250) print \$1"\\t"\$2-1"\\t"\$2}' ${depth} | bedtools merge | awk '{print \$3-\$2+1}' | awk 'BEGIN {sum=0} {sum+=\$1} END {print sum}'`
    depth250=`echo "scale=4; \${sample250}/\${design}" | bc | awk '{print \$1*100}'`
    sample100=`awk '{if (\$3>=100) print \$1"\\t"\$2-1"\\t"\$2}' ${depth} | bedtools merge | awk '{print \$3-\$2+1}' | awk 'BEGIN {sum=0} {sum+=\$1} END {print sum}'`
    depth100=`echo "scale=4; \${sample100}/\${design}" | bc | awk '{print \$1*100}'`
    sample500=`awk '{if (\$3>=500) print \$1"\\t"\$2-1"\\t"\$2}' ${depth} | bedtools merge | awk '{print \$3-\$2+1}' | awk 'BEGIN {sum=0} {sum+=\$1} END {print sum}'`
    depth500=`echo "scale=4; \${sample500}/\${design}" | bc | awk '{print \$1*100}'`

    Rscript ${projectDir}/bin/inc_calculate_IQR_and_coverage_plot_v1.0.R ${depth} ${sample}_R_coverage_summary ${sample}_coverage_plot.pdf ${sample}

    Rscript ${projectDir}/bin/inc_af_v1.0.R ${vcf} .

    first_quantile=`tail -1 ${sample}_R_coverage_summary | awk '{print \$2}'`
    third_quantile=`tail -1 ${sample}_R_coverage_summary | awk '{print \$5}'`
    median_coverage=`tail -1 ${sample}_R_coverage_summary | awk '{print \$3}'`

    echo "${header}" > ${sample}_mapping_stats.csv

    if [ "${params.xy_stats}" == "true" ]
    then
        chrX=`grep chrX ${sample}_idxstats | cut -f3`
        chrY=`grep chrY ${sample}_idxstats | cut -f3`

        if [ \${chrX} == 0 ]
        then
            chrY_chrX="No chrX reads"
        else
            chrY_chrX=\$(awk "BEGIN { print (\$chrY / \$chrX) }")
        fi

        echo "${sample},\${raw_data},\${prinseq_per},\${all_mapped},\${per_read},\${duplicates},\${usable_reads},\${specificity},\${first_quantile},\${third_quantile},\${median_coverage},\${depth100},\${depth250},\${depth500},\${mean_insert},\${deviation},\${all_vars},\${chrY_chrX},\${chrX},\${chrY}" >> ${sample}_mapping_stats.csv

        Rscript ${projectDir}/bin/inc_chr_mapped_reads_plots.R ${sample}_idxstats .
    else
        echo "${sample},\${raw_data},\${prinseq_per},\${all_mapped},\${per_read},\${duplicates},\${usable_reads},\${specificity},\${first_quantile},\${third_quantile},\${median_coverage},\${depth100},\${depth250},\${depth500},\${mean_insert},\${deviation},\${all_vars}" >> ${sample}_mapping_stats.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        r-base: \$(Rscript --version 2>&1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: STATS_BY_GENE
// ==========================================
// inc_gene_cover_v1.1.py → <corrida>_stats_by_gene.tsv (todas las muestras en una tabla).
// Se le da el mismo <sample>_depth que usa rutina (samtools depth sobre target.bed ordenado y
// fusionado), así que el script no vuelve a calcular la profundidad.
process STATS_BY_GENE {
    tag "${params.study_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(ids), path(depths)
    path global_stats
    path panel_design
    path inc_db

    output:
    path "analysis/QC/mapping/${params.study_name}_stats_by_gene.tsv", emit: tsv
    path "versions.yml"                                                , emit: versions

    script:
    def name  = params.study_name
    def links = [ids, [depths].flatten()].transpose().collect { id, depth -> "ln -s \$PWD/${depth} tmp/mapping/${id}_depth" }.join('\n    ')
    """
    mkdir -p analysis/QC/mapping analysis/mapping tmp/mapping
    cp ${global_stats} analysis/QC/mapping/
    ${links}

    python3 ${projectDir}/bin/inc_gene_cover_v1.1.py \\
        -a \$PWD/analysis \\
        -t \$PWD/tmp \\
        -p \$PWD/${panel_design} \\
        -c ${params.cov} \\
        -d \$PWD/${inc_db} \\
        -s samtools \\
        -b bedtools \\
        -n ${name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: NON_COVERED_REGIONS
// ==========================================
// inc_non_covered_regions_v1.1.py → <corrida>_common_non_covered_regions_<cov>.bed
// (posiciones con profundidad <= cov en >= 70% de las muestras con Nt_100x >= 75).
process NON_COVERED_REGIONS {
    tag "${params.study_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(ids), path(bams), path(bais)
    path global_stats
    path target
    path inc_db

    output:
    path "analysis/QC/mapping/${params.study_name}_common_non_covered_regions_${params.cov}.bed", emit: bed
    path "versions.yml"                                                                         , emit: versions

    script:
    def name  = params.study_name
    def links = [ids, [bams].flatten(), [bais].flatten()].transpose().collect { id, bam, bai ->
        "ln -s \$PWD/${bam} analysis/mapping/${id}.bam\n    ln -s \$PWD/${bai} analysis/mapping/${id}.bam.bai"
    }.join('\n    ')
    """
    mkdir -p analysis/QC/mapping analysis/mapping tmp/mapping
    cp ${global_stats} analysis/QC/mapping/
    ${links}

    python3 ${projectDir}/bin/inc_non_covered_regions_v1.1.py \\
        -g analysis/QC/mapping/${global_stats} \\
        -a \$PWD/analysis \\
        -t \$PWD/tmp \\
        -r \$PWD/${target} \\
        -c ${params.cov} \\
        -d \$PWD/${inc_db} \\
        -s samtools \\
        -b bedtools \\
        -n ${name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}


// ==========================================
// PROCESS: STATS_BY_AMPLICON
// ==========================================
// Sólo paneles de amplicones. Igual que rutina.sh:
//   - BRCA (params.brca): inc_amplicons_BRCA_universal_v1.1.py → <corrida>_BRCA_stats_by_amplicon.tsv
//     (v1.1 = v1.0 de rutina con la misma salida, pero recorre el SAM una vez por muestra y no una vez por amplicón)
//   - TP53 (params.panel == 'TP53'): inc_amplicons_TP53_universal_v1.1.py → <corrida>_TP53_stats_by_amplicon.tsv
//     (v1.1 = v1.0 de rutina con la misma salida: misma optimización y el diseño del panel como argumento)
//   - resto:              inc_amplicons_cover_v1.2_mod.py      → <corrida>_<panel>_stats_by_amplicon.tsv
// Todos los scripts son python2.
// Los scripts escriben la cabecera con os.system("echo '...\\t...'"): en el servidor de rutina /bin/sh es
// dash y su echo convierte \t en tabulador; aquí /bin/sh es bash, así que se activa xpg_echo para que haga lo mismo.
process STATS_BY_AMPLICON {
    tag "${params.study_name}"
    label 'process_low'

    conda "${moduleDir}/environment_python2.yml"

    input:
    tuple val(ids), path(bams), path(bais)
    path panel_design

    output:
    path "${params.study_name}_*_stats_by_amplicon.tsv", emit: tsv
    path "versions.yml"            , emit: versions

    script:
    def name = params.study_name
    // Orden alfabético, como "ls ${analysis_dir}/mapping/*bam" en rutina: la primera muestra lleva las
    // columnas del amplicón y el resto se pegan detrás (paste de los ficheros ordenados).
    def samples = [ids, [bams].flatten(), [bais].flatten()].transpose().sort { it[0] }
    if (params.brca || params.panel == 'TP53') {
        // rutina.sh: BRCA = yes → script de BRCA; si no, panel TP53 → script de TP53
        def type = params.brca ? 'BRCA' : 'TP53'
        def links = samples.collect { id, bam, bai ->
            "ln -s \$PWD/${bam} ${name}/analysis/mapping/${id}.bam\n    ln -s \$PWD/${bai} ${name}/analysis/mapping/${id}.bam.bai"
        }.join('\n    ')
        // El script saca el nombre de la corrida de la ruta: <...>/<corrida>/<analysis>/
        """
        mkdir -p ${name}/analysis/mapping ${name}/analysis/QC/mapping
        ${links}

        env BASHOPTS=xpg_echo python2 ${projectDir}/bin/inc_amplicons_${type}_universal_v1.1.py \\
            \$PWD/${name}/analysis \\
            ${samples[0][0]} \\
            ${params.cov} \\
            \$PWD/${panel_design}

        cp ${name}/analysis/QC/mapping/${name}_${type}_stats_by_amplicon.tsv .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python2 --version 2>&1 | sed 's/Python //')
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """
    } else {
        def links = samples.collect { id, bam, bai ->
            "ln -s \$PWD/${bam} mapping/${id}.bam\n    ln -s \$PWD/${bai} mapping/${id}.bam.bai"
        }.join('\n    ')
        def bam_list = samples.collect { id, bam, bai -> "\$PWD/mapping/${id}.bam" }.join(',')
        """
        mkdir -p mapping QC/mapping
        ${links}

        env BASHOPTS=xpg_echo python2 ${projectDir}/bin/inc_amplicons_cover_v1.2_mod.py \\
            \$PWD/QC/mapping \\
            \$PWD/${panel_design}/capture_amplicons.bed \\
            ${params.cov} \\
            ${bam_list} \\
            ${name} \\
            --panel_name ${params.panel}

        cp QC/mapping/${name}_${params.panel}_stats_by_amplicon.tsv .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python2 --version 2>&1 | sed 's/Python //')
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """
    }
}


// ==========================================
// PROCESS: AF_MERGE_PLOT
// ==========================================
// inc_af_merge_v1.0.R → <corrida>_af_plot_<n>.png (densidad de AF, 6 muestras por gráfico).
process AF_MERGE_PLOT {
    tag "${params.study_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(ids), path(vcfs)
    path global_stats

    output:
    path "QC/mapping/*_af_plot_*.png", emit: png, optional: true
    path "versions.yml"              , emit: versions

    script:
    def links = [ids, [vcfs].flatten()].transpose().collect { id, vcf -> "ln -s \$PWD/${vcf} variant_calling/${id}.vcf.gz" }.join('\n    ')
    """
    mkdir -p variant_calling QC/mapping
    cp ${global_stats} QC/mapping/
    ${links}

    Rscript ${projectDir}/bin/inc_af_merge_v1.0.R \$PWD/variant_calling \$PWD/QC/mapping ${params.study_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """
}

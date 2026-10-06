// ==========================================
// PROCESS: GLOBAL_MAPPING_STATS
// ==========================================
// Este proceso concentra TODAS las métricas globales de mapeo
// y genera un CSV por muestra con valores consolidados.
// 
// Entradas:
//   - CRAM final
//   - stats de Samtools (.sorted.cram.stats)
//   - cobertura por base y thresholds de Mosdepth
//   - Fastp JSON y FASTQ (para Rawdata y calidad)
//   - bcftools stats (Mutect2)
//   - FASTA + FAI de referencia
//   - BED de captura
//
// Output: <sample>_mapping_stats.csv
// Contiene: Rawdata, %LowQReads, %MappedReads, %Paired_reads,
// %DuplicateReads, %OnTargetNoDupReads, Kit_specificity, cuartiles,
// cobertura ≥100/250/500x, Mean insert, Insert SD, NumVars.
process GLOBAL_MAPPING_STATS {
    tag "${meta.sample}"
    container "${projectDir}/containers/statistics.sif"

    //container "${params.container_global_stats}"

    input:
    tuple val(meta), path(cram), path(sorted_stats), path(perbase), path(thr), path(fastp), path(bcftools), path(fastq1)
    val fasta
    val fai
    path bed

    output:
    path "${meta.sample}_mapping_stats.csv"

    script:
    """
    echo "▶️ [GLOBAL_MAPPING_STATS] Procesando muestra ${meta.sample}" >&2

    # Nombre real de la muestra desde FASTQ R1
    sample=${meta.sample}

    # 1. Estadísticas básicas (incluyen Rawdata ahora)
    summarize_fastq_cram_stats.sh ${fastp} ${sorted_stats} ${fastq1} ${meta.sample}.fastq_summary.txt
    rawdata=\$(grep "^Rawdata" ${meta.sample}.fastq_summary.txt | cut -f2)
    lowq=\$(grep "^%LowQReads" ${meta.sample}.fastq_summary.txt | cut -f2)
    mapped=\$(grep "^%MappedReads" ${meta.sample}.fastq_summary.txt | cut -f2)
    paired=\$(grep "^%Paired_reads" ${meta.sample}.fastq_summary.txt | cut -f2)
    dup=\$(grep "^%DuplicateReads" ${meta.sample}.fastq_summary.txt | cut -f2)

    # 2. Cuartiles de cobertura
    extract_coverage_quartiles.py -p ${perbase} -c ${bed} -o ${meta.sample}.coverage_quartiles.txt
    cov_q1=\$(grep "^Cov_1stQ" ${meta.sample}.coverage_quartiles.txt | cut -f2)
    cov_q3=\$(grep "^Cov_3rdQ" ${meta.sample}.coverage_quartiles.txt | cut -f2)
    cov_med=\$(grep "^Cov_Median" ${meta.sample}.coverage_quartiles.txt | cut -f2)

    # 3. Cobertura ≥100x/250x/500x
    coverage_depth_thresholds_100_250_500.py -i ${thr} -o ${meta.sample}.coverage_thresholds.txt
    nt100=\$(grep "^Nt_100x" ${meta.sample}.coverage_thresholds.txt | cut -f2)
    nt250=\$(grep "^Nt_250x" ${meta.sample}.coverage_thresholds.txt | cut -f2)
    nt500=\$(grep "^Nt_500x" ${meta.sample}.coverage_thresholds.txt | cut -f2)

    # 4. Inserto medio y SD (ajuste para nuevo script)
    calculate_insert_metrics.py -j ${fastp} -o ${meta.sample}.insert_metrics.txt
    insert_mean=\$(grep "^Mean_insert" ${meta.sample}.insert_metrics.txt | cut -f2)
    insert_sd=\$(grep "^Insert_SD" ${meta.sample}.insert_metrics.txt | cut -f2)

    # 5. Número de variantes
    count_filtered_variants.sh ${bcftools} ${meta.sample}.variant_count.txt
    numvars=\$(awk 'NR==1{print \$1}' ${meta.sample}.variant_count.txt)

    # 6. %OnTargetNoDupReads
    calculate_ontarget_picard.sh ${cram} ${fastq1} ${fasta} ${meta.sample} .
    ontarget=\$(awk 'NR==1{print \$NF}' ${meta.sample}.ontarget_result.txt)

    # 7. Kit Specificity
    calculate_kit_specificity_picard.sh ${cram} ${fasta} ${bed} ${meta.sample} .
    kit=\$(awk 'NR==1{print \$NF}' ${meta.sample}.kit_specificity_result.txt)

    # 8. Escribir CSV final
    echo "Sample,Rawdata,%LowQReads,%MappedReads,%Paired_reads,%DuplicateReads,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,Nt_100x,Nt_250x,Nt_500x,Mean_insert,Insert_SD,NumVars" > ${meta.sample}_mapping_stats.csv
    echo "${meta.sample},\${rawdata},\${lowq},\${mapped},\${paired},\${dup},\${ontarget},\${kit},\${cov_q1},\${cov_q3},\${cov_med},\${nt100},\${nt250},\${nt500},\${insert_mean},\${insert_sd},\${numvars}" >> ${meta.sample}_mapping_stats.csv
    """
}

// ==========================================
// PROCESS: STATS_BY_GENE
// ==========================================
// Calcula métricas por gen (cobertura y variantes).
// Usa pybedtools + pysam sobre per-base de Mosdepth + BED de genes.
// Opcionalmente añade variantes no cubiertas desde un VCF.
//
// Output: <sample>_stats_by_gene.tsv
process STATS_BY_GENE {
    tag "${meta.sample}"
    // container "${params.container_stats}"
    container "${projectDir}/containers/statistics.sif"

      
    input:
    tuple val(meta), path(perbase), path(bed), path(incdb)

    output:
    path "${meta.sample}_stats_by_gene.tsv"

    script:
    """
    python /opt/statistics/bin/create_stats_by_gene.py \
      -p ${perbase} \
      -c ${bed} \
      -i ${incdb} \
      -s ${meta.sample} \
      -o ${meta.sample}_stats_by_gene.tsv
    """
}


// ==========================================
// PROCESS: STATS_BY_AMPLICON
// ==========================================
// Calcula métricas de cobertura por amplicón a partir de Mosdepth
// y un archivo BED de amplicones.
//
// Output: <sample>_stats_by_amplicon.tsv
process STATS_BY_AMPLICON {
    tag "${meta.sample}"
    //container "${params.container_stats}"
    container "${projectDir}/containers/statistics.sif"

    
    input:
    tuple val(meta), path(perbase), path(amplicon_bed)

    output:
    path "${meta.sample}_stats_by_amplicon.tsv"

    script:
    """
    python /opt/statistics/bin/create_stats_by_amplicon.py \
      -p ${perbase} \
      -c ${amplicon_bed} \
      -s ${meta.sample} \
      -o ${meta.sample}_stats_by_amplicon.tsv
    """
}


// ==========================================
// PROCESS: STATS_ANNOTATED (opcional)
// ==========================================
// Convierte un VCF anotado por VEP en un TSV con campos legibles.
// Incluye: Gene, CHROM, POS, REF, ALT, Consequence, IMPACT,
// posiciones de proteína, codones, variantes existentes, SIFT, PolyPhen...
//
// Output: <sample>_annotated.tsv
process STATS_ANNOTATED {
    tag "${meta.sample}"
    //container "${params.container_stats}"
    container "${projectDir}/containers/statistics.sif"


    input:
    tuple val(meta), path(ann_vcf)

    output:
    path "${meta.sample}_annotated.tsv"

    script:
    """
    python /opt/statistics/bin/create_annotated_tsv.py \
      -i ${ann_vcf} \
      -o ${meta.sample}_annotated.tsv \
      -s ${meta.sample}
    """
}

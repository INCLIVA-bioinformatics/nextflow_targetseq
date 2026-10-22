// ===============================
// SUBWORKFLOW: STATISTICS
// ===============================
// Reproduce las estadísticas de rutina.sh (función statistics() y final de annotate_variants()):
//   - GLOBAL_MAPPING_STATS → <corrida>_global_mapping_stats.csv + ficheros por muestra en QC/mapping
//   - STATS_BY_GENE        → <corrida>_stats_by_gene.tsv            (inc_gene_cover_v1.1.py)
//   - NON_COVERED_REGIONS  → <corrida>_common_non_covered_regions_<cov>.bed (inc_non_covered_regions_v1.1.py)
//   - STATS_BY_AMPLICON    → <corrida>_<panel>_stats_by_amplicon.tsv (sólo amplicones)
//   - AF_MERGE_PLOT        → <corrida>_af_plot_<n>.png              (inc_af_merge_v1.0.R)
//
// Pendiente (depende de la anotación, que aún no está en el pipeline):
//   - inc_samples_contamination_v2.0.R (lee <corrida>_annotated.csv)

include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS as PICARD_CASM_SORTED             } from '../modules/incliva/picard/collectalignmentsummarymetrics/main'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS as PICARD_CASM_ONTARGET_SORTED    } from '../modules/incliva/picard/collectalignmentsummarymetrics/main'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS as PICARD_CASM_ONTARGET_PROCESSED } from '../modules/incliva/picard/collectalignmentsummarymetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS                                         } from '../modules/incliva/picard/collectinsertsizemetrics/main'
include { GLOBAL_MAPPING_STATS                                                    } from '../modules/incliva/statistics/main'
include { STATS_BY_GENE                                                           } from '../modules/incliva/statistics/main'
include { NON_COVERED_REGIONS                                                     } from '../modules/incliva/statistics/main'
include { STATS_BY_AMPLICON                                                       } from '../modules/incliva/statistics/main'
include { AF_MERGE_PLOT                                                           } from '../modules/incliva/statistics/main'

workflow STATISTICS {
    take:
    raw_reads              // [meta, fastq_1, fastq_2]  FASTQ crudos (Rawdata)
    preprocessed_reads     // [meta, [r1, r2]]          FASTQ tras prinseq (%LowQReads)
    sorted_bam             // [meta, bam]               <sample>_sorted.bam
    ontarget_sorted_bam    // [meta, bam]               ontarget_<sample>_sorted.bam
    ontarget_processed_bam // [meta, bam]               ontarget_<sample>_processed.bam
    dedup_metrics          // [meta, metrics]           MarkDuplicates (vacío en amplicones)
    postprocessed_bam      // [meta, bam]               BAM final (analysis/mapping/<sample>.bam)
    postprocessed_bai      // [meta, bai]
    depth                  // [meta, tsv]               <sample>_depth
    non_covered            // [meta, bed]               <sample>_non_covered_<cov>x.bed
    vcf                    // [meta, vcf.gz]            VCF combinado por muestra
    target_merged_bed      // [meta, bed]               target.bed ordenado y fusionado
    fasta
    fasta_fai
    dict
    versions

    main:
    def inc_db = file(params.inc_db, checkIfExists: true)
    def panel_design = file(params.panel_design, checkIfExists: true)
    def target_bed = file(params.target_intervals, checkIfExists: true)

    // Métricas de Picard (mismos BAM que rutina.sh)
    PICARD_CASM_SORTED(sorted_bam, fasta, fasta_fai, dict)
    PICARD_CASM_ONTARGET_SORTED(ontarget_sorted_bam, fasta, fasta_fai, dict)
    PICARD_CASM_ONTARGET_PROCESSED(ontarget_processed_bam, fasta, fasta_fai, dict)
    PICARD_COLLECTINSERTSIZEMETRICS(ontarget_processed_bam)

    versions = versions
        .mix(PICARD_CASM_SORTED.out.versions.first())
        .mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions.first())

    def dup = params.amplicon
        ? postprocessed_bam.map { meta, bam -> [ meta.id, [] ] }
        : byId(dedup_metrics)

    stats_input = postprocessed_bam.map { meta, bam -> [ meta.id, meta ] }
        .join(raw_reads.map { meta, fastq_1, fastq_2 -> [ meta.id, fastq_1 ] })
        .join(preprocessed_reads.map { meta, reads -> [ meta.id, reads[0] ] })
        .join(byId(PICARD_CASM_SORTED.out.metrics))
        .join(byId(PICARD_CASM_ONTARGET_SORTED.out.metrics))
        .join(byId(PICARD_CASM_ONTARGET_PROCESSED.out.metrics))
        .join(byId(PICARD_COLLECTINSERTSIZEMETRICS.out.metrics))
        .join(dup)
        .join(byId(postprocessed_bam))
        .join(byId(postprocessed_bai))
        .join(byId(depth))
        .join(byId(non_covered))
        .join(byId(vcf))
        .map { row -> row[1..-1] } // quitar la clave: [meta, raw_r1, prinseq_r1, ...]

    GLOBAL_MAPPING_STATS(
        stats_input,
        target_merged_bed.map { meta, bed -> bed },
        inc_db
    )
    versions = versions.mix(GLOBAL_MAPPING_STATS.out.versions.first())

    global_csv = GLOBAL_MAPPING_STATS.out.csv
        .map { meta, csv -> csv }
        .collectFile(
            name: "${params.study_name}_global_mapping_stats.csv",
            keepHeader: true,
            sort: true,
            storeDir: "${params.outdir}/QC/mapping"
        )

    bams = postprocessed_bam.join(postprocessed_bai)
        .map { meta, bam, bai -> [ meta.id, bam, bai ] }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows -> [ rows*.get(0), rows*.get(1), rows*.get(2) ] }

    STATS_BY_GENE(collectById(depth), global_csv, panel_design, inc_db)
    versions = versions.mix(STATS_BY_GENE.out.versions)

    NON_COVERED_REGIONS(bams, global_csv, target_bed, inc_db)
    versions = versions.mix(NON_COVERED_REGIONS.out.versions)

    if (params.amplicon) {
        STATS_BY_AMPLICON(bams, panel_design)
        versions = versions.mix(STATS_BY_AMPLICON.out.versions)
    }

    AF_MERGE_PLOT(collectById(vcf), global_csv)
    versions = versions.mix(AF_MERGE_PLOT.out.versions)

    emit:
    global_csv
    versions
}

// [meta, fichero] -> [meta.id, fichero], para unir canales por muestra
def byId(ch) {
    return ch.map { meta, f -> [ meta.id, f ] }
}

// [ [ids], [ficheros] ] ordenado por id, para los procesos que trabajan con todas las muestras
def collectById(ch) {
    return ch
        .map { meta, f -> [ meta.id, f ] }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows -> [ rows*.get(0), rows*.get(1) ] }
}

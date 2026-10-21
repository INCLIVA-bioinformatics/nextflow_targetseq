/*
 * Mapping and post-processing subworkflow for TARGETSEQ.
 *
 * Compatible target:
 *   Nextflow 24.x
 *   DSL2 classic parser
 *
 * This file still imports modules from the current central Sarek/INCLIVA
 * installation. Once the workflow is stable, copy/version these modules
 * inside this repository and replace absolute paths with relative paths.
 */


include { BWAMEM2_MEM as BWAMEM2 } from '../modules/nf-core/bwamem2/mem/main'

include { PICARD_CLEANSAM } from '../modules/incliva/picard/cleansam/main'

include { PICARD_MARKDUPLICATES } from '../modules/incliva/picard/markduplicates/main'

include { BAMCLIPPER } from '../modules/incliva/bamclipper/main'

include { SAMTOOLS_INDEX as INDEX_SORTED_BAM } from '../modules/nf-core/samtools/index/main'

include { SAMTOOLS_DEPTH } from '../modules/incliva/samtools/depth/main'

include { SAMTOOLS_INDEX as INDEX_PROCESSED_BAM } from '../modules/nf-core/samtools/index/main'

include { SAMTOOLS_INDEX as INDEX_POSTPROCESSED_BAM } from '../modules/nf-core/samtools/index/main'

include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_ONTARGET } from '../modules/nf-core/samtools/view/main'

include { GATK4_BASERECALIBRATOR } from '../modules/nf-core/gatk4/baserecalibrator/main'

include { GATK4_APPLYBQSR } from '../modules/nf-core/gatk4/applybqsr/main'

include { LOFREQ_INDELQUAL } from '../modules/incliva/lofreq/indelqual/main'

include { BEDTOOLS_SORT } from '../modules/incliva/bedtools/sort/main'

include { BEDTOOLS_MERGE } from '../modules/incliva/bedtools/merge/main'

include { GET_COVERED_BED } from '../incliva_modules'


workflow MAPPING_POSTPROCESSING {

    take:
    preprocessed_reads
    bwamem2
    intervals
    target
    fasta
    fasta_fai
    dict
    known_sites
    known_sites_tbi
    versions

    main:
    /*
     * PASO 1: Mapeado con BWA-MEM2. MAPEADO DE LAS READS CONTRA UN GENOMA DE REFERENCIA
     */

    sort_bam = true // reads will be sorted

    index = bwamem2.map { index_files ->
        [
            [id: 'index'],
            index_files
        ]
    }

    fasta_reference = fasta.map { fasta_file ->
        [
            [id: fasta_file.baseName],
            fasta_file
        ]
    }

    fasta_meta           = fasta.map { f -> [ [id:'fasta'], f ] }
    fasta_fai_meta       = fasta_fai.map { f -> [ [id:'fai'], f ] }
    dict_meta            = dict.map { f -> [ [id:'dict'], f ] }
    known_sites_meta     = known_sites.map { f -> [ [id:'known_sites'], f ] }
    known_sites_tbi_meta = known_sites_tbi.map { f -> [ [id:'known_sites_tbi'], f ] }

    // Mapeado contra genoma de referencia
    BWAMEM2(
        preprocessed_reads,
        index,
        fasta_reference,
        sort_bam
    )

    sorted_bam = BWAMEM2.out.bam.map { meta, bam ->
        [
            meta + [data_type: 'bam'] - meta.subMap('read_group'),
            bam
        ]
    }

    versions = versions.mix(BWAMEM2.out.versions.first())

    /*
     * PASO 2: Postprocesado del BAM
     */

    // TODO - UMIs: rutina.sh usa umi_tools dedup cuando has_umi = yes
    if (params.amplicon) {
        // Amplicones: soft-clipping de primers con BAMClipper (no se eliminan duplicados)
        INDEX_SORTED_BAM(sorted_bam)
        versions = versions.mix(INDEX_SORTED_BAM.out.versions.first())

        BAMCLIPPER(
            sorted_bam.join(INDEX_SORTED_BAM.out.bai, failOnMismatch: true),
            file(params.amplicons_bedpe, checkIfExists: true)
        )

        deduplicated_bam = BAMCLIPPER.out.bam
        versions = versions.mix(BAMCLIPPER.out.versions.first())
    } else {
        // Captura: eliminación de duplicados con Picard MarkDuplicates
        PICARD_MARKDUPLICATES(sorted_bam)

        deduplicated_bam = PICARD_MARKDUPLICATES.out.bam
        versions = versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
    }

    PICARD_CLEANSAM(deduplicated_bam)

    processed_bam = PICARD_CLEANSAM.out.bam
    versions = versions.mix(PICARD_CLEANSAM.out.versions.first())
    // Indexamos el archivo bam procesado
    INDEX_PROCESSED_BAM(processed_bam)

    processed_bai = INDEX_PROCESSED_BAM.out.bai // Actualizar meta.data_type ya que ahora tenemos bam
    versions = versions.mix(INDEX_PROCESSED_BAM.out.versions.first())

    processed_bam_bai = processed_bam.join(
        processed_bai,
        failOnMismatch: true
    )
    // Combine bam and intervals for spread and gather strategy
    processed_bam_bai_intervals = processed_bam_bai.combine(intervals)
    // Calcular tabla de recalibración
    GATK4_BASERECALIBRATOR(
        processed_bam_bai_intervals,
        fasta_meta,
        fasta_fai_meta,
        dict_meta,
        known_sites_meta,
        known_sites_tbi_meta
    )

    recalibration_table = GATK4_BASERECALIBRATOR.out.table
    versions = versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

    // ApplyBQSR sin intervalos (igual que rutina.sh): se recalibra el BAM completo
    bam_to_recalibrate = processed_bam_bai
        .join(recalibration_table, failOnMismatch: true)
        .map { meta, bam, bai, bqsr_table ->
            [
                meta,
                bam,
                bai,
                bqsr_table,
                []
            ]
        }
    // Aplicamos recalibración
    GATK4_APPLYBQSR(
        bam_to_recalibrate,
        fasta,
        fasta_fai,
        dict
    )

    recalibrated_bam = GATK4_APPLYBQSR.out.bam
    versions = versions.mix(GATK4_APPLYBQSR.out.versions.first())
    /// Add indel qualities
    LOFREQ_INDELQUAL(
        recalibrated_bam,
        fasta_reference
    )

    postprocessed_bam = LOFREQ_INDELQUAL.out.bam
    versions = versions.mix(LOFREQ_INDELQUAL.out.versions.first())
    // Index final bam
    INDEX_POSTPROCESSED_BAM(postprocessed_bam)

    postprocessed_bai = INDEX_POSTPROCESSED_BAM.out.bai
    versions = versions.mix(INDEX_POSTPROCESSED_BAM.out.versions.first())
   // En este punto ya se ha creado el bam recalibrado sobre el que llevar a cabo el variant calling (postprocessed_bam)
    qname = Channel.value([]) // channel required for readnames in samtools view
    index_format = Channel.value('bai')

    SAMTOOLS_VIEW_ONTARGET(
        processed_bam_bai,
        fasta_reference,
        qname,
        index_format
    )

    ontarget_processed_bam = SAMTOOLS_VIEW_ONTARGET.out.bam
    versions = versions.mix(SAMTOOLS_VIEW_ONTARGET.out.versions.first())
    // Se ordena el archivo target.bed
    target_to_sort = target.map { target_bed ->
        [
            [id: target_bed.baseName],
            target_bed
        ]
    }

    BEDTOOLS_SORT(target_to_sort)

    target_sorted = BEDTOOLS_SORT.out.sorted.map { meta, bed ->
        [
            [id: bed.baseName],
            bed
        ]
    }

    versions = versions.mix(BEDTOOLS_SORT.out.versions.first())
    // Se mergean los intervalos de target_sorted.bed para crear target_sorted_merged.bed
    BEDTOOLS_MERGE(target_sorted)
    
    target_sorted_merged = BEDTOOLS_MERGE.out.bed.map { meta, bed ->
        [
            [id: meta.id],
            bed
        ]
    }

    versions = versions.mix(BEDTOOLS_MERGE.out.versions.first())

    def bed_path = target_sorted_merged.first().last()  // Extraer el path

    // Calculamos la cobertura de las posiciones del panel en la muestra
    SAMTOOLS_DEPTH(
        postprocessed_bam,
        bed_path
    )


    sample_intervals_depth = SAMTOOLS_DEPTH.out.tsv
    versions = versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    // Se obtiene el covered_regions.bed de cada muestra (archivo bed con aquellas regiones del panel que superan la cobertura deseada)
    GET_COVERED_BED(sample_intervals_depth)

    covered_regions = GET_COVERED_BED.out.covered_regions
    versions = versions.mix(GET_COVERED_BED.out.versions.first())

    emit:
    postprocessed_bam
    postprocessed_bai
    ontarget_processed_bam
    covered_regions
    versions
}
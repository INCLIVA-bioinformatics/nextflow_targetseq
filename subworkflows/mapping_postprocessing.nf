include { BWAMEM2_MEM as BWAMEM2                    } from "${params.modules_path}/bwamem2/mem/main"
include { PICARD_CLEANSAM                           } from "${params.modules_path}/picard/cleansam/main"
include { SAMTOOLS_DEPTH                            } from "${params.modules_path}/samtools/depth/main"
include { SAMTOOLS_INDEX as INDEX_PROCESSED_BAM     } from "${params.modules_path}/samtools/index/main"
include { SAMTOOLS_INDEX as INDEX_POSTPROCESSED_BAM } from "${params.modules_path}/samtools/index/main"
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_ONTARGET   } from "${params.modules_path}/samtools/view/main"
include { GATK4_BASERECALIBRATOR                    } from "${params.modules_path}/gatk4/baserecalibrator/main"
include { GATK4_APPLYBQSR                           } from "${params.modules_path}/gatk4/applybqsr/main"
include { LOFREQ_INDELQUAL                          } from "${params.modules_path}/lofreq/indelqual/main"
include { BEDTOOLS_SORT                             } from "${params.modules_path}/bedtools/sort/main"
include { BEDTOOLS_MERGE                            } from "${params.modules_path}/bedtools/merge/main"
include { BEDTOOLS_INTERSECT                        } from "${params.modules_path}/bedtools/intersect/main"
include { GET_COVERED_BED                           } from "${params.incliva_modules_path}"

workflow MAPPING_POSTPROCESSING {
    take:
    preprocessed_reads // channel: [mandatory] [meta, [fastq1, fastq2]]
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
    // PASO 1: MAPEADO DE LAS READS CONTRA UN GENOMA DE REFERENCIA
    sort_bam = true // reads will be sorted
    index = bwamem2.map{ it -> [ [ id:'index' ], it ] }

    // Mapeado contra genoma de referencia
    BWAMEM2(preprocessed_reads, index, sort_bam)

    sorted_bam = BWAMEM2.out.bam.map{ meta, bam -> [ meta + [ data_type:'bam' ] - meta.subMap('read_group'), bam ] } // Actualizar meta.data_type ya que ahora tenemos bam y quitar read_group ya que ya se ha añadido al archivo del alineamiento
    versions = versions.mix(BWAMEM2.out.versions.first())

    // PASO 2: POSTPROCESADO
    // TODO - PROCESADO DE DUPLICADOS: UMIs = umitools, no UMIs + amplicon = bamclipper, no UMIs + sondas = picard Markduplicates

    PICARD_CLEANSAM(sorted_bam)

    processed_bam = PICARD_CLEANSAM.out.bam
    versions = versions.mix(PICARD_CLEANSAM.out.versions.first())

    // Indexamos el archivo bam procesado
    INDEX_PROCESSED_BAM(processed_bam)

    processed_bai = INDEX_PROCESSED_BAM.out.bai // Actualizar meta.data_type ya que ahora tenemos bam
    versions = versions.mix(INDEX_PROCESSED_BAM.out.versions.first())

    processed_bam_bai = processed_bam.join(processed_bai, failOnMismatch: true)

    // Combine bam and intervals for spread and gather strategy
    processed_bam_bai_intervals = processed_bam_bai.combine(intervals)

    // Calcular tabla de recalibración
    GATK4_BASERECALIBRATOR(processed_bam_bai_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)

    recalibration_table = GATK4_BASERECALIBRATOR.out.table
    versions = versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

    bam_to_recalibrate = processed_bam_bai_intervals.join(recalibration_table, failOnMismatch: true).map{ meta, bam, bai, intervals, bqsr_table -> [ meta, bam, bai, bqsr_table, intervals ]}

    // Aplicamos recalibración
    GATK4_APPLYBQSR(bam_to_recalibrate, fasta, fasta_fai, dict)

    recalibrated_bam = GATK4_APPLYBQSR.out.bam
    versions = versions.mix(GATK4_APPLYBQSR.out.versions.first())

    /// Add indel qualities
    LOFREQ_INDELQUAL(recalibrated_bam, fasta)

    postprocessed_bam = LOFREQ_INDELQUAL.out.bam
    versions = versions.mix(LOFREQ_INDELQUAL.out.versions.first())

    // Index final bam
    INDEX_POSTPROCESSED_BAM(postprocessed_bam)

    postprocessed_bai = INDEX_POSTPROCESSED_BAM.out.bai
    versions = versions.mix(INDEX_POSTPROCESSED_BAM.out.versions.first())

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    En este punto ya se ha creado el bam recalibrado sobre el que llevar a cabo el variant calling (postprocessed_bam)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // Se crea ontarget_sample_processed.bam (se usa más adelante) utilizando el bam antes del recalibrado (porque si no se pierden reads off-target)
    fasta_reference = fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] }
    qname = Channel.value([]) // channel required for readnames in samtools view

    SAMTOOLS_VIEW_ONTARGET(processed_bam_bai, fasta_reference, qname)

    ontarget_processed_bam = SAMTOOLS_VIEW_ONTARGET.out.bam
    versions = versions.mix(SAMTOOLS_VIEW_ONTARGET.out.versions.first())

    // Se ordena el archivo target.bed
    genome_file = Channel.value([]) // channel required for bedtools sort
    target_to_sort = target.map{ target -> [ [ id:target.baseName ], target ] }

    BEDTOOLS_SORT(target_to_sort, genome_file)

    target_sorted = BEDTOOLS_SORT.out.sorted.map{ meta, bed -> [ [ id:bed.baseName ], bed ] }
    versions = versions.mix(BEDTOOLS_SORT.out.versions.first())

    // Se mergean los intervalos de target_sorted.bed para crear target_sorted_merged.bed
    BEDTOOLS_MERGE(target_sorted)

    target_sorted_merged = BEDTOOLS_MERGE.out.bed.last().map{ meta, bed -> [ [ id:bed.baseName ], bed ] }
    versions = versions.mix(BEDTOOLS_MERGE.out.versions.first())

    // Calculamos la cobertura de las posiciones del panel en la muestra
    SAMTOOLS_DEPTH(postprocessed_bam, target_sorted_merged)

    sample_intervals_depth = SAMTOOLS_DEPTH.out.tsv
    versions = versions.mix(BEDTOOLS_MERGE.out.versions.first())

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

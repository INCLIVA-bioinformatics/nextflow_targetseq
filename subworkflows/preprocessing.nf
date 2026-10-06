include { FASTQC          as FASTQC_RAW                } from '../modules/nf-core/fastqc/main'
include { FASTQC          as FASTQC_PREPROCESSED       } from '../modules/nf-core/fastqc/main'

// include { UMITOOLS_EXTRACT                             } from "${params.modules_path}/umitools/extract/main"

include { CUTADAPT        as CUTADAPT_TRIMMING         } from '../modules/incliva/cutadapt/main'
include { CUTADAPT        as CUTADAPT_FILTER_BY_PRIMER } from '../modules/incliva/cutadapt/main'
include { PRINSEQ as PRINSEQ                   } from '../modules/incliva/prinseq/main' // TODO - Estudiar si utilizar FASTP en su lugar, ya que prinseq puede estar anticuado

workflow PREPROCESSING {
    take:
    samples // channel: [mandatory] [meta, fastq1, fastq2]
    reports // channel: []
    versions // channel: []

    main:
    // PASO 0: CONTROL DE CALIDAD & TRIMEADO
    fastq_samples = samples.map{ meta, fastq1, fastq2 -> [meta, [fastq1, fastq2]] } // reordenar canal para matchear cardinalidad

    // Raw QC 1
    FASTQC_RAW(fastq_samples)

    reports = reports.mix(FASTQC_RAW.out.zip.collect{ meta, logs -> logs })
    versions = versions.mix(FASTQC_RAW.out.versions.first())

    // TODO - raw QC 2: fastq_screen. Mirar #3157 de nfcore/modules.
    // FASTQ_SCREEN()

    // TODO - UMIs. ¿umitools o fgbio?

    // Trimeado de adaptadores y descarte de reads demasiado cortas
    CUTADAPT_TRIMMING(fastq_samples)

    trimmed_reads = CUTADAPT_TRIMMING.out.reads
    reports = reports.mix(CUTADAPT_TRIMMING.out.log.collect{ meta, log -> log })
    versions = versions.mix(CUTADAPT_TRIMMING.out.versions.first())

    // Filtrado por primer
    if (params.filter_by_primer) {
        CUTADAPT_FILTER_BY_PRIMER(trimmed_reads)

        reads_for_prinseq = CUTADAPT_FILTER_BY_PRIMER.out.reads
        reports = reports.mix(CUTADAPT_FILTER_BY_PRIMER.out.log.collect{ meta, log -> log })
        versions = versions.mix(CUTADAPT_FILTER_BY_PRIMER.out.versions.first())

    } else {
        reads_for_prinseq = trimmed_reads
    }

    // Trimeado por calidad de las reads
    PRINSEQ(reads_for_prinseq)

    preprocessed_reads = PRINSEQ.out.good_reads
    reports = reports.mix(PRINSEQ.out.log.collect{ meta, log -> log })
    versions = versions.mix(PRINSEQ.out.versions.first())

    // Preprocessed QC 1
    FASTQC_PREPROCESSED(preprocessed_reads)

    reports = reports.mix(FASTQC_PREPROCESSED.out.zip.collect{ meta, logs -> logs })
    versions = versions.mix(FASTQC_PREPROCESSED.out.versions.first())

    // TODO - preprocessed QC 2: fastq_screen. Mirar #3157 de nfcore/modules.
    // FASTQ_SCREEN()

    emit:
    preprocessed_reads // reads ya preprocesadas
    reports
    versions
}

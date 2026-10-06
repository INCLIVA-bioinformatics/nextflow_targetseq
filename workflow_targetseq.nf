/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    workflow_targetseq.nf

    Pipeline de rutina adaptado a nextflow para el analisis de datos centrado exclusivamente en los paneles de genes
    basados en captura mediante sondas o en amplificación por PCR. Intenta seguir la estructura de archivos y
    directorios utilizada en nf-core. También utiliza, cuando es posible, módulos desarrollados para éste. Además, toma
    algunas de las funciones desarrolladas para sarek.

    Versión: 1.0.0
    Fecha: 06/10/2026
    Autores: Unidad de Bioinformatica de INCLIVA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOG INICIAL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// Para obtener información para la fecha y el timestamp:
// import java.text.SimpleDateFormat
// def date = new Date()
// def sdf  = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss")

// params.pipeline_version = "1.0.0"

// // Parámetros para el log:
// params.date = sdf.format(date)

// // Información de log mostrada en pantalla:
// log.info """
//     T A R G E T S E Q  W O R K F L O W
//     ==================================
//     Panel             = ${params.panel}
//     Nextflow version  = ${nextflow.version}
//     Pipeline version  = ${params.pipeline_version}

//     Start timestamp   = ${params.date}
//     User name         = ${workflow.userName}

//     Project directory = ${projectDir}
//     Output directory  = ${params.outdir}
//     Work directory    = ${workflow.workDir}
//     ==================================
//     """
//     .stripIndent()

// // Información de ayuda:
// params.help = false

// if (params.help) { // Se muestra en pantalla cuando se introduce la flag "--help".
//     log.info "Texto de ayuda." // TODO - Incluir texto de ayuda
//     exit 1
// }

// workflow.onComplete {
//     def end_date = new Date()
//     params.end_date = sdf.format(end_date)
//     log.info "\nCompletion timestamp: ${params.end_date} " // Better than ${workflow.complete}
//     log.info(workflow.success ? "\nLa ejecución del pipeline ha terminado correctamente.\n" : "Oops, algo ha fallado.")
// }

// workflow.onError {
//     log.info "\nOops... La ejecución del pipeline se ha detenido con el siguiente mensaje: ${workflow.errorMessage}"
// }

// TODO - Incluir subworkflows utilizados
include { PREPROCESSING          } from './subworkflows/preprocessing'
include { MAPPING_POSTPROCESSING } from './subworkflows/mapping_postprocessing'
include { VARIANT_CALLING        } from './subworkflows/variant_calling'

workflow {

    if (params.help) {
        log.info '''
        TARGETSEQ

        Uso:
          nextflow run workflow_targetseq.nf -profile singularity \
            --samplesheet samplesheet.csv \
            --fasta reference.fa \
            --fasta_fai reference.fa.fai \
            --dict reference.dict \
            --bwamem2 /path/to/bwamem2/index \
            --intervals capture.bed \
            --target_intervals target.bed \
            --outdir results
        '''.stripIndent()

        return
    }

    log.info """
        T A R G E T S E Q  W O R K F L O W
        ==================================
        Panel             = ${params.panel}
        Nextflow version  = ${nextflow.version}
        Pipeline version  = 1.0.0

        User name         = ${workflow.userName}

        Project directory = ${projectDir}
        Output directory  = ${params.outdir}
        Work directory    = ${workflow.workDir}
        ==================================
    """.stripIndent()

    /*
     * Canales de recursos de referencia.
     *
     * collect() mantiene un único item con todos los ficheros de cada
     * recurso, compatible con los módulos que requieren listas de paths.
     */



    def fasta = params.fasta
        ? Channel.fromPath(params.fasta, checkIfExists: true).collect()
        : Channel.empty()

    def fasta_fai = params.fasta_fai
        ? Channel.fromPath(params.fasta_fai, checkIfExists: true).collect()
        : Channel.empty()

    def dict = params.dict
        ? Channel.fromPath(params.dict, checkIfExists: true).collect()
        : Channel.empty()

    def bwamem2 = params.bwamem2
        ? Channel.fromPath(params.bwamem2 + '*', checkIfExists: true).collect()
        : Channel.empty()

    def intervals = params.intervals
        ? Channel.fromPath(params.intervals, checkIfExists: true).collect()
        : Channel.empty()

    def target = params.target_intervals
        ? Channel.fromPath(params.target_intervals, checkIfExists: true)
        : Channel.empty()

    def dbsnp = params.dbsnp
        ? Channel.fromPath(params.dbsnp, checkIfExists: true)
        : Channel.empty()

    def gatk_indels = params.gatk_indels
        ? Channel.fromPath(params.gatk_indels, checkIfExists: true)
        : Channel.empty()

    def dbsnp_tbi = params.dbsnp && params.dbsnp_tbi
        ? Channel.fromPath(params.dbsnp_tbi, checkIfExists: true)
        : Channel.empty()

    def gatk_indels_tbi = params.gatk_indels && params.gatk_indels_tbi
        ? Channel.fromPath(params.gatk_indels_tbi, checkIfExists: true)
        : Channel.empty()

    def known_sites = dbsnp
        .concat(gatk_indels)
        .collect()

    def known_sites_tbi = dbsnp_tbi
        .concat(gatk_indels_tbi)
        .collect()

    def fasta_meta           = fasta.map { f -> [ [id:'fasta'], f ] }
    def fasta_fai_meta       = fasta_fai.map { f -> [ [id:'fai'], f ] }
    def dict_meta            = dict.map { f -> [ [id:'dict'], f ] }
    def known_sites_meta     = known_sites.map { f -> [ [id:'known_sites'], f ] }
    def known_sites_tbi_meta = known_sites_tbi.map { f -> [ [id:'known_sites_tbi'], f ] }
    /*
     * Canales acumuladores de reports y versiones.
     */

    def reports = Channel.empty()
    def versions = Channel.empty()

    /*
     * Samplesheet.
     *
     * extract_csv devuelve:
     *   [meta, fastq_1, fastq_2]
     */

    def samples = extract_csv(params.samplesheet)

    PREPROCESSING(
        samples,
        reports,
        versions
    )

    def preprocessed_reads = PREPROCESSING.out.preprocessed_reads
    reports = PREPROCESSING.out.reports
    versions = PREPROCESSING.out.versions

    MAPPING_POSTPROCESSING(
        preprocessed_reads,
        bwamem2,
        intervals,
        target,
        fasta,
        fasta_fai,
        dict,
        known_sites,
        known_sites_tbi,
        versions
    )

    def postprocessed_bam = MAPPING_POSTPROCESSING.out.postprocessed_bam
    def postprocessed_bai = MAPPING_POSTPROCESSING.out.postprocessed_bai
    def covered_regions = MAPPING_POSTPROCESSING.out.covered_regions

    versions = MAPPING_POSTPROCESSING.out.versions

    VARIANT_CALLING(
        postprocessed_bam,
        postprocessed_bai,
        fasta,
        fasta_fai,
        dict,
        covered_regions,
        versions
    )
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS (taken from nf-core/sarek)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
def extract_csv(csv_file) {

    // -----------------------------------------------------------------
    // 1. Validaciones iniciales: duplicados y consistencia patient/sample/lane
    // -----------------------------------------------------------------
    def patient_sample_lane_combinations = []
    def sample2patient = [:]

    def validation_ch = Channel.fromPath(csv_file, checkIfExists: true).splitCsv(header: true)

    validation_ch.map { row ->
        if (!row.lane) {
            error(
                'The sample sheet should specify a lane for patient "' +
                row.patient.toString() +
                '" and sample "' +
                row.sample.toString() +
                '".'
            )
        }

        def patient_sample_lane = [
            row.patient.toString(),
            row.sample.toString(),
            row.lane.toString()
        ]

        if (patient_sample_lane in patient_sample_lane_combinations) {
            error(
                'The patient-sample-lane combination "' +
                row.patient.toString() +
                '", "' +
                row.sample.toString() +
                '", and "' +
                row.lane.toString() +
                '" is present multiple times in the sample sheet.'
            )
        }

        patient_sample_lane_combinations.add(patient_sample_lane)

        if (!sample2patient.containsKey(row.sample.toString())) {
            sample2patient[row.sample.toString()] = row.patient.toString()
        }
        else if (sample2patient[row.sample.toString()] != row.patient.toString()) {
            error(
                'The sample "' +
                row.sample.toString() +
                '" is registered for both patient "' +
                row.patient.toString() +
                '" and "' +
                sample2patient[row.sample.toString()] +
                '" in the sample sheet.'
            )
        }

        return row
    }

    // -----------------------------------------------------------------
    // 2. Canal principal: agrupar por patient+sample, contar lanes,
    //    construir tuplas [meta, fastq_1, fastq_2] / [meta, bam, bai] / [meta, vcf]
    // -----------------------------------------------------------------
    def samples_ch = Channel.fromPath(csv_file, checkIfExists: true).splitCsv(header: true)

    def grouped_ch = samples_ch.map { row ->
        if (!(row.patient && row.sample)) {
            error(
                "Missing field in csv file header. " +
                "The csv file must have fields named 'patient' and 'sample'."
            )
        }

        if (row.patient.contains(' ') || row.sample.contains(' ')) {
            error(
                "Invalid value in csv file. " +
                "Values for 'patient' and 'sample' can not contain spaces."
            )
        }

        return [
            [row.patient.toString(), row.sample.toString()],
            row
        ]
    }

    def lanes_ch = grouped_ch
        .groupTuple()
        .map { sample_key, rows ->
            return [rows, rows.size()]
        }
        .transpose()

    def output_ch = lanes_ch.map { row, num_lanes ->

        def meta = [:]

        meta.patient = row.patient.toString()
        meta.sample  = row.sample.toString()
        meta.sex     = row.sex ? row.sex.toString() : 'NA'
        meta.status  = row.status ? row.status.toInteger() : 0
        meta.num_lanes = num_lanes.toInteger()
        meta.size = 1

        // -----------------------------------------------------------------
        // FASTQ
        // -----------------------------------------------------------------
        if (row.lane && row.fastq_1 && row.fastq_2) {
            meta.id = "${row.sample}-${row.lane}".toString()
            meta.data_type = 'fastq'

            def fastq_1 = file(row.fastq_1, checkIfExists: true)
            def fastq_2 = file(row.fastq_2, checkIfExists: true)

            def cn = params.seq_center ? "CN:${params.seq_center}\\t" : ''
            def flowcell = flowcellLaneFromFastq(fastq_1)

            meta.read_group =
                "\"@RG\\tID:${flowcell}.${row.sample}.${row.lane}\\t${cn}" +
                "PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}" +
                "\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

            return [meta, fastq_1, fastq_2]
        }

        // -----------------------------------------------------------------
        // BAM
        // -----------------------------------------------------------------
        if (row.lane && row.bam) {
            meta.id = "${row.sample}-${row.lane}".toString()
            meta.data_type = 'bam'

            def bam = file(row.bam, checkIfExists: true)
            def bai = row.bai ? file(row.bai, checkIfExists: true) : []

            return [meta, bam, bai]
        }

        // -----------------------------------------------------------------
        // VCF (solo para annotate)
        // -----------------------------------------------------------------
        if (row.vcf) {
            meta.id = meta.sample
            meta.data_type = 'vcf'
            meta.variantcaller = row.variantcaller ?: ''

            if (params.step != 'annotate') {
                error(
                    "Samplesheet contains VCF files but step is `${params.step}`. " +
                    "Use --step annotate or provide FASTQ/BAM input."
                )
            }

            def vcf = file(row.vcf, checkIfExists: true)
            return [meta, vcf]
        }

        // -----------------------------------------------------------------
        // Caso no reconocido
        // -----------------------------------------------------------------
        error(
            "Missing or unknown input fields in samplesheet. " +
            "Expected fastq_1/fastq_2, bam, or vcf."
        )
    }

    return output_ch
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}



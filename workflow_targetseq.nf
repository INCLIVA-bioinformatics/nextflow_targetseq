/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    workflow_targetseq.nf

    Pipeline de rutina adaptado a nextflow para el analisis de datos centrado exclusivamente en los paneles de genes
    basados en captura mediante sondas o en amplificación por PCR. Intenta seguir la estructura de archivos y
    directorios utilizada en nf-core. También utiliza, cuando es posible, módulos desarrollados para éste. Además, toma
    algunas de las funciones desarrolladas para sarek.

    Versión: 0.0.1
    Fecha: 21/08/2023
    Autores: Unidad de Bioinformatica de INCLIVA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOG INICIAL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// Para obtener información para la fecha y el timestamp:
import java.text.SimpleDateFormat
def date = new Date()
def sdf  = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss")

params.pipeline_version = "0.0.1"

// Parámetros para el log:
params.date = sdf.format(date)

// Información de log mostrada en pantalla:
log.info """
    T A R G E T S E Q  W O R K F L O W
    ==================================
    Panel             = ${params.panel}
    Nextflow version  = ${nextflow.version}
    Pipeline version  = ${params.pipeline_version}

    Start timestamp   = ${params.date}
    User name         = ${workflow.userName}

    Project directory = ${projectDir}
    Output directory  = ${params.outdir}
    Work directory    = ${workflow.workDir}
    ==================================
    """
    .stripIndent()

// Información de ayuda:
params.help = false

if (params.help) { // Se muestra en pantalla cuando se introduce la flag "--help".
    log.info "Texto de ayuda." // TODO - Incluir texto de ayuda
    exit 1
}

workflow.onComplete {
    def end_date = new Date()
    params.end_date = sdf.format(end_date)
    log.info "\nCompletion timestamp: ${params.end_date} " // Better than ${workflow.complete}
    log.info(workflow.success ? "\nLa ejecución del pipeline ha terminado correctamente.\n" : "Oops, algo ha fallado.")
}

workflow.onError {
    log.info "\nOops... La ejecución del pipeline se ha detenido con el siguiente mensaje: ${workflow.errorMessage}"
}

// TODO - Incluir subworkflows utilizados
include { PREPROCESSING          } from "${params.subworkflows_path}/preprocessing.nf"
include { MAPPING_POSTPROCESSING } from "${params.subworkflows_path}/mapping_postprocessing.nf"
include { VARIANT_CALLING        } from "${params.subworkflows_path}/variant_calling.nf"

workflow {
    // Initialize file channels based on params
    fasta       = params.fasta            ? Channel.fromPath(params.fasta).collect()            : null
    fasta_fai   = params.fasta_fai        ? Channel.fromPath(params.fasta_fai).collect()        : null
    dict        = params.dict             ? Channel.fromPath(params.dict).collect()             : null
    bwamem2     = params.bwamem2          ? Channel.fromPath(params.bwamem2).collect()          : null
    intervals   = params.intervals        ? Channel.fromPath(params.intervals).collect()        : null // capture.bed
    target      = params.target_intervals ? Channel.fromPath(params.target_intervals)           : null // target.bed
    dbsnp       = params.dbsnp            ? Channel.fromPath(params.dbsnp).collect()            : null
    gatk_indels = params.gatk_indels      ? Channel.fromPath(params.gatk_indels).collect()      : null

    // To gather all QC reports for MultiQC
    reports = Channel.empty()
    // To gather used softwares versions for MultiQC
    versions = Channel.empty()

    // Tabix indexed vcf files
    dbsnp_tbi       = params.dbsnp       ? params.dbsnp_tbi       ? Channel.fromPath(params.dbsnp_tbi).collect()       : null : Channel.value([])
    gatk_indels_tbi = params.gatk_indels ? params.gatk_indels_tbi ? Channel.fromPath(params.gatk_indels_tbi).collect() : null : Channel.value([])

    // known_sites is made by grouping both the dbsnp and the known indels (gatk_indels) resources
    known_sites     = dbsnp.concat(gatk_indels).collect()
    known_sites_tbi = dbsnp_tbi.concat(gatk_indels_tbi).collect()

    // Parsear muestras de la samplesheet y cargarlas a un canal
    csv_file = params.samplesheet
    samples  = extract_csv(file(csv_file))

    // TODO - Añadir procesado de UMIs

    PREPROCESSING(samples,
                  reports,
                  versions)

    preprocessed_reads = PREPROCESSING.out.preprocessed_reads
    reports = PREPROCESSING.out.reports
    versions = PREPROCESSING.out.versions

    MAPPING_POSTPROCESSING(preprocessed_reads,
                           bwamem2,
                           intervals,
                           target,
                           fasta,
                           fasta_fai,
                           dict,
                           known_sites,
                           known_sites_tbi,
                           versions)

    postprocessed_bam = MAPPING_POSTPROCESSING.out.postprocessed_bam
    postprocessed_bai = MAPPING_POSTPROCESSING.out.postprocessed_bai
    covered_regions = MAPPING_POSTPROCESSING.out.covered_regions

    VARIANT_CALLING(postprocessed_bam,
                    postprocessed_bai,
                    fasta,
                    fasta_fai,
                    dict,
                    covered_regions,
                    versions)
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS (taken from nf-core/sarek)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
def extract_csv(csv_file) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, samplesheet_line_count = 0;
        while ((line = reader.readLine()) != null) {samplesheet_line_count ++}
            if (samplesheet_line_count < 2) {
                error("Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines.")
            }
    }

    // Additional check of sample sheet:
    // 1. Each row should specify a lane and the same combination of patient, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different patients.
    def patient_sample_lane_combinations = []
    def sample2patient = [:]

    Channel.of(csv_file).splitCsv(header: true)
        .map{ row ->
            if ( !row.lane ) {  // This also handles the case where the lane is left as an empty string
                error('The sample sheet should specify a lane for patient "' + row.patient.toString() + '" and sample "' + row.sample.toString() + '".')
            }
            def patient_sample_lane = [row.patient.toString(), row.sample.toString(), row.lane.toString()]
            if (patient_sample_lane in patient_sample_lane_combinations) {
                error('The patient-sample-lane combination "' + row.patient.toString() + '", "' + row.sample.toString() + '", and "' + row.lane.toString() + '" is present multiple times in the sample sheet.')
            } else {
                patient_sample_lane_combinations.add(patient_sample_lane)
            }

            if (!sample2patient.containsKey(row.sample.toString())) {
                sample2patient[row.sample.toString()] = row.patient.toString()
            } else if (sample2patient[row.sample.toString()] != row.patient.toString()) {
                error('The sample "' + row.sample.toString() + '" is registered for both patient "' + row.patient.toString() + '" and "' + sample2patient[row.sample.toString()] + '" in the sample sheet.')
            }
        }

    sample_count_all = 0
    sample_count_normal = 0
    sample_count_tumor = 0

    Channel.of(csv_file).splitCsv(header: true)
        // Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            sample_count_all ++
            if (!(row.patient && row.sample)) {
                error("Missing field in csv file header. The csv file must have fields named 'patient' and 'sample'.")
            }
            else if (row.patient.contains(" ") || row.sample.contains(" ")) {
                error("Invalid value in csv file. Values for 'patient' and 'sample' can not contain space.")
            }
            [ [ row.patient.toString(), row.sample.toString() ], row ]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [ rows, size ]
        }.transpose()
        .map{ row, num_lanes -> // from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no sex specified, sex is not considered
        // sex is only mandatory for somatic CNV
        if (row.sex) meta.sex = row.sex.toString()
        else meta.sex = 'NA'

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        if (meta.status == 1) sample_count_tumor ++
        else sample_count_normal ++

        // start from FASTQ
        if (row.lane && row.fastq_2) { // TODO - añadir parseo de UMIs?
            meta.id         = "${row.sample}-${row.lane}".toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''

            def flowcell    = flowcellLaneFromFastq(fastq_1)
            // Don't use a random element for ID, it breaks resuming. read_group "@RG" is needed for aligner.
            def read_group  = "\"@RG\\tID:${flowcell}.${row.sample}.${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

            meta.num_lanes  = num_lanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'fastq'

            meta.size       = 1 // default number of splitted fastq

            return [ meta, fastq_1, fastq_2 ]
        }

        // start from BAM
        if (row.lane && row.bam) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def bai         = row.bai ? file(row.bai,   checkIfExists: true) : []

            meta.num_lanes  = num_lanes.toInteger()
            meta.data_type  = 'bam'

            meta.size       = 1 // default number of splitted fastq

            return [ meta, bam, bai ]
        }

        // start from VCF (annotation)
        else if (row.vcf) {
            meta.id = meta.sample
            def vcf = file(row.vcf, checkIfExists: true)

            meta.data_type     = 'vcf'
            meta.variantcaller = row.variantcaller ?: ''

            if (params.step == 'annotate') return [ meta, vcf ]
            else {
                error("Samplesheet contains vcf files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }
        } else {
            error("Missing or unknown field in csv file header. Please check your samplesheet")
        }
    }
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

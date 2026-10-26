include { INC_VEP                                                } from '../modules/incliva/annotation/main'
include { INC_ANNOTATE_HOTSPOTS                                  } from '../modules/incliva/annotation/main'
include { INC_SNPSIFT_ANNOTATE as SNPSIFT_COSMIC                 } from '../modules/incliva/annotation/main'
include { INC_MV_COSMIC_ID                                       } from '../modules/incliva/annotation/main'
include { INC_SNPSIFT_ANNOTATE as SNPSIFT_MIXED                  } from '../modules/incliva/annotation/main'
include { INC_SNPSIFT_ANNOTATE as SNPSIFT_INC_DB                 } from '../modules/incliva/annotation/main'
include { INC_ADD_CONSERVATION                                   } from '../modules/incliva/annotation/main'
include { INC_SNPSIFT_ANNOTATE as SNPSIFT_ONCOKB                 } from '../modules/incliva/annotation/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_ANNOTATED         } from '../modules/nf-core/tabix/bgziptabix/main'
include { VCFLIB_VCFUNIQ as VCFLIB_VCFUNIQ_ANNOTATED             } from '../modules/incliva/vcflib/vcfuniq/main'
include { INC_VCF_TO_CSV                                         } from '../modules/incliva/annotation/main'
include { INC_AA_3TO1                                            } from '../modules/incliva/annotation/main'
include { INC_SAMPLES_CONTAMINATION                              } from '../modules/incliva/annotation/main'
include { INC_SPLIT_ANNOTATED_CSV                                } from '../modules/incliva/annotation/main'

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATION: anotación del vcf multi-muestra y csv de variantes (rutina.sh: annotate_variants)
        VEP + SpliceAI → cancerhotspots → COSMIC → mv_cosmic_id → mixed_db → INCLIVA somatic → conservación → OncoKB
        → vcfuniq → csv anotado (vcf_to_csv + aa_3to1) → contaminación entre muestras y csv por muestra
    (pangolin no se incluye: rutina.sh sólo lo ejecuta con pangolin_splicing="yes" y ningún panel lo usa)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
workflow ANNOTATION {
    take:
    run_vcf  // channel: [meta, vcf.gz] vcf multi-muestra
    fasta
    versions

    main:
    def vep_cache       = file(params.vep_cache, checkIfExists: true)
    def spliceai_snv    = [ file(params.spliceai_snv, checkIfExists: true), file("${params.spliceai_snv}.tbi", checkIfExists: true) ]
    def spliceai_indel  = [ file(params.spliceai_indel, checkIfExists: true), file("${params.spliceai_indel}.tbi", checkIfExists: true) ]
    def hotspots        = [ file(params.hotspots, checkIfExists: true), file("${params.hotspots}.tbi", checkIfExists: true) ]
    def cosmic          = [ file(params.cosmic, checkIfExists: true), file("${params.cosmic}.tbi", checkIfExists: true) ]
    def mixed_db        = [ file(params.mixed_db, checkIfExists: true), file("${params.mixed_db}.tbi", checkIfExists: true) ]
    def inc_db          = [ file(params.inc_db, checkIfExists: true), [] ]
    def oncokb          = [ file(params.oncokb, checkIfExists: true), file("${params.oncokb}.tbi", checkIfExists: true) ]
    def conservation    = file("${params.panel_design}/conservation.bed", checkIfExists: true)
    // rutina.sh: -r ${panel}/list_of_transcripts si el panel tiene lista de tránscritos
    def transcripts     = file("${params.panel_design}/list_of_transcripts")
    def list_of_transcripts = transcripts.exists() ? transcripts : []

    INC_VEP(run_vcf, fasta, vep_cache, spliceai_snv, spliceai_indel)
    versions = versions.mix(INC_VEP.out.versions)

    INC_ANNOTATE_HOTSPOTS(INC_VEP.out.vcf, hotspots)
    versions = versions.mix(INC_ANNOTATE_HOTSPOTS.out.versions)

    SNPSIFT_COSMIC(INC_ANNOTATE_HOTSPOTS.out.vcf, cosmic)
    versions = versions.mix(SNPSIFT_COSMIC.out.versions)

    INC_MV_COSMIC_ID(SNPSIFT_COSMIC.out.vcf)
    versions = versions.mix(INC_MV_COSMIC_ID.out.versions)

    SNPSIFT_MIXED(INC_MV_COSMIC_ID.out.vcf, mixed_db)
    versions = versions.mix(SNPSIFT_MIXED.out.versions)

    SNPSIFT_INC_DB(SNPSIFT_MIXED.out.vcf, inc_db)
    versions = versions.mix(SNPSIFT_INC_DB.out.versions)

    INC_ADD_CONSERVATION(SNPSIFT_INC_DB.out.vcf, conservation)
    versions = versions.mix(INC_ADD_CONSERVATION.out.versions)

    SNPSIFT_ONCOKB(INC_ADD_CONSERVATION.out.vcf, oncokb)
    versions = versions.mix(SNPSIFT_ONCOKB.out.versions)

    TABIX_BGZIPTABIX_ANNOTATED(SNPSIFT_ONCOKB.out.vcf)
    VCFLIB_VCFUNIQ_ANNOTATED(TABIX_BGZIPTABIX_ANNOTATED.out.gz_index)
    versions = versions.mix(VCFLIB_VCFUNIQ_ANNOTATED.out.versions)

    annotated_vcf = VCFLIB_VCFUNIQ_ANNOTATED.out.vcf

    INC_VCF_TO_CSV(annotated_vcf, list_of_transcripts)
    versions = versions.mix(INC_VCF_TO_CSV.out.versions)

    INC_AA_3TO1(INC_VCF_TO_CSV.out.csv)
    versions = versions.mix(INC_AA_3TO1.out.versions)

    annotated_csv = INC_AA_3TO1.out.csv

    INC_SAMPLES_CONTAMINATION(annotated_csv)
    versions = versions.mix(INC_SAMPLES_CONTAMINATION.out.versions)

    INC_SPLIT_ANNOTATED_CSV(annotated_csv)
    versions = versions.mix(INC_SPLIT_ANNOTATED_CSV.out.versions)

    emit:
    vcf = annotated_vcf
    csv = annotated_csv
    versions
}

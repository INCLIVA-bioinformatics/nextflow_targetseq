include { GATK4_MUTECT2                                          } from '../modules/nf-core/gatk4/mutect2/main'
include { GATK4_FILTERMUTECTCALLS                                } from '../modules/nf-core/gatk4/filtermutectcalls/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MUTECT2                 } from '../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_LOFREQ                  } from '../modules/nf-core/bcftools/norm/main'
include { INC_SELECT_VARIANTS_GATK                               } from '../incliva_modules'
include { INC_SELECT_VARIANTS_LOFREQ                             } from '../incliva_modules'
include { INC_ANNOTATE_BAM_INFO as INC_ANNOTATE_BAM_INFO_MUTECT2 } from '../incliva_modules'
include { INC_ANNOTATE_BAM_INFO as INC_ANNOTATE_BAM_INFO_LOFREQ  } from '../incliva_modules'
include { VCFLIB_VCFUNIQ as VCFLIB_VCFUNIQ_MUTECT2               } from '../modules/incliva/vcflib/vcfuniq/main'
include { VCFLIB_VCFUNIQ as VCFLIB_VCFUNIQ_LOFREQ                } from '../modules/incliva/vcflib/vcfuniq/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_MUTECT2                 } from '../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_LOFREQ                  } from '../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_MUTECT2           } from '../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_LOFREQ            } from '../modules/nf-core/tabix/bgziptabix/main'
//include { TABIX_BGZIPTABIX as TABIX_BGZIP                                            } from "/mnt/zonahpc/home/bioinformatica/comun/software/sarek_fork_INCLIVA/nf-core/tabix/bgziptabix/main"
include { TABIX_TABIX                                            } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as INDEX_MUTECT2_VCF                       } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as INDEX_LOFREQ_VCF                        } from '../modules/nf-core/tabix/tabix/main'
include { LOFREQ_CALLPARALLEL                                    } from '../modules/nf-core/lofreq/callparallel/main'
include { INC_ANNOTATE_FORMAT_LOFREQ                             } from '../incliva_modules'
include { BCFTOOLS_CONCAT                                        } from '../modules/nf-core/bcftools/concat/main'


workflow VARIANT_CALLING {
    take:
    postprocessed_bam
    postprocessed_bai
    fasta
    fasta_fai
    dict
    covered_regions
    versions

    main:
    // Se preparan canales para los callers
    bam_vcf = postprocessed_bam.join(postprocessed_bai).join(covered_regions)
    fasta_vc = fasta.map { f ->
        [ [ id: f.baseName ], f ]
    }
    
    fasta_fai_vc = fasta_fai.map { fai ->
    [ [ id: fai.baseName ], fai ]
    }

    dict_vc = dict.map { d ->
        [ [ id: d.baseName ], d ]
    }
    germline_resource = Channel.value([]) // channel required for mutect2
    germline_resource_tbi = Channel.value([]) // channel required for mutect2
    panel_of_normals = Channel.value([]) // channel required for mutect2
    panel_of_normals_tbi = Channel.value([]) // channel required for mutect2

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PASO 3-A: VARIANT CALLING (MUTECT2) // TODO - Crear subworkflows para cada uno de los callers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    GATK4_MUTECT2(bam_vcf, fasta_vc, fasta_fai_vc, dict_vc, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)

    raw_mutect2_vcf = GATK4_MUTECT2.out.vcf
    raw_mutect2_tbi = GATK4_MUTECT2.out.tbi
    raw_mutect2_stats = GATK4_MUTECT2.out.stats
    versions = versions.mix(GATK4_MUTECT2.out.versions.first())

    // Se utiliza FilterMutectCalls siguiendo la guía de buenas prácticas
    mutect2_vcf_to_filter = raw_mutect2_vcf.join(raw_mutect2_tbi, failOnMismatch: true).join(raw_mutect2_stats, failOnMismatch: true).map{ meta, vcf, tbi, stats -> [ meta + [ data_type:'vcf' ], vcf, tbi, stats, [], [], [], [] ] } // Actualizar meta.data_type ya que ahora tenemos vcf. Los 4 campos vacíos corresponden a oriantion_bias, segmentation, table, estimate (necesarios para el módulo GATK4_FILTERMUTECTCALLS)

    GATK4_FILTERMUTECTCALLS(mutect2_vcf_to_filter, fasta_vc, fasta_fai_vc, dict_vc)

    mutect2_vcf = GATK4_FILTERMUTECTCALLS.out.vcf
    mutect2_tbi = GATK4_FILTERMUTECTCALLS.out.tbi
    versions = versions.mix(GATK4_FILTERMUTECTCALLS.out.versions.first())

    // Se normalizan las variantes (indels)
    mutect2_variants = mutect2_vcf.join(mutect2_tbi, failOnMismatch: true)

    BCFTOOLS_NORM_MUTECT2(mutect2_variants, fasta_vc)

    normalized_mutect2_variants = BCFTOOLS_NORM_MUTECT2.out.vcf
    versions = versions.mix(BCFTOOLS_NORM_MUTECT2.out.versions.first())

    // Filtrar variantes por AD y AF y se añade "BAM = Clipped" INFO field
    INC_SELECT_VARIANTS_GATK(normalized_mutect2_variants)
    INC_ANNOTATE_BAM_INFO_MUTECT2(INC_SELECT_VARIANTS_GATK.out.selected_variants)

    versions = versions.mix(INC_ANNOTATE_BAM_INFO_MUTECT2.out.versions.first())

    TABIX_BGZIPTABIX_MUTECT2(INC_ANNOTATE_BAM_INFO_MUTECT2.out.bam_info)

    selected_mutect2_variants = TABIX_BGZIPTABIX_MUTECT2.out.gz_index
    //versions = versions.mix(TABIX_BGZIPTABIX_MUTECT2.out.versions.first())
    versions = versions
        .mix(TABIX_BGZIPTABIX_MUTECT2.out.versions_tabix)
        .mix(TABIX_BGZIPTABIX_MUTECT2.out.versions_bgzip)

    // Se eliminan variantes duplicadas y se ordena el vcf resultante final de Mutect2
    VCFLIB_VCFUNIQ_MUTECT2(selected_mutect2_variants)
    versions = versions.mix(VCFLIB_VCFUNIQ_MUTECT2.out.versions.first())

    BCFTOOLS_SORT_MUTECT2(VCFLIB_VCFUNIQ_MUTECT2.out.vcf)

    sorted_mutect2_vcf = BCFTOOLS_SORT_MUTECT2.out.vcf
    //versions = versions.mix(BCFTOOLS_SORT_MUTECT2.out.versions.first())

    INDEX_MUTECT2_VCF(sorted_mutect2_vcf)

    mutect2_vcf = sorted_mutect2_vcf.join(INDEX_MUTECT2_VCF.out.tbi)

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PASO 3-B: VARIANT CALLING (LOFREQ) // TODO - Crear subworkflows para cada uno de los callers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    def fasta_meta = fasta.map { f -> [ [id: 'genome'], f ] }
    def fai_meta   = fasta_fai.map { f -> [ [id: 'genome'], f ] }

    LOFREQ_CALLPARALLEL(bam_vcf, fasta_meta, fai_meta)
    //LOFREQ_CALLPARALLEL(bam_vcf, fasta, fasta_fai)

    raw_lofreq_vcf = LOFREQ_CALLPARALLEL.out.vcf.map{ meta, vcf -> [ meta + [ data_type:'vcf' ], vcf] } // Actualizar meta.data_type ya que ahora tenemos vcf.
    versions = versions.mix(LOFREQ_CALLPARALLEL.out.versions.first())

    // Se normalizan las variantes
    TABIX_TABIX(raw_lofreq_vcf)

    BCFTOOLS_NORM_LOFREQ(raw_lofreq_vcf.join(TABIX_TABIX.out.tbi, failOnMismatch: true), fasta_vc)

    normalized_lofreq_variants = BCFTOOLS_NORM_LOFREQ.out.vcf

    // Filtrar variantes por AD y AF y se añade FORMAT y "BAM = Clipped" en INFO field
    INC_SELECT_VARIANTS_LOFREQ(normalized_lofreq_variants)
    INC_ANNOTATE_FORMAT_LOFREQ(INC_SELECT_VARIANTS_LOFREQ.out.selected_variants)
    INC_ANNOTATE_BAM_INFO_LOFREQ(INC_ANNOTATE_FORMAT_LOFREQ.out.formated_variants)

    TABIX_BGZIPTABIX_LOFREQ(INC_ANNOTATE_BAM_INFO_LOFREQ.out.bam_info)

    selected_lofreq_variants = TABIX_BGZIPTABIX_LOFREQ.out.gz_index
    versions = versions
        .mix(TABIX_BGZIPTABIX_LOFREQ.out.versions_tabix)
        .mix(TABIX_BGZIPTABIX_LOFREQ.out.versions_bgzip)

    // Se eliminan variantes duplicadas y se ordena el vcf resultante final de Lofreq
    VCFLIB_VCFUNIQ_LOFREQ(selected_lofreq_variants)
    versions = versions.mix(VCFLIB_VCFUNIQ_LOFREQ.out.versions.first())

    BCFTOOLS_SORT_LOFREQ(VCFLIB_VCFUNIQ_LOFREQ.out.vcf)

    sorted_lofreq_vcf = BCFTOOLS_SORT_LOFREQ.out.vcf
    versions = versions.mix(BCFTOOLS_SORT_LOFREQ.out.versions.first())

    INDEX_LOFREQ_VCF(sorted_lofreq_vcf)

    lofreq_vcf = sorted_lofreq_vcf.join(INDEX_LOFREQ_VCF.out.tbi)

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PASO 4: Juntar VCFs creados // TODO - acabar esta parte de amplicones. Sacar bedtools del script de python. // TODO - Mirar sarek github issue #738
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if (params.amplicon) {
        // TODO - python inc_gatk_realiable Low_confidence
        // TODO - python inc_gatk_realiable High_confidence
        calls_ch = mutect2_vcf.join(lofreq_vcf, failOnMismatch: true).map { m, vcf_m, tbi_m, vcf_l, tbi_l ->
        [ m, [ vcf_m, vcf_l ], [ tbi_m, tbi_l ] ]
    }
    
        BCFTOOLS_CONCAT(calls_ch)

        isected_vcf = BCFTOOLS_CONCAT.out.vcf
        versions = versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    } else {
    calls_ch = mutect2_vcf.join(lofreq_vcf, failOnMismatch: true).map { m, vcf_m, tbi_m, vcf_l, tbi_l ->
            [ m, [ vcf_m, vcf_l ], [ tbi_m, tbi_l ] ]
        }
        
        BCFTOOLS_CONCAT(calls_ch)

        isected_vcf = BCFTOOLS_CONCAT.out.vcf
        versions = versions.mix(BCFTOOLS_CONCAT.out.versions.first())
    }

    // python inc_reformat_genotype
    // python inc_check_empty_vcf

    emit:
    vcf = isected_vcf // VCF combinado por muestra (Mutect2 + LoFreq)
    versions
}

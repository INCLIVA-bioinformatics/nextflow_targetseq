include { GATK4_MUTECT2                                          } from '../modules/nf-core/gatk4/mutect2/main'
include { GATK4_FILTERMUTECTCALLS                                } from '../modules/nf-core/gatk4/filtermutectcalls/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MUTECT2                 } from '../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_LOFREQ                  } from '../modules/nf-core/bcftools/norm/main'
include { INC_RESCUE_VARIANTS                                    } from '../incliva_modules'
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
include { INC_REFORMAT_GENOTYPE                                  } from '../incliva_modules'
include { INC_TAG_RELIABLE as INC_TAG_RELIABLE_MUTECT2           } from '../incliva_modules'
include { INC_TAG_RELIABLE as INC_TAG_RELIABLE_LOFREQ            } from '../incliva_modules'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_TAG_MUTECT2       } from '../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_TAG_LOFREQ        } from '../modules/nf-core/tabix/bgziptabix/main'
include { INC_ADD_FMC_TAG                                        } from '../incliva_modules'
include { INC_ADD_CSB_INFO                                       } from '../incliva_modules'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_FINAL                   } from '../modules/nf-core/bcftools/sort/main'
include { INC_BCFTOOLS_MERGE                                     } from '../incliva_modules'


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

    // Se guarda el FILTER original en FORMAT (FMC) y se rescatan las variantes con la etiqueta "haplotype"
    INC_RESCUE_VARIANTS(normalized_mutect2_variants)
    versions = versions.mix(INC_RESCUE_VARIANTS.out.versions.first())

    // Filtrar variantes por AD y AF y se añade "BAM = Clipped" INFO field
    INC_SELECT_VARIANTS_GATK(INC_RESCUE_VARIANTS.out.rescued_variants)
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
    PASO 4: Juntar VCFs creados // TODO - Mirar sarek github issue #738
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if (params.amplicon) {
        // Variantes en zona de primers: Low_confidence para Mutect2 y High_confidence para Lofreq (REL en INFO)
        def primers = file("${params.panel_design}/primers.bed", checkIfExists: true)

        INC_TAG_RELIABLE_MUTECT2(sorted_mutect2_vcf, primers)
        INC_TAG_RELIABLE_LOFREQ(sorted_lofreq_vcf, primers)
        versions = versions
            .mix(INC_TAG_RELIABLE_MUTECT2.out.versions.first())
            .mix(INC_TAG_RELIABLE_LOFREQ.out.versions.first())

        TABIX_BGZIPTABIX_TAG_MUTECT2(INC_TAG_RELIABLE_MUTECT2.out.tagged_variants)
        TABIX_BGZIPTABIX_TAG_LOFREQ(INC_TAG_RELIABLE_LOFREQ.out.tagged_variants)

        mutect2_vcf = TABIX_BGZIPTABIX_TAG_MUTECT2.out.gz_index
        lofreq_vcf = TABIX_BGZIPTABIX_TAG_LOFREQ.out.gz_index
    }

    calls_ch = mutect2_vcf.join(lofreq_vcf, failOnMismatch: true).map { m, vcf_m, tbi_m, vcf_l, tbi_l ->
        [ m, [ vcf_m, vcf_l ], [ tbi_m, tbi_l ] ]
    }

    BCFTOOLS_CONCAT(calls_ch)
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    // Se recalculan GT y AF a partir de AD y se filtra por AF y profundidad mínima
    INC_REFORMAT_GENOTYPE(BCFTOOLS_CONCAT.out.vcf)
    versions = versions.mix(INC_REFORMAT_GENOTYPE.out.versions.first())

    // python inc_check_empty_vcf: no se incluye. Sólo actúa si el vcf no tiene ninguna línea, y tras el concat
    // siempre queda al menos la cabecera.

    // Se añaden las etiquetas originales de FilterMutectCalls (FMC, -1 si no la llamó Mutect2)
    INC_ADD_FMC_TAG(INC_REFORMAT_GENOTYPE.out.reformatted_variants.join(INC_RESCUE_VARIANTS.out.rescued_variants, failOnMismatch: true))
    versions = versions.mix(INC_ADD_FMC_TAG.out.versions.first())

    // Se añade el strand bias de cada caller (cSB)
    csb_input = INC_ADD_FMC_TAG.out.fmc_variants
        .join(mutect2_vcf.map { meta, vcf, tbi -> [ meta, vcf ] }, failOnMismatch: true)
        .join(lofreq_vcf.map { meta, vcf, tbi -> [ meta, vcf ] }, failOnMismatch: true)

    INC_ADD_CSB_INFO(csb_input)
    versions = versions.mix(INC_ADD_CSB_INFO.out.versions.first())

    BCFTOOLS_SORT_FINAL(INC_ADD_CSB_INFO.out.csb_variants)

    sample_vcf = BCFTOOLS_SORT_FINAL.out.vcf
    versions = versions.mix(BCFTOOLS_SORT_FINAL.out.versions.first())

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PASO 5: VCF multi-muestra (rutina.sh: bcftools merge -m none; con una sola muestra se usa su vcf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // Las muestras se unen en el orden del samplesheet (rutina.sh: orden de "samples" en el config)
    def sample_order = file(params.samplesheet).splitCsv(header: true)*.sample.unique()

    all_vcfs = sample_vcf
        .join(BCFTOOLS_SORT_FINAL.out.tbi, failOnMismatch: true)
        .toSortedList { a, b -> sample_order.indexOf(a[0].sample) <=> sample_order.indexOf(b[0].sample) }
        .map { files -> [ [ id: params.study_name ], files*.get(1), files*.get(2) ] }
        .branch { meta, vcfs, tbis ->
            multi:  vcfs.size() > 1
            single: vcfs.size() == 1
        }

    INC_BCFTOOLS_MERGE(all_vcfs.multi)
    versions = versions.mix(INC_BCFTOOLS_MERGE.out.versions)

    run_vcf = INC_BCFTOOLS_MERGE.out.vcf
        .mix(all_vcfs.single.map { meta, vcfs, tbis -> [ meta, vcfs[0] ] })

    emit:
    vcf     = sample_vcf // VCF combinado por muestra (Mutect2 + LoFreq)
    run_vcf              // VCF multi-muestra de la carrera
    versions
}

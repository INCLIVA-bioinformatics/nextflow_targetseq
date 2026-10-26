process GET_COVERED_BED {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/bedtools:2.30.0----he941832_1' :
    //     'biocontainers/bedtools:2.30.0----he941832_1' }"
    //container 'biocontainers/bedtools:2.30.0----he941832_1'

    input:
    tuple val(meta), path(depth_tsv)

    output:
    tuple val(meta), path("*_covered_regions_${params.cov}x.bed"), emit: covered_regions
    tuple val(meta), path("*_non_covered_${params.cov}x.bed"), emit: non_covered

    path "versions.yml"                                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '{if (\$3 >= '${params.cov}') print \$1"\t"\$2-1"\t"\$2}' ${depth_tsv} | bedtools merge ${args} > ${prefix}_covered_regions_${params.cov}x.bed
    bedtools subtract -a ${params.target_intervals} -b ${prefix}_covered_regions_${params.cov}x.bed > ${prefix}_non_covered_${params.cov}x.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_covered_regions_${params.cov}x.bed
    touch ${prefix}_non_covered_${params.cov}x.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}

process INC_SELECT_VARIANTS_GATK {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyvcf3=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.3--pyhdfd78af_0' :
        'biocontainers/pyvcf3:1.0.3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_inc_filtered.vcf"), emit: selected_variants
    path "versions.yml"                        , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    '''
    inc_select_variants_gatk.py

    Selecciona aquellas variantes cuya AF y profundidad de lectura es mayor que unas definidas en el config.
    Se utiliza con las variantes encontradas por gatk.

    Versión: 1.0
    Autor: UB INCLIVA
    Fecha: 11/12/2018
    '''
    import vcf
    import sys
    import os

    vcf_reader = vcf.Reader(open("${vcf}",'rb')) # Abrir binary mode.
    vcf_writer = vcf.Writer(open("${prefix}_inc_filtered.vcf",'w'), template = vcf_reader)

    # Como inc_select_variants_v1.4.py (hotfix 29/11/2024): no se filtra por profundidad, solo se avisa.
    # Tras normalizar, la suma de AD de un alelo puede quedar < cov aunque la posición esté cubierta a cov.
    for record in vcf_reader:
        for sample_index, sample in enumerate(record.samples):

            if sample.called:
                sample_vaf = round(float(sample.data.AD[1]) / (float(sample.data.AD[0]) + float(sample.data.AD[1])), 4)
                total_depth = sample.data.AD[1] + sample.data.AD[0]

                if sample_vaf >= ${params.af}:
                    vcf_writer.write_record(record)

                if total_depth < ${params.cov}:
                    print('WARNING: la variante %s no supera el corte de cobertura (%s).' % (record, ${params.cov}))

    # Versions
    command_string = '''cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    '''
    os.system(command_string)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_inc_filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python v//g")
    END_VERSIONS
    """
}

process INC_ANNOTATE_BAM_INFO {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyvcf3=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.3--pyhdfd78af_0' :
        'biocontainers/pyvcf3:1.0.3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_bam_info.vcf"), emit: bam_info
    path "versions.yml"                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    '''
    inc_anno_bam_info.py

    Añade a un archivo vcf la información sobre el bam utilizado para conseguir este archivo vcf. Añade tanto la línea
    en la cabecera como el campo correspondiente en la entrada de cada variante.

    Versión: 1.0
    Autor: UB INCLIVA
    Fecha: 07/01/2019
    '''
    import gzip
    import os

    with open("${vcf}", 'rt') as file:
        f = file.readlines()

    with open("${prefix}_bam_info.vcf", 'w') as outfile:
        sample = "${vcf}"
        sample = sample.split("/")[-1]
        sample = sample.split(".")[0]

        bam = "Clipped"

        for i in range(len(f)):
            line = f[i].rstrip()

            if line.startswith("#"):

                if line.startswith("##INFO") and not f[i+1].startswith("##INFO"):
                    outfile.write(f'{line}\\n')
                    outfile.write(f'##INFO=<ID=BAM,Number=1,Type=String,Description="BAM file used to get VCF">\\n')

                else:
                    outfile.write(f'{line}\\n')

            else:
                values = line.split("\\t")

                info = values[7].split(";")
                info = ";".join(info) + ";BAM=" + bam

                line = (f'{values[0]}\\t{values[1]}\\t{values[2]}\\t{values[3]}\\t{values[4]}\\t{values[5]}\\t{values[6]}\\t{info}\\t{values[8]}\\t{values[9]}\\n')

                outfile.write(line)

    # Versions
    command_string = '''cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    '''
    os.system(command_string)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_bam_info.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python v//g")
    END_VERSIONS
    """
}

process INC_SELECT_VARIANTS_LOFREQ {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyvcf3=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.3--pyhdfd78af_0' :
        'biocontainers/pyvcf3:1.0.3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_inc_filtered.vcf"), emit: selected_variants
    path "versions.yml"                        , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    '''
    inc_select_variants_lofreq.py

    Selecciona aquellas variantes cuya AF y profundidad de lectura es mayor que unas definidas en el config.
    Se utiliza con las variantes encontradas por lofreq.

    Versión: 1.0
    Autor: UB INCLIVA
    Fecha: 11/12/2018
    '''
    import vcf
    import sys
    import os

    vcf_reader = vcf.Reader(open("${vcf}",'rb')) # Abrir binary mode.
    vcf_writer = vcf.Writer(open("${prefix}_inc_filtered.vcf",'w'), template = vcf_reader)

    # Como inc_select_variants_v1.4.py (hotfix 29/11/2024): no se filtra por profundidad, solo se avisa.
    for record in vcf_reader:
        sample_vaf = record.INFO['AF']
        total_depth = record.INFO['DP']

        if sample_vaf >= ${params.af}:
            vcf_writer.write_record(record)

        if total_depth < ${params.cov}:
            print('WARNING: la variante %s no supera el corte de cobertura (%s).' % (record, ${params.cov}))

    # Versions
    command_string = '''cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    '''
    os.system(command_string)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_inc_filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}

process INC_ANNOTATE_FORMAT_LOFREQ {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyvcf3=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.3--pyhdfd78af_0' :
        'biocontainers/pyvcf3:1.0.3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_inc_formated.vcf"), emit: formated_variants
    path "versions.yml"                                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    '''
    inc_anno_format_lofreq.py

    Añade información a un vcf sobre campos creados por lofreq. Concretamente, añade información sobre AD, AF, DP y GT
    a la cabecera del archivo vcf. Tamnién añade esta información correctamente a las diferentes entradas de las
    variantes encontradas.

    Versión: 1.0
    Autor: UB INCLIVA
    Fecha: 09/05/2022
    '''
    import os

    with open("${vcf}", 'rt') as file:
        f = file.readlines()

    with open("${prefix}_inc_formated.vcf", 'w') as outfile:
        sample = "${meta.sample}"
        sample = sample.split("/")[-1]
        sample = sample.split(".")[0]

        for i in range(len(f)):
            line = f[i].rstrip()

            if line.startswith("##"):

                if line.startswith("##INFO") and f[i+1].startswith("##FILTER"):
                    outfile.write(line + "\\n")
                    outfile.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\\n')
                    outfile.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">\\n')
                    outfile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\\n')
                    outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\\n')

                else:
                    outfile.write(line + "\\n")

            elif line.startswith("#"):
                outfile.write(line + "\\t" + "FORMAT" + "\\t" + sample + "\\n")

            else:
                values = line.split("\\t")

                info = values[7].split(";")
                dp4 = info[3].split("=")[1]
                dps = dp4.split(",")
                ref = int(dps[0]) + int(dps[1])
                alt = int(dps[2]) + int(dps[3])
                ad = str(ref) + "," + str(alt)
                af = info[1].split("=")[1]
                dp = info[0].split("=")[1]

                outfile.write(line + "\\t" + "GT:AD:AF:DP" + "\\t" + ".:" + ad + ":" + af + ":" + dp + "\\n")

    # Versions
    command_string = '''cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    '''
    os.system(command_string)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_inc_formated.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}

process INC_RESCUE_VARIANTS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyvcf3=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.3--pyhdfd78af_0' :
        'biocontainers/pyvcf3:1.0.3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_rescued.vcf.gz"), emit: rescued_variants
    path "versions.yml"                      , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def labels_to_rescue = task.ext.args ?: ''
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    '''
    inc_rescue_variants_v1.0.py

    Añade la información del filtro añadido por FilterMutectCalls al campo de FORMAT (FMC). Además, si una variante
    ha sido descartada por FilterMutectCalls por alguna de las etiquetas indicadas en ext.args (separadas por comas),
    se quita esa etiqueta de FILTER (si era la única, pasa a "PASS").

    Correspondencia código-etiqueta en FMC: 0.PASS, 1.base_qual, 2.clustered_events, 3.contamination, 4.duplicate,
    5.fragment, 6.germline, 7.haplotype, 8.low_allele_frac, 9.map_qual, 10.multiallelic, 11.n_ratio,
    12.normal_artifact, 13.orientation, 14.panel_of_normals, 15.position, 16.possible_numt, 17.slippage,
    18.strand_bias, 19.strict_strand, 20.weak_evidence. El código -1 (no llamada por Mutect2) se asigna más adelante.

    Versión: 1.0
    Autor: UB INCLIVA - Francisco Martínez Picó, Roxana Andreea Moldovan
    Fecha: 10/10/2023
    '''
    import gzip
    import os

    labels_to_rescue = "${labels_to_rescue}"

    # Diccionario con la correspondencia etiqueta - código que se asigna en FMC
    codigo_etiqueta = {
        'PASS': 0,
        'base_qual': 1,
        'clustered_events': 2,
        'contamination': 3,
        'duplicate': 4,
        'fragment': 5,
        'germline': 6,
        'haplotype': 7,
        'low_allele_frac': 8,
        'map_qual': 9,
        'multiallelic': 10,
        'n_ratio': 11,
        'normal_artifact': 12,
        'orientation': 13,
        'panel_of_normals': 14,
        'position': 15,
        'possible_numt': 16,
        'slippage': 17,
        'strand_bias': 18,
        'strict_strand': 19,
        'weak_evidence': 20
    }

    with gzip.open("${vcf}", 'rt') as input_vcf, gzip.open("${prefix}_rescued.vcf.gz", 'wt') as output_vcf:
        for line in input_vcf:
            line = line.strip()

            if line.startswith('##'): # líneas de la cabecera en general
                output_vcf.write(line + '\\n')

            elif line.startswith('#'): # línea con las cabeceras de las columnas (sólo hay una)
                output_vcf.write('##FORMAT=<ID=FMC,Number=.,Type=Integer,Description="Coded tags originally assigned by FilterMutectCalls (except -1, which means variant was not called by Mutect2 but other caller) in FILTER: -1.no_mutect2, 0.PASS, 1.base_qual, 2.clustered_events, 3.contamination, 4.duplicate, 5.fragment, 6.germline, 7.haplotype, 8.low_allele_frac, 9.map_qual, 10.multiallelic, 11.n_ratio, 12.normal_artifact, 13.orientation, 14.panel_of_normals, 15.position, 16.possible_numt, 17.slippage, 18.strand_bias, 19.strict_strand, 20.weak_evidence">\\n')
                output_vcf.write(line + '\\n')

            else: # líneas con las variantes. Campos: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, MUESTRA
                campos = line.split('\\t')

                filter_tag = campos[6] # Tag asignada originalmente por FilterMutectCalls. La movemos a FORMAT.

                campos[8] = campos[8] + ':FMC'
                campos[9] = campos[9] + ':' + ','.join(str(codigo_etiqueta[fmc_tag]) for fmc_tag in filter_tag.split(';'))

                if labels_to_rescue:
                    current_labels = filter_tag.split(';')

                    for label in labels_to_rescue.split(','):
                        if label in current_labels:
                            current_labels.remove(label)

                    campos[6] = ';'.join(current_labels) or 'PASS' # si se ha borrado la única etiqueta, pasa a PASS

                output_vcf.write('\\t'.join(campos) + '\\n')

    # Versions
    command_string = '''cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    '''
    os.system(command_string)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_rescued.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}

process INC_REFORMAT_GENOTYPE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyvcf3=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.3--pyhdfd78af_0' :
        'biocontainers/pyvcf3:1.0.3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_reformatted.vcf"), emit: reformatted_variants
    path "versions.yml"                       , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def max_het_vaf = task.ext.max_het_vaf ?: 0.90 // rutina.sh: -m 0.90
    def min_depth = task.ext.min_depth ?: 20       // rutina.sh: -d 20
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    '''
    inc_reformat_genotype_v1.1.py

    Corrige la información sobre genotipo, frecuencia alélica y allelic depth de un archivo vcf. Recalcula AF a partir
    de AD, asigna GT según AF (0/0 < af <= 0/1 < max_het_vaf <= 1/1) y descarta las variantes con AF < af o con
    profundidad (AD ref + alt) < min_depth. FORMAT queda como GT:AD:AF.

    Versión: 1.1
    Autor: UB INCLIVA
    Fecha: 03/07/2020
    '''
    import gzip
    import os

    min_vaf = ${params.af}
    max_het_vaf = ${max_het_vaf}
    min_depth = ${min_depth}

    def get_genotype_position(format_field):
        gt = ''
        ad = ''

        for field_index, field in enumerate(format_field.split(':')):
            if 'GT' == field:
                gt = field_index

            if 'AD' == field:
                ad = field_index

        return [gt, ad]

    def set_genotype(vaf):
        if vaf < min_vaf:
            return '0/0'

        elif vaf < max_het_vaf:
            return '0/1'

        return '1/1'

    with gzip.open("${vcf}", 'rt') as vcf_reader, open("${prefix}_reformatted.vcf", 'w') as vcf_writer:
        for record in vcf_reader:
            record = record.rstrip()

            if not record.startswith('#'):
                record_fields = record.split('\\t')
                info = get_genotype_position(record_fields[8])

                if record_fields[9].split(':')[info[1]] != '.':
                    ad = record_fields[9].split(':')[info[1]].split(',')
                    vaf = round(float(ad[1]) / (float(ad[0]) + float(ad[1])), 4)
                    gt = set_genotype(vaf)
                    total_depth = int(ad[1]) + int(ad[0])

                    if vaf >= min_vaf and total_depth >= min_depth:
                        record_fields[8] = 'GT:AD:AF'
                        record_fields[9] = gt + ':' + ','.join(ad) + ':' + str(vaf)
                        vcf_writer.write('\\t'.join(record_fields) + '\\n')

            elif record.startswith('##FORMAT=<ID='):
                if record.startswith('##FORMAT=<ID=AF') or record.startswith('##FORMAT=<ID=GT') or record.startswith('##FORMAT=<ID=AD'):
                    vcf_writer.write(record + '\\n')

            else:
                vcf_writer.write(record + '\\n')

    # Versions
    command_string = '''cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    '''
    os.system(command_string)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_reformatted.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}


process INC_TAG_RELIABLE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(vcf)
    path primers

    output:
    tuple val(meta), path("*_tagged.vcf"), emit: tagged_variants
    path "versions.yml"                  , emit: versions

    script:
    // inc_tag_reliable_v1.2.py (UB INCLIVA, 03/09/2020) reescrito con bedtools + awk, ya que no hay un contenedor con
    // python y bedtools. Añade REL en INFO: ext.args (Low_confidence para Mutect2, High_confidence para Lofreq) si la
    // variante cae en zona de primers y "High_Confidence" en caso contrario. La línea de cabecera de REL se escribe
    // justo después de la primera línea ##INFO.
    def prefix = task.ext.prefix ?: "${meta.id}"
    def confidence_tag = task.ext.args ?: ''
    """
    bedtools intersect -header -c -a ${vcf} -b ${primers} \\
        | awk -F '\\t' -v OFS='\\t' -v tag='${confidence_tag}' \\
            -v rel_header='##INFO=<ID=REL,Number=1,Type=String,Description="Variant Confidence">' '
            /^#/ { print; if (!p && /^##INFO/) { print rel_header; p = 1 }; next }
            { rel = (\$NF > 0) ? tag : "High_Confidence"; print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8 ";REL=" rel, \$9, \$10 }
        ' > ${prefix}_tagged.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_tagged.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}

process INC_ADD_FMC_TAG {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyvcf3=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.3--pyhdfd78af_0' :
        'biocontainers/pyvcf3:1.0.3' }"

    input:
    tuple val(meta), path(vcf), path(mutect2_vcf)

    output:
    tuple val(meta), path("*_fmc.vcf"), emit: fmc_variants
    path "versions.yml"               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    '''
    inc_add_FMC_tag_v1.0.py

    Añade al vcf combinado (Mutect2 + Lofreq) el campo FMC de FORMAT con las etiquetas que FilterMutectCalls había
    añadido previamente (se toman del vcf de Mutect2 tras INC_RESCUE_VARIANTS). Las variantes no llamadas por Mutect2
    (y por tanto llamadas por lofreq) aparecerán como -1. Ver INC_RESCUE_VARIANTS para la correspondencia código-etiqueta.

    Versión: 1.0
    Autor: UB INCLIVA - Francisco Martínez Picó
    Fecha: 16/10/2023
    '''
    import vcf
    import gzip
    import os

    # Leemos el archivo de Mutect2 para extraer información de FMC
    vcf_reader = vcf.Reader(open("${mutect2_vcf}", 'rb'))

    for line in vcf_reader._header_lines: # Parseamos la línea del campo FMC de la cabecera ya que también es necesario volver a añadirla
        if line.startswith('##FORMAT=<ID=FMC'):
            fmc_header_line = line

    variants_fmc_tags = {} # id CHROM_POS_REF_ALT -> valor del campo FMC

    for record in vcf_reader:
        v2_id = f'{record.CHROM}_{record.POS}_{record.REF}_{record.ALT[0]}' # vcf unialélicos (separados con bcftools norm)

        sample = record.samples[0] # Sólo hay una muestra (tumoral)

        # FMC puede ser un valor entero (sólo una etiqueta) o una lista de enteros (varias etiquetas). Nos aseguramos que se lee como lista.
        fmc_value = sample.data.FMC
        if isinstance(fmc_value, int):
            fmc_tags = [fmc_value]
        else:
            fmc_tags = fmc_value

        variants_fmc_tags[v2_id] = ','.join(str(fmc_tag) for fmc_tag in fmc_tags)

    open_vcf = gzip.open if "${vcf}".endswith('.gz') else open

    with open_vcf("${vcf}", 'rt') as input_vcf, open("${prefix}_fmc.vcf", 'w') as output_vcf:
        for line in input_vcf:
            line = line.strip()

            if line.startswith('##'): # líneas de la cabecera
                output_vcf.write(line + '\\n')

            elif line.startswith('#'): # línea con las cabeceras de las columnas (sólo hay una)
                output_vcf.write(fmc_header_line + '\\n')
                output_vcf.write(line + '\\n')

            else: # líneas de las variantes
                splitted_line = line.split()

                v_id = f'{splitted_line[0]}_{splitted_line[1]}_{splitted_line[3]}_{splitted_line[4]}' # id es CHROM_POS_REF_ALT

                splitted_line[8] = splitted_line[8] + ':FMC'
                splitted_line[9] = splitted_line[9] + ':' + variants_fmc_tags.get(v_id, '-1') # -1 si no ha sido llamada por Mutect2

                output_vcf.write('\\t'.join(splitted_line) + '\\n')

    # Versions
    command_string = '''cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    '''
    os.system(command_string)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fmc.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}

process INC_ADD_CSB_INFO {
    tag "$meta.id"
    label 'process_single'

    // Imagen de GATK 4.6 (la misma que GATK4_MUTECT2): es la que tenemos con python3 + PyVCF + numpy + scipy.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
    tuple val(meta), path(vcf), path(mutect2_vcf), path(lofreq_vcf)

    output:
    tuple val(meta), path("*_csb.vcf"), emit: csb_variants
    path "versions.yml"               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    '''
    inc_add_cSB_info_v1.0.py

    Añade al vcf combinado el campo cSB de FORMAT con el strand bias (phred score del p-valor del test exacto de Fisher)
    calculado para cada caller: primero Mutect2 (a partir de SB) y después Lofreq (a partir de DP4). Si la variante no ha
    sido llamada por alguno de los callers, aparece un '.'.

    Phred Score (Q) = -10 * log10(P). Un Phred Score por encima de 20 significa p-valor < 0.01 (13 para p-valor < 0.05 y
    28 para p-valor < 0.001).

    Versión: 1.0
    Autor: UB INCLIVA - Francisco Martínez Picó
    Fecha: 25/10/2023
    '''
    import gzip
    import math
    import os
    import numpy as np
    import vcf
    from scipy.stats import fisher_exact

    def create_mutect2_variant_ids(vcf_reader):
        # id CHROM_POS_REF_ALT -> SB (REF forward, REF reverse, ALT forward, ALT reverse)
        variants_id = {}

        for record in vcf_reader:
            v1_id = f'{record.CHROM}_{record.POS}_{record.REF}_{record.ALT[0]}' # vcf unialélicos
            variants_id[v1_id] = record.samples[0].data.SB

        return variants_id

    def create_lofreq_variant_ids(vcf_reader):
        # id CHROM_POS_REF_ALT -> DP4 (REF forward, REF reverse, ALT forward, ALT reverse)
        variants_id = {}

        for record in vcf_reader:
            v1_id = f'{record.CHROM}_{record.POS}_{record.REF}_{record.ALT[0]}' # vcf unialélicos
            variants_id[v1_id] = record.INFO['DP4']

        return variants_id

    def calculate_phred_score(variants_caller):
        caller_phred_scores = {}

        for key in variants_caller.keys():
            sb_values = variants_caller[key]

            contingency_table = np.array([[sb_values[0], sb_values[1]], [sb_values[2], sb_values[3]]])
            odds_ratio, p_value = fisher_exact(contingency_table)

            if p_value == 0.0: # log10(0) no existe: asignamos un valor muy cercano a 0
                p_value = 0.0000000001

            caller_phred_scores[key] = str(round(math.log10(p_value) * -10))

        return caller_phred_scores

    def get_csb_value(variant_id, callers_sb_info):
        csb_value = []

        for caller in callers_sb_info.keys():
            csb_value.append('.')

            if variant_id in callers_sb_info[caller][0]:
                csb_value[callers_sb_info[caller][1]] = callers_sb_info[caller][0][variant_id]

        return csb_value

    mutect2_phred = calculate_phred_score(create_mutect2_variant_ids(vcf.Reader(open("${mutect2_vcf}", 'rb'))))
    lofreq_phred = calculate_phred_score(create_lofreq_variant_ids(vcf.Reader(open("${lofreq_vcf}", 'rb'))))

    # Posición del valor SB de cada caller en el campo cSB. Actualizar si se añaden más callers.
    callers_sb_info = { 'mutect2': [mutect2_phred, 0],
                        'lofreq': [lofreq_phred, 1]}

    sb_header_line = '##FORMAT=<ID=cSB,Number=2,Type=Integer,Description="Phred-scaled strand bias at this position. Calculated using Fisher\\'s exact test. This contains info from all the callers used in the pipeline. In order: 0.mutect2, 1.lofreq">'

    open_vcf = gzip.open if "${vcf}".endswith('.gz') else open

    with open_vcf("${vcf}", 'rt') as input_vcf, open("${prefix}_csb.vcf", 'w') as output_vcf:
        for line in input_vcf:
            line = line.strip()

            if line.startswith('##'): # líneas de la cabecera
                output_vcf.write(line + '\\n')

            elif line.startswith('#'): # línea con las cabeceras de las columnas (sólo hay una)
                output_vcf.write(sb_header_line + '\\n')
                output_vcf.write(line + '\\n')

            else: # líneas de las variantes
                splitted_line = line.split()
                splitted_line[8] = splitted_line[8] + ':cSB'

                variant_id = f'{splitted_line[0]}_{splitted_line[1]}_{splitted_line[3]}_{splitted_line[4]}' # id es CHROM_POS_REF_ALT

                splitted_line[9] = splitted_line[9] + ':' + ','.join(get_csb_value(variant_id, callers_sb_info))

                output_vcf.write('\\t'.join(splitted_line) + '\\n')

    # Versions
    command_string = '''cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed -e "s/Python //g")
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    '''
    os.system(command_string)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_csb.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}

process INC_BCFTOOLS_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data' :
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    path "versions.yml"                      , emit: versions

    script:
    // A diferencia del módulo nf-core BCFTOOLS_MERGE, no reordena los vcfs por nombre: se unen en el orden recibido
    // (el del samplesheet, como la lista "samples" de rutina.sh). El orden decide las columnas de muestra y de qué
    // muestra se copian QUAL e INFO en las variantes compartidas.
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools merge \\
        ${args} \\
        --output ${prefix}.vcf.gz \\
        ${vcfs.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

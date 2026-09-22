process GET_COVERED_BED {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(depth_tsv)

    output:
    tuple val(meta), path("*_covered_regions_${params.cov}x.bed"), emit: covered_regions
    path "versions.yml"                                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '{if (\$3 >= '${params.cov}') print \$1"\t"\$2-1"\t"\$2}' ${depth_tsv} | bedtools merge ${args} > ${prefix}_covered_regions_${params.cov}x.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_covered_regions_${params.cov}x.bed

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

    for record in vcf_reader:
        for sample_index, sample in enumerate(record.samples):

            if sample.called:
                sample_vaf = round(float(sample.data.AD[1]) / (float(sample.data.AD[0]) + float(sample.data.AD[1])), 4)
                total_depth = sample.data.AD[1] + sample.data.AD[0]

                if sample_vaf >= ${params.af} and total_depth >= ${params.cov}:
                    vcf_writer.write_record(record)

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

    for record in vcf_reader:
        sample_vaf = record.INFO['AF']
        total_depth = record.INFO['DP']

        if sample_vaf >= ${params.af} and total_depth >= ${params.cov}:
            vcf_writer.write_record(record)

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

// process INC_TAG_RELIABLE { // TODO - acabar este módulo
//     tag "$meta.id"
//     label 'process_single'

//     conda "bioconda::pyvcf3=1.0.3"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/pyvcf3:1.0.3--pyhdfd78af_0' :
//         'biocontainers/pyvcf3:1.0.3' }"

//     input:
//     tuple val(meta), path(vcf)

//     output:
//     tuple val(meta), path("*_conf.vcf"), emit: conf_variants
//     path "versions.yml"                , emit: versions

//     script:
//     def prefix = task.ext.prefix ?: "${meta.id}"
//     """
//     #!/usr/bin/env python3
//     # -*- coding: utf-8 -*-
//     '''
//     inc_gatk_reliable.py

//     Añade información sobre la confianza de una variable. La tag sobre la confianza se pasa como argumento. En general,
//     al utilizarlo se considera gatk Low_confidence y lofreq High_confidence. Antiguamente llamado
//     "inc_gatk_reliable.py", pero renombrado porque ya no solo se aplica a vcfs generados por gatk.

//     Versión: 1.0
//     Autor: UB INCLIVA
//     Fecha: 03/09/2020
//     '''
//     output = subprocess.check_output("bedtools intersect -a " + args.file +
//                                     " -b "+ args.primers + " | awk -F " + '"\t" ' +
//                                     "'{print $1" + '"_"$2"_"$4"_"$5}' + "'", shell= True)

//     output = output.decode() # bytes object since python3, convert to string
//     vars_ = output.split("\n")

//     file = gzip.open(args.file, 'rt')
//     f = file.readlines()
//     file.close()

//     outfile = open(args.outfile, 'w')

//     p = 0

//     for line in f:
//         line = line.rstrip()

//         if line.startswith("#"):

//             if line.startswith('##INFO') and p == 0:
//                 outfile.write(line + "\n")
//                 outfile.write('##INFO=<ID=REL,Number=1,Type=String,Description="Variant Confidence">' + "\n")
//                 p = 1

//             else:
//                 outfile.write(line + "\n")

//             continue

//         values = line.split("\t")
//         line = "\t".join(values[0:7]) + "\t" + values[7] + ";REL=High_Confidence" + "\t" + values[8] + "\t" + values[9]

//         id_ = values[0] + "_" + values[1] + "_" + values[3] + "_" + values[4]

//         for var in vars_:
//             if var != "":
//                 if id_ == var:
//                     line = "\t".join(values[0:7]) + "\t" + values[7] + ";REL=" + args.tag +  "\t" + values[8] + "\t" +  values[9]
//         outfile.write(line + "\n")

//     outfile.close()

//     # Versions
//     command_string = '''cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         python: \$(python --version | sed -e "s/Python //g")
//     END_VERSIONS
//     '''
//     os.system(command_string)
//     """

//     stub:
//     def prefix = task.ext.prefix ?: "${meta.id}"
//     """
//     touch ${prefix}_inc_formated.vcf

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         python: \$(python --version | sed -e "s/Python //g")
//     END_VERSIONS
//     """
// }

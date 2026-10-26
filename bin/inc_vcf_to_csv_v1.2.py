#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
inc_vcf_to_csv_v1.2.py

Este rograma convierte el vcf.gz anotado en una archivo csv con la información distribuida por columnas. En el csv final
sólo aparecerán aquellas variantes del vcf que (1) hayan sido anotadas (tengan campo ANN en la columna de INFO),
(2) sean de biotipo "protein_coding" y que (3) tengan un transcrito asociado en la anotación que aparezca en el
"list_of_transcripts" proporcionado.

Versión: 1.2
Autor: UB INCLIVA - Francisco Martínez Picó, Roxana Andreea Moldovan
Fecha: 11/10/2023
"""
import vcf
import re
import argparse
import sys

def get_sample_names(sample_dict):
    """
    This is just in case the list of sample names in the vcf_header accessed by vcf_hearder.samples does not mantain the
    right sample order. Obtiene el nombre de las muestras.

    INPUT
    - sample_dict (dict): diccionario cuyas claves son el nombre de las muestras y los valores su posición en las
    columnas del vcf.

    RETURN
    - samples (list): lista cuyo cada elemnto es el nombre de una muestra.
    """
    samples = []

    for sample_index in range(0, len(sample_dict)):
        for keys, values in sample_dict.items():
            if sample_index == values:
                samples.append(keys)

    return samples

def parse_record_info(record_info):
    """
    Parsea la información de cada "record" (variante) del vcf.

    INPUT
    - record_info (dict): del paquete PyVCF, tiene la información de cada variante.

    RETURN
    - new_info (dict): contiene la anotación del campo "ANN" del vcf.
    """
    csq = []
    new_info = {}

    is_annotated = False

    for keys, values in record_info.items():
        if keys == 'ANN':
            is_annotated = True

    if is_annotated:
        for annotation in record_info['ANN']:
            csq.append(annotation)

        new_info['ANN'] = ','.join(csq)

    return new_info

def parse_sample_records(record_samples, min_heterozygous, max_heterozygous):
    """
    Sirve para obtener la información del campo FORMAT de las muestras.

    INPUTS
    - record_samples (list): para cada variante, su distinto valor a través de las diferentes muestras.
    - min_heterozygous (float): valor mínimo de af para considerar una variante heterozigota.
    - max_hterozygous (float): valor máximo de af para considerar una variante heterozigota.

    RETURN
    - sample_format (lista):lista que contien cada uno de los campos FORMAT para cada una de las variantes por muestra.
    """
    sample_format = []

    for sample in record_samples:
        try: # la muestra puede tener la variante o no (sample.called True|False)
            allele_freq = round(float(sample.data.AD[1]) / (float(sample.data.AD[0]) + float(sample.data.AD[1])), 4)

            depth = sample.data.AD[1] + sample.data.AD[0]

        except: # Caso en el que una variante no está en una muestra
            genotype = '.'
            sample_format.append(genotype)

            continue

        # Se parsea FMC, aquí siempre habrá un número porque, si no lo ha llamado mutect, tendrá '-1' (y si no lo llama ningún caller pues obviamente no estaría la variante):
        sample_fmc = sample.data.FMC  # puede ser int o lista

        # si no es una lista, lo convertimos en lista de un solo elemento
        if not isinstance(sample_fmc, (list, tuple)):
            sample_fmc = [sample_fmc]

        # convertimos a str y unimos
        sample_fmc = [str(numero) for numero in sample_fmc]
        sample_fmc = '_'.join(sample_fmc)

        # Se parsea cSB (en el caso de las variantes llamadas sólo por un caller aparece un punto '.' en el otro caller, que se parsea como None y da problemas):
        for index, value in enumerate(sample.data.cSB):
            if value is None: # Sustituimos valores 'None' por un caracter
                sample.data.cSB[index] = '-'

        sample_csb = sample.data.cSB # lista que puede tener o no varios elementos
        sample_csb = [str(numero) for numero in sample_csb] # Se convierte lista de in a lista de str
        sample_csb = '_'.join(sample_csb)

        # FORMAT se devolverá como GT:AF:DP:FMC:cSB

        # Se induce el genotipo y se crea el campo FORMAT:
        if allele_freq > min_heterozygous and allele_freq <= max_heterozygous:
            genotype = '0/1'
            sample_format.append(genotype + ':' + str(allele_freq) + ':' + str(depth) + ':' + sample_fmc+ ':' + sample_csb)

        elif allele_freq <= min_heterozygous:
            sample_format.append('.')

        else:
            genotype = '1/1'
            sample_format.append(genotype + ':' + str(allele_freq) + ':' + str(depth) + ':' + sample_fmc + ':' + sample_csb)

    return sample_format

def get_ann_information(vcf_reader):
    """
    Coge la línea ANN del header y obtiene el orden en el que se pone la información.

    INPUT
    - vcf_reader (vcf.parser.Reader): objeto del paquete vcf. Reader del archivo vcf.

    RETURN
    - anntotation (dict)
    """
    annotation = {}

    for line in vcf_reader._header_lines:
        if "INFO=<ID=ANN" in line:
            line = re.sub(r'#.*Format: ', '', line)
            line = re.sub(r'\">', '', line)

            # Se crea un diccionario que tiene { campo_ANN: su índice, ... }
            for index, field in enumerate(line.split('|')):
                annotation[str(field)] = index

    return annotation

def write_csv(output_file_handler, record_chrom, record_pos, record_ref, allele, allele_info, samples_format, ann_info,
              refseq_ids, pathogenic_coordinates, record_id, mixed, origin, sig, conser, rel, variant_in_oncokb, treatment,
              oncogenic, hsl, pangolin, fmc_tag, called_by_mutect2, csb_tag):
    """
    Escribe el archivo csv. Cada línea corresponde con una variante. Sigue formato csv.

    INPUTS
    - output_file_handler ():
    - record_chrom ():
    - record_pos ():
    - record_ref ():
    - allele ():
    - allele_info ():
    - samples_format ():
    - ann_info ():
    - refseq_ids ()
    - pathogenic_coordinates ():
    - record_id ():
    - mixed ():
    - origin ():
    - sig ():
    - conser():
    - rel ():
    - variant_in_oncokb ():
    - treatment ():
    - oncogenic ():
    - hsl ():
    - pangolin ():
    - fmc_tag ():
    - called_by_mutect2 ():
    - csb_tag ():
    """
    record_allele_ann = {}

    is_annotated = False

    for keys, values in allele_info.items():
        if keys == 'ANN':
            is_annotated = True

    if is_annotated:
        formatted_samples = []

        for sample in samples_format:
            formatted_samples.append(parse_sample_info_to_csv(sample))

        for annotation in allele_info['ANN'].split(','): # Todo el campo de ANN
            for annotation_field, annotation_index in ann_info.items():

                # Se va rellenando el diccionario record_allele_ann con la información del campo ANN (si algo no está se pone "-"):
                for index, field in enumerate(annotation.split('|')):
                    if annotation_index == index:
                        if field == '':
                            record_allele_ann[annotation_field] = '-'
                        else:
                            record_allele_ann[annotation_field] = field

            there_is_a_variant = 0

            for one_sample in formatted_samples: # Muestras formateadas como "AF,DP"
                if one_sample.split(",")[0] != "-": # Si no hay variante tiene valor "-"
                    there_is_a_variant += 1

            if there_is_a_variant > 0 and record_allele_ann['BIOTYPE'] == "protein_coding":
                # Inicialmente asumimos no-patogénica y frecuencia 0, iremos iterando para comrprobarlo
                is_pathogenic = "-"
                allele_freq = 0
                add_frequencies = []

                if record_id is not None:
                    is_pathogenic = "yes"

                for sample in formatted_samples:
                    if sample.split(',')[0] != "-":
                        allele_freq += 1
                        add_frequencies.append(float(sample.split(',')[0]))

                max_frequency = max(add_frequencies)

                record_allele_ann['GAF'] = allele_freq # En cuantas muestras está la variante

                if refseq_ids is not None:
                    if record_allele_ann['SYMBOL'] in refseq_ids.keys() and refseq_ids[record_allele_ann['SYMBOL']] in record_allele_ann['RefSeq']:
                        # Sometimes the refseq transcript is within a list of transcripts separated with &
                        if record_allele_ann['HGVSc'] != "-":
                            record_allele_ann['HGVSc'] = refseq_ids[record_allele_ann['SYMBOL']] + ":" + \
                                                            record_allele_ann['HGVSc'].split(':')[1]
                        else:
                            record_allele_ann['HGVSc'] = refseq_ids[record_allele_ann['SYMBOL']]

                        output_file_handler.write(
                            format_output(record_allele_ann, record_chrom, record_pos, record_ref, allele.sequence,
                                            formatted_samples, pathogenic_coordinates, is_pathogenic,
                                            max_frequency, mixed, origin, sig, conser, rel, variant_in_oncokb,
                                            treatment, oncogenic, hsl, pangolin, fmc_tag, called_by_mutect2, csb_tag))

                else:
                    output_file_handler.write(
                        format_output(record_allele_ann, record_chrom, record_pos, record_ref, allele.sequence,
                                      formatted_samples, pathogenic_coordinates, is_pathogenic, max_frequency,
                                      mixed, origin, sig, conser, rel, variant_in_oncokb, treatment, oncogenic, hsl,
                                      pangolin, fmc_tag, called_by_mutect2, csb_tag))

def format_output(record_allele_ann, record_chrom, record_pos, record_ref, allele_sequence, formatted_samples,
                  pathogenic_coordinates, is_pathogenic, max_frequency, mixed, origin, sig, conser, rel,
                  variant_in_oncokb, treatment, oncogenic, hsl, pangolin, fmc_tag, called_by_mutect2, csb_tag):
    """
    Devuelve la string a escibir en el csv final. Utilizado en la función write_csv.

    INPUTS

    RETURN
    - output_string (str): línea que se escribirá en el csv final. Cada variante una línea. Está formateada para
    coincidir con las cabeceras del archivo csv.
    """
    if ":" in record_allele_ann['HGVSp']:
        record_allele_ann['HGVSp'] = record_allele_ann['HGVSp'].split(':')[1]
        if "%3D" in record_allele_ann['HGVSp']:
            record_allele_ann['HGVSp'] = record_allele_ann['HGVSp'].replace("%3D", "=")
    if pathogenic_coordinates == "HS":
        pathogenic_coordinates = "yes"

    # SpliceAI
    if 'SpliceAI_cutoff' in record_allele_ann and record_allele_ann['SpliceAI_cutoff'] == 'PASS':
        spliceai_pred = (record_allele_ann['SpliceAI_pred_DP_AG'] + ':' +
                         record_allele_ann['SpliceAI_pred_DS_AG'] + '|' +
                         record_allele_ann['SpliceAI_pred_DP_AL'] + ':' +
                         record_allele_ann['SpliceAI_pred_DS_AL'] + '|' +
                         record_allele_ann['SpliceAI_pred_DP_DG'] + ':' +
                         record_allele_ann['SpliceAI_pred_DS_DG'] + '|' +
                         record_allele_ann['SpliceAI_pred_DP_DL'] + ':' +
                         record_allele_ann['SpliceAI_pred_DS_DL'])
    else:
        spliceai_pred = '-'

    # Called by mutect2 or not
    if called_by_mutect2:
        mutect2 = 'Yes'

    else:
        mutect2 = 'No'

    output_string = (record_allele_ann['SYMBOL'] + ',' +                 # GENE_SYMBOL
                     record_chrom + ':' + str(+record_pos) + ',' +       # POSITION
                     record_ref + ',' +                                  # REF
                     allele_sequence + ',' +                             # ALT
                     str(record_allele_ann['GAF']) + ',' +               # GAF
                     str(max_frequency) + ',' +                          # MAX_FREQ
                     record_allele_ann['Consequence'] + ',' +            # CONSEQUENCE
                     record_allele_ann['IMPACT'] + ',' +                 # IMPACT
                     record_allele_ann['EXON'].replace("/", "|") + ',' + # EXON
                     record_allele_ann['HGVSc'] + ',' +                  # HGVSc
                     record_allele_ann['HGVSp'] + ',' +                  # HGVSp
                     is_pathogenic + ',' +                               # INCLIVA_PATH_MUTATION
                     pathogenic_coordinates + ',' +                      # HOTSPOT
                     mixed + "," +                                       # MIXED_DB
                     origin + "," +                                      # MIXED_ORIGIN
                     sig + "," +                                         # MIXED_CLIN_SIG
                     record_allele_ann['Existing_variation'] + ',' +     # VARIANT_ID
                     record_allele_ann['AF'] + ',' +                     # 1000G_AF
                     record_allele_ann['MAX_AF'] + ',' +                 # MAX_AF
                     record_allele_ann['gnomADe_AF'] + ',' +             # GNOMAD
                     record_allele_ann['CLIN_SIG'] + ',' +               # CLINVAR
                     record_allele_ann['SIFT'] + ',' +                   # SIFT
                     record_allele_ann['PolyPhen'] + ',' +               # POLYPHEN
                     conser + ',' +                                      # CONSERVATION
                     spliceai_pred + ',' +                               # SPLICEAI_PRED
                     pangolin + ',' +                                    # PANGOLIN_PRED
                     record_allele_ann['PUBMED'] + ',' +                 # PUBMED
                     rel +                                               # CONFIDENCE - (comma included in variable value assignment, since this field is only used in amplicons)
                     variant_in_oncokb + ',' +                           # REVIEWED_BY_ONCOKB
                     oncogenic + ',' +                                   # ONCOGENIC
                     ''.join(treatment) + ',' +                          # TREATMENT
                     hsl + "," +                                         # EVIDENCE_LEVEL
                     fmc_tag + ',' +                                     # FILTERMUTECTCALLS
                     mutect2 + ',' +                                     # MUTECT2
                     csb_tag + ',' +                                     # STRAND_BIAS_BY_CALLER
                     ','.join(formatted_samples) +                       # MUESTRAS
                     '\n')

    return output_string

def parse_sample_info_to_csv(sample):
    """
    Extrae la información de frecuencia de la variante (AF) y profundidad (DP). En caso de no estar la variante en la muestra, se devuelve "-" separado por ",".

    INPUT
    - sample (str): cadena del tipo "0/1:0.0306:6007" que contiene información sobre GT:AF:AD.

    RETURN
    - formatted_sample (str): cadena formateada que tiene la información de AF y DP para una variante en una muestra.
    """
    if sample == '.': # Muestras en las que la variante no está
        formatted_sample = '-,-'

    else: # Se extrae la información
        sample_information = sample.split(':')
        allele_freq = float(sample_information[1])
        formatted_sample = (str(allele_freq) + ',' + sample_information[2])

    return formatted_sample

def main(arguments):
    vcf_reader = vcf.Reader(open(arguments.input_vcf, 'rb'))
    output_file_handler = open(arguments.output_file, 'w')

    samples_in_header = []
    refseq_ids = {}

    if arguments.select_refseq:
        # list_of_transcripts es un archivo con 3 columnas, cada línea es un tránscrito.
        for transcript in open(arguments.select_refseq, 'r'):
            # Se crea cada tránscrito como una clave del diccionario refseq_ids, cuyos valores serán los id de los tránscritos.
            refseq_ids[transcript.rstrip().split("\t")[0]] = transcript.rstrip().split("\t")[2]

    else:
        refseq_ids = None

    sample_names = get_sample_names(vcf_reader._sample_indexes)

    # Por cada muestra, tendremos _FREQUENCY (AF) y DP (profundidad de lectura)
    for sample in sample_names:
        samples_in_header.append(sample + '_FREQUENCY,' + sample + '_DP')

    # Estas son las columnas que tendrá el csv. Se añaden las de las muestras.
    csv_header_list = ['GENE_SYMBOL',
                       'POSITION',
                       'REF',
                       'ALT',
                       'GAF',
                       'MAX_FREQ',
                       'CONSEQUENCE',
                       'IMPACT',
                       'EXON',
                       'HGVSc',
                       'HGVSp',
                       'INCLIVA_PATH_MUTATION',
                       'HOTSPOT',
                       'MIXED_DB',
                       'MIXED_ORIGIN',
                       'MIXED_CLIN_SIG',
                       'VARIANT_ID',
                       '1000G_AF',
                       'MAX_AF',
                       'GNOMAD',
                       'CLINVAR',
                       'SIFT',
                       'POLYPHEN',
                       'CONSERVATION',
                       'SPLICEAI_PRED',
                       'PANGOLIN_PRED',
                       'PUBMED',
                       'CONFIDENCE',
                       'REVIEWED_BY_ONCOKB',
                       'ONCOGENIC',
                       'TREATMENT',
                       'EVIDENCE_LEVEL',
                       'FILTERMUTECTCALLS',
                       'MUTECT2',
                       'STRAND_BIAS_BY_CALLER'] + samples_in_header

    # Se escribe la cabecera del csv. Esta depende de se trata de amplicones o no:
    if arguments.amplicon:
        csv_header = ','.join(csv_header_list) + '\n'

    else: # Si no es de amplicones, no se incluye la columna de "CONFIDENCE"
        csv_header_list.remove('CONFIDENCE')
        csv_header = ','.join(csv_header_list) + '\n'

    output_file_handler.write(csv_header)

    for record in vcf_reader: # Se itera sobre las variantes del vcf
        allele_info = parse_record_info(record.INFO)

        samples_format = parse_sample_records(record.samples, arguments.min_frequency, arguments.max_frequency)

        # Se inicializan valores de diferentes campos
        variant_in_oncokb = "False"
        oncogenic = "-"
        treatment = "-"
        hsl = "-"
        pathogenic_coordinates = "-"
        mixed = "-"
        origin = "-"
        sig = "-"
        pangolin = "-"
        rel = '' # Vacío porque si NO estamos en amplicones no se incluye esta columna
        fmc_tag = '-'
        called_by_mutect2 = False # Si es falso escribe '-' y si es True escribe 'Yes'
        csb_tag = '-'

        # ONCOKB
        if 'VARIANT_IN_ONCOKB' in record.INFO.keys(): # record.INFO.keys() es un diccionario que contiene las claves del campo INFO
            variant_in_oncokb = "True"

        # ONCOGENIC
        if 'Oncogenic' in record.INFO.keys():
            oncogenic = record.INFO["Oncogenic"]

        # TREATMENT
        if 'Treatment' in record.INFO.keys():
            treatment = record.INFO["Treatment"]

        # HSL
        if "hSL" in record.INFO.keys():
            hsl = record.INFO["hSL"]

        # HS
        if 'HS' in record.INFO.keys():
            pathogenic_coordinates = "yes"

        # ORIGIN
        if "ORIGIN" in record.INFO.keys():
            mixed = "yes"
            origin = record.INFO["ORIGIN"]
            sig = record.INFO["SIG"]

        # PANGOLIN
        if 'Pangolin' in record.INFO.keys():
            pangolin = record.INFO["Pangolin"]
            pangolin = pangolin[0].split("|")[1:-1]
            if len(pangolin) == 0:
                pangolin = "-"
            else:
                pangolin = "|".join(pangolin)

        # CONSERVATION
        conser = str(record.INFO["Conservation"][0])

        # CONFIDENCE
        if arguments.amplicon: # Si estamos con amplicones se incluye "rel". Este valor aparecerá en la columna "CONFIDENCE".
            rel = '-,'
            if record.INFO["REL"] == "Low_confidence":
                rel = "Low_confidence,"
            else:
                rel = '-' + ','

        # FILTERMUTECTCALLS FILTER TAG
        fmc_included = False

        for line in vcf_reader._header_lines: # Comprobamos las líneas del header para ver si hay campo FMC en FORMAT.
            if '##FORMAT=<ID=FMC' in line:
                fmc_included = True
                break

        if fmc_included:
            fmc_tag = [] # La etiqueta de este campo será FMC asignadas a cada muestra (concatenadas), cada muestra separada por '|' y cada etiqueta separada por '_'

            for sample in samples_format:
                campos = sample.split(':')

                if campos == ['.']: # casos en los que una variante no se encuentre en una muestra
                    campos = ['0/0', '0.0', '0', '.', '.'] # emulamos campo FORMAT

                fmc_tag.append(campos[3])

            fmc_tag = '|'.join(fmc_tag)

        # MUTECT
        for tag in fmc_tag.split('|'): # Muestras separadas por '|'
            if tag != '-1' and tag != '.': # -1 era código para no llamado por mutect2 y . es que no se ha llamado la variante en esa muestra
                called_by_mutect2 = True
                break

        # # STRAND_BIAS_BY_CALLER
        csb_included = False

        for line in vcf_reader._header_lines: # Comprobamos las líneas del header para ver si hay campo FMC en FORMAT.
            if '##FORMAT=<ID=cSB' in line:
                csb_included = True
                break

        if csb_included:
            csb_tag = [] # La etiqueta de este campo será el valor de strand bias asignado por cada caller a cada muestra (concatenados), cada muestra separada por '|' y cada caller separada por '_'

            for sample in samples_format:
                campos = sample.split(':')

                if campos == ['.']: # casos en los que una variante no se encuentre en una muestra
                    campos = ['0/0', '0.0', '0', '-', '-_-'] # emulamos campo FORMAT

                csb_tag.append(campos[4])

            csb_tag = '|'.join(csb_tag)

        ann_info = get_ann_information(vcf_reader)

        # Pasamos a escribir el cuerpo del csv
        write_csv(output_file_handler, record.CHROM, record.POS, record.REF, record.alleles[1], allele_info,
                  samples_format, ann_info, refseq_ids, pathogenic_coordinates, record.ID, mixed, origin, sig, conser,
                  rel, variant_in_oncokb, treatment, oncogenic, hsl, pangolin, fmc_tag, called_by_mutect2, csb_tag)

    output_file_handler.close()

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', action = 'store', type = str, dest = 'input_vcf',
                    default = None, help = 'Full path to the input vcf file')
parser.add_argument('-o', '--output_file', action = 'store', type = str, dest = 'output_file',
                    default = None, help = 'Full path to the output file')
parser.add_argument('-m', '--min_frequency', action='store', type=float, dest='min_frequency',
                    default = None, help = 'Minimum variant frequency to consider a variant as heterozygous')
parser.add_argument('-n', '--max_frequency', action = 'store', type = float, dest = 'max_frequency',
                    default = None, help = 'Maximum variant frequency to consider a variant as heterozygous')
parser.add_argument('-r', '--select_refseq', action = 'store', type = str, dest = 'select_refseq',
                    default = False, help = 'List of RefSeq transcript IDs to filter')
parser.add_argument('-a', '--amplicon', action = 'store_true', dest = 'amplicon',
                    default = False, help = 'True for amplicons, False otherwise. If True, CONFIDENCE column will appear.')

arguments = parser.parse_args()

if not arguments.input_vcf:
    parser.error('ERROR: missing input file!! Please check usage')
    sys.exit()
elif not arguments.output_file:
    parser.error('ERROR: missing output file!! Please check usage')
    sys.exit()

if __name__ == '__main__':
    main(arguments)

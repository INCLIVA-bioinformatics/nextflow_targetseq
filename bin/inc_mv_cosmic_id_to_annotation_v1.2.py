#!/usr/bin/python
import vcf
import argparse
import sys
import re


def get_ann_information(vcf_header):
    annotation = dict()
    for line in vcf_header._header_lines:
        if "INFO=<ID=ANN" in line:
            line = re.sub(r'#.*Format: ', '', line)
            line = re.sub(r'\">', '', line)
            for index, field in enumerate(line.split('|')):
                annotation[str(field)] = index
    return annotation

def parse_vep_annotation(field_position, vep_annotation, cosmic_id):
    if cosmic_id is not None:
        variants = cosmic_id.split(';')
    new_annotation = []
    annotation_ids = vep_annotation[0].split('|')[field_position['Existing_variation']].split('&')
    if annotation_ids == ['']:
        annotation_ids = []
    annotation_ids2 = []

    for an in annotation_ids:
        if an.startswith("COSM") or an.startswith("COSN"):
            pass
        else:
            annotation_ids2.append(an)

    if cosmic_id is not None:
        for id in variants:
            if id not in annotation_ids2:
                    annotation_ids2.append(id)

    all_ids = list(set(annotation_ids2))

    for annotation in vep_annotation:
        one_annotation = []
        for field_index, field_value in enumerate(annotation.split('|')):
            if field_index == field_position['Existing_variation']:
                one_annotation.append('&'.join(all_ids))
            else:
                one_annotation.append(field_value)
        new_annotation.append('|'.join(one_annotation))

    return new_annotation


def main(arguments):
    vcf_reader=vcf.Reader(open(arguments.input_vcf,'r'), compressed=False)
    vcf_writer=vcf.Writer(open(arguments.output_vcf,'w'), template=vcf_reader)
    for record in vcf_reader:
        new_record_ann = parse_vep_annotation(get_ann_information(vcf_reader), record.INFO['ANN'], record.ID)
        record.INFO['ANN'] = new_record_ann
        record.ID = None
        vcf_writer.write_record(record)


parser = argparse.ArgumentParser(description="This script adds the cosmic ID in the VCF ID column to the vep annotation")
parser.add_argument('-i', '--file', action='store', type=str, dest='input_vcf',
                            default=None, help='Full path to the input vcf file')
parser.add_argument('-o', '--output_directory', action='store', type=str, dest='output_vcf',
                            default=None, help='Full path and name of the output file')


arguments = parser.parse_args()
if not arguments.input_vcf or not arguments.output_vcf:
    parser.error('ERROR: missing input file!! Please check usage')
    sys.exit()

main(arguments)
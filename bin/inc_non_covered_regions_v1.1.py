#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
inc_non_covered_regions_v1.1.py

Crea QC/mapping/...non_covered_regions_100.bed.

Versión: 1.1
Autor: UB INCLIVA
Fecha: 16/11/2023
"""
import os
import argparse


def check_dir(dir):
    '''
    This function check if the directory given by the user
    ends with the properly character
    '''
    if dir[-1] == "/":
        dir = dir
    else:
        dir = dir + "/"

    return dir

def dir_exists(dir):
    '''
    This function checks if the directory given by the user
    exists
    '''
    if os.path.exists(dir):
        pass
    else:
        print("ERROR: The entry directory path is wrong, please check the configuration file")

def open_file(g_stats):
    '''
    This function open the global stats file and
    return a list of samples which have more than 75%
    of reads covered at 250x
    '''

    file = open(g_stats, 'rb')

    f = file.readlines()
    file.close()
    f = f[1:]
    sample_list = []
    for line in f:
        line = line.strip()
        values = line.decode('ascii').split(',')
        try:
            if float(values[11]) >= 75:
                sample_list.append(values[0])
        except:
            print(values[11])

    return (sample_list)

def depth_calc(sample, run_path_analysis, run_path_tmp, bed):
    '''
    This function calculate the sample depths in the target region
    '''
    os.system(args.samtools + " depth -d 9898989898 -a -b " + bed + " " + run_path_analysis + "mapping/" +
              sample + ".bam > " + run_path_tmp + "mapping/" + sample + "_alldepth")

def non_cov(run_path_tmp, min_depth):
    '''
    This function get and generate various files
     with non-covered regions
    '''
    for file in os.listdir(run_path_tmp + "mapping/"):
        if file.endswith("alldepth"):
            name = file.split("_")[0]
            os.system("awk '$3<=" + str(min_depth) + " {print $0}' " +
                      run_path_tmp + "mapping/" + file +
                      " > " + run_path_tmp + "mapping/" + name + "_non_cov")

def intersect(run_path_tmp):
    '''
    This function compare each bed file for sample, and select
    those mutations which appear in more than 70% of samples, creating an uniq bed file
    '''
    non_cov_dic = {}
    total_samples = 0
    for file in os.listdir(run_path_tmp + "mapping/"):
        if file.endswith("non_cov"):
            non_file = open(run_path_tmp + "mapping/" + file, 'rb')
            f = non_file.readlines()
            non_file.close()
            total_samples += 1
            for line in f:
                line = line.strip()
                values = line.decode('ascii').split('\t')
                key = values[0] + "_" + values[1]
                if key not in non_cov_dic.keys():
                    non_cov_dic[key] = 1
                else:
                    non_cov_dic[key] += 1

    inter_file = open(run_path_tmp + "mapping/" + "intersect_all.bed", 'wb')
    for key, value in non_cov_dic.items():
        percentaje = float(value) / total_samples
        if percentaje >= 0.7:
            chr = key.split("_")[0]
            ini = key.split("_")[1]
            ini = str(int(ini) - 1)
            fin = key.split("_")[1]
            tab = "\t"
            salto = "\n"
            inter_file.write(chr.encode() + tab.encode() +
                             ini.encode() + tab.encode() + fin.encode() + salto.encode())

    inter_file.close()

def annotate(run_path_tmp, target_bed, min_depth, run_path_analysis):
    '''
    This function sort and merge the bed file created with the coincidence mutations, and annotate it.
    '''
    os.system(args.bedtools + " sort -i " + run_path_tmp + "mapping/" +
              "intersect_all.bed > " + run_path_tmp + "mapping/" +
              "intersect_all_sorted.bed")
    os.system(args.bedtools + " merge -i " + run_path_tmp + "mapping/" +
              "intersect_all_sorted.bed > " +
              run_path_tmp + "mapping/" + "intersect_all_merge.bed")

    tab = '"\\t"'
    quote = "--"

    os.system(args.bedtools + " intersect -wb -a " + run_path_tmp + "mapping/" +
              "intersect_all_merge.bed" +
              " -b " + target_bed + " | awk '{print $1" +
              '"\\t"$2"\\t"$3"\\t"$7}' +
              "' > " + run_path_tmp + "mapping/" + "intersect_all_merge_gene.bed")

    name = args.name if args.name else run_path_tmp.split("/")[3]
    os.system(args.bedtools + " intersect -a " + run_path_tmp + "mapping/" +
              "intersect_all_merge_gene.bed" +
              " -b " + args.inc_db + " -wao | " +
              "awk '{if ($5 ~ /chr/) print $12" +
              tab + "$1" + tab + "$2" + tab + "$3" + tab + "$4; else print " + '"' +
              quote + '\\t"$1"\\t"$2"\\t"' + "$3" + tab + "$4}'" +
              " | sed 's/.*CDS=//' | sed 's/;CNT.*\\tchr/\\tchr/' | awk '{print $2" +
              tab + "$3" + tab + "$4" +
              tab + "$5" + tab + "$1}' > " + run_path_tmp + "mapping/" + name +
              "_common_non_covered_regions_" + str(min_depth) + ".bed")
    os.system("cp " + run_path_tmp + "mapping/" + name +
              "_common_non_covered_regions_" + str(min_depth) + ".bed " + run_path_analysis +
              "QC/mapping/")


def main():
    g_stats = args.g_stats
    run_path_analysis = check_dir(args.run_path_analysis)
    run_path_tmp = check_dir(args.run_path_tmp)
    dir_exists(run_path_analysis)
    min_depth = str(args.min_depth)

    sample_list = open_file(g_stats)

    for sample in sample_list:
        depth_calc(sample, run_path_analysis, run_path_tmp, args.target_bed)

    non_cov(run_path_tmp, min_depth)

    intersect(run_path_tmp)

    annotate(run_path_tmp, args.target_bed, min_depth, run_path_analysis)

# Several arguments are needed to run the python script
parser = argparse.ArgumentParser()
# Arguments are added to the parser
parser.add_argument("-g","--global_stats", dest = "g_stats", help = "Global stats path", type = str)
parser.add_argument("-a", "--analysis_path", dest = "run_path_analysis",help = "Run analysis dir path.", type = str)
parser.add_argument("-t", "--tmp_path", dest = "run_path_tmp", help = "Run tmp dir path.", type = str)
parser.add_argument("-r", "--target", dest = "target_bed", help="Target bed file path", type=str)
parser.add_argument("-c", "--coverage", dest = "min_depth", help = "Minimun depth value to consider a region covered", type = str)
parser.add_argument("-d", "--inc_db", dest = "inc_db", help = "Full path to inc_db.vcf.", type = str)
parser.add_argument("-s", "--samtools", dest = "samtools", help = "Full path a samtools.", type = str)
parser.add_argument("-b", "--bedtools", dest = "bedtools", help = "Full path a bedtools.", type = str)
# Nextflow: nombre de la corrida (en rutina se deduce de la ruta de tmp)
parser.add_argument("-n", "--name", dest = "name", help = "Run name.", type = str, default = None)
args = parser.parse_args()

if __name__ == '__main__':
    main()

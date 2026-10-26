#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
inc_gene_cover_v1.1.py

Crea stats_by_gene.tsv.

Versión: 1.1
Autor: UB INCLIVA
Fecha: 16/11/2023
"""
import glob
import os
import argparse
import math

def check_dir(dir):
    """
    This function check if the directory given by the user
    ends with the properly character
    """
    if dir[-1] == "/":
        dir = dir
    else:
        dir = dir + "/"

    return dir

def dir_exists(dir):
    """
    This function checks if the directory given by the user
    exists
    """
    if os.path.exists(dir):
        pass
    else:
        print("ERROR: The entry directory path is wrong, please check the configuration file.")

def main():
    inc_db = args.inc_db
    samtools = args.samtools
    bedtools = args.bedtools

    run_path_tmp = check_dir(args.run_path_tmp)
    panel_path = check_dir(args.panel_path)

    dir_exists(run_path_tmp)
    dir_exists(panel_path)

    run_path_analysis = check_dir(args.run_path_analysis)
    dir_exists(run_path_analysis)

    for file_stats in os.listdir(f'{run_path_analysis}QC/mapping/'):
        if file_stats.endswith("_mapping_stats.csv"):

            with open(f'{run_path_analysis}QC/mapping/{file_stats}', 'r') as file:
                f = file.readlines()
                f = f[1:] # quitamos cabecera

                sample_list = []

                for line in f: # creamos una lista con el nombre de las muestras
                    line = line.strip()
                    values = line.split(',')
                    sample_list.append(values[0])

                for sample in sample_list: # comprobamos si existe el depth de cada muestra
                    exists = os.path.isfile(f'{run_path_tmp}mapping/{sample}_depth')

                    if exists:
                        continue

                    else:
                        comando1 = f'{samtools} depth -d 9898989898 -a -b {panel_path}target.bed {run_path_analysis}mapping/{sample}.bam > {run_path_tmp}mapping/{sample}_depth'
                        os.system(comando1)

    name = args.name if args.name else run_path_tmp.split("/")[3] # TODO - susceptible de dar problemas

    with open(f'{run_path_tmp}mapping/{name}_stats_by_gene.tsv', 'w') as gene_stats:
        checked = False
        lista_samples = []
        lista_genes = []

        for file_gene in os.listdir(f'{panel_path}by_gene/'):

            os.chdir(f'{run_path_tmp}/mapping/')

            for file_all in glob.glob("*_depth"):

                if file_gene.endswith("_coord.bed") and file_all.endswith("_depth"):
                    sample = file_all.split("_")[0]
                    gene = file_gene.split("_")[0]

                    with open(f'{run_path_tmp}mapping/{file_all}', 'r') as f_all:
                        f = f_all.readlines()

                    with open(f'{run_path_tmp}mapping/{file_all}.bed', 'w') as f_out:
                        for line in f:
                            line = line.strip()
                            values = line.split('\t')
                            f_out.write(f'{values[0]}\t{str(int(values[1]) - 1)}\t{values[1]}\t{values[2]}\n')

                    os.system(bedtools + " intersect -wb -a " + panel_path + "by_gene/" +
                            file_gene + " -b " + run_path_tmp + "mapping/" + file_all + ".bed | awk '" +
                            '{print $1"\\t"$2"\\t"$3"\\t"' + "$7}' > " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + ".bed")

                    with open(f'{run_path_tmp}mapping/{str(sample)}_{str(gene)}.bed', 'r') as f_int:
                        f = f_int.readlines()

                    if len(f) <= 1:
                        os.system(f'rm {run_path_tmp}mapping/{str(sample)}_{str(gene)}.bed')
                        continue

                    os.system("awk '$4 < " + args.cov + " {print $1" + '"\\t"$2"\\t"$3' + "}' " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + ".bed > " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_non.bed")

                    os.system(bedtools + " sort -i " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_non.bed > " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_sorted.bed")
                    os.system(bedtools + " merge -i " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_sorted.bed > " +
                            run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_merge.bed")

                    os.system(bedtools + " intersect -wb -a " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_merge.bed" +
                            " -b " + panel_path + "target.bed | awk '{print $1" +
                            '"\\t"$2"\\t"$3"\\t"$7}' +
                            "' > " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_inter_merge.bed")

                    tab = '"\\t"'
                    quote = "--"
                    os.system(bedtools + " intersect -a " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_inter_merge.bed" +
                            " -b " + inc_db +" -wao | " +
                            "awk '{if ($5 ~ /chr/) print $12" +
                            tab + "$1" + tab + "$2" + tab + "$3" + tab + "$4; else print " + '"' +
                            quote + '\\t"$1"\\t"$2"\\t"' + "$3" + tab + "$4}'" +
                            " | sed 's/;CNT.*\\tchr/\\tchr/' | awk '{print $2" +
                            tab + "$3" + tab + "$4" +
                            tab + "$5" + tab + "$1}' > " + run_path_tmp + "mapping/" +
                            str(sample) + "_" + str(gene) + "_anno.bed")

                    with open(f'{run_path_tmp}mapping/{str(sample)}_{str(gene)}_anno.bed', 'r') as file_a:
                        f_a = file_a.readlines()

                    with open(f'{run_path_tmp}mapping/{str(sample)}_{str(gene)}_anno2.bed', 'w') as outfile_a:
                        for line in f_a:
                            line = line.rstrip()
                            values = line.split("\t")

                            if values[4] != "--":
                                cds = values[4].split(";")[1]
                                aa = values[4].split(";")[0]
                                c = cds.split("=")[1]
                                a = aa.split("=")[1]
                                line = f'{values[0]}\t{values[1]}\t{values[2]}\t{values[3]}\t{c} ({a})'
                            outfile_a.write(line + "\n")

                    with open(f'{run_path_tmp}mapping/{str(sample)}_{str(gene)}_anno2.bed', 'r') as vars_f:
                        vars_lines = vars_f.readlines()

                    total_vars = 0
                    vars_s = ""

                    for line in vars_lines:
                        values = line.split("\t")

                        if len(values[4]) <= 3:
                            leng = 0
                            v = ""

                        else:
                            leng = 1
                            v = values[4].strip("\n") + "|"
                        total_vars = total_vars + leng
                        vars_s = vars_s + v

                    total = 0
                    cov = 0

                    for line in f:
                        line = line.strip()
                        values = line.split('\t')

                        if int(values[3]) >= int(args.cov):
                            cov = cov + 1
                        total = total + 1

                    if cov == 0:
                        perc = str(0)

                    else:
                        perc = float(cov) / total
                        perc = str(math.floor(perc * 10000) / 100)

                    if not checked:
                        gene_stats.write("#GENE\t")
                        checked = True

                    line = f'%COVERED ({str(sample)})\t# NON COVERED VARIANTS ({str(sample)})\tNON COVERED VARIANTS ({str(sample)})'

                    if line not in lista_samples:
                        lista_samples.append(line)

                    line = f'{gene}\t{perc}\t{str(total_vars)}\t{vars_s[:-1]}'

                    lista_genes.append(line)

                    os.system(f'rm -r {run_path_tmp}mapping/{str(sample)}_{str(gene)}_inter_merge.bed')
                    os.system(f'rm -r {run_path_tmp}mapping/{str(sample)}_{str(gene)}_merge.bed')
                    os.system(f'rm -r {run_path_tmp}mapping/{str(sample)}_{str(gene)}_sorted.bed')
                    os.system(f'rm -r {run_path_tmp}mapping/{str(sample)}_{str(gene)}_non.bed')

        gene_stats.write("\t".join(lista_samples))
        lista_gen = []

        for stat in lista_genes:
            gen = stat.split("\t")[0]
            if gen not in lista_gen:
                gene_stats.write("\n")
                gene_stats.write(stat)
                lista_gen.append(gen)
            else:
                resto = stat.split("\t")[1:]
                gene_stats.write("\t" + "\t".join(resto))

    os.system(f'cp {run_path_tmp}mapping/{name}_stats_by_gene.tsv {run_path_analysis}QC/mapping/')

# Several arguments are needed to run the python script
parser = argparse.ArgumentParser()
# Arguments are added to the parser
parser.add_argument("-a", "--analysis_path", dest = "run_path_analysis",help = "Run analysis dir path.", type = str)
parser.add_argument("-t", "--tmp_path", dest = "run_path_tmp", help = "Run tmp dir path.", type = str)
parser.add_argument("-p", "--panel_path", dest = "panel_path", help = "Panel path.", type = str)
parser.add_argument("-c", "--coverage", dest = "cov", help = "Minimum coverage to consider a region covered.", type = str)
parser.add_argument("-d", "--inc_db", dest = "inc_db", help = "Full path to inc_db.vcf.", type = str)
parser.add_argument("-s", "--samtools", dest = "samtools", help = "Full path a samtools.", type = str)
parser.add_argument("-b", "--bedtools", dest = "bedtools", help = "Full path a bedtools.", type = str)
# Nextflow: nombre de la corrida (en rutina se deduce de la ruta de tmp)
parser.add_argument("-n", "--name", dest = "name", help = "Run name.", type = str, default = None)

args = parser.parse_args()

if __name__ == '__main__':
    main()

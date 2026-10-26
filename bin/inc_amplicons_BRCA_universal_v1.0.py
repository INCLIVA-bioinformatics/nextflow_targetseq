import os
import re
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
        print ("ERROR: The entry directory path is wrong," \
            " please check the configuration file")


# Several arguments are needed to run the python script
parser = argparse.ArgumentParser()
# Arguments are added to the parser
parser.add_argument("dir", help="Input directory", type=str)
parser.add_argument("name", help="First sample name", type=str)
parser.add_argument("cov", help="Minimum coverage value to consider a region covered", type=str)
parser.add_argument("panel", help="Panel path", type=str)
args = parser.parse_args()

dir = args.dir
dir = check_dir(dir)
dir_exists(dir)
panel = check_dir(args.panel)
os.system("mkdir -p " + dir + "amplicons")
os.system("mkdir -p " + dir + "amplicons/def")

i = 0
for bam in os.listdir(dir + "/mapping/"):
    if bam.endswith(".bam"):
        name = bam.split(".")[0]
        os.system("samtools view -h " + dir + "/mapping/" + bam + " > " + dir + "amplicons/" + name + ".sam")
    else:
        continue

    file = open(dir + "amplicons/" + name + ".sam", 'rb')
    f = file.readlines()
    file.close()

    p_file = open(panel + "capture_amplicons.bed", 'rb')
    p_f = p_file.readlines()
    p_file.close()



    for line in p_f:
        line = line.rstrip()
        values = line.split("\t")
        amp = int(values[2]) - int(values[1]) + 1

        outfile = open(dir + "amplicons/" + name + "_" + values[3] + ".sam", 'wb')

        for l in f:
            l = l.rstrip()
            if l.startswith("@"):
                outfile.write(l + "\n")
            else:
                values2 = l.split("\t")

                if values2[2] == values[0]:
                    start = int(values[1]) - 40
                    end = int(values[2]) + 40
                    if int(values2[8]) > 0:

                        start2 = int(values2[3])
                        stri = values2[5]

                        long = [int(s) for s in re.findall(r'-?\d+\.?\d*', stri)]
                        if len(long) == 0:
                            continue
                        lon = max(long)
                        if lon < (amp*0.60 ):
                            continue

                        end2 = int(values2[3]) + lon

                    else:
                        start2 = int(values2[3])
                        stri = values2[5]

                        long = [int(s) for s in re.findall(r'-?\d+\.?\d*', stri)]
                        if len(long) == 0:
                            continue
                        lon = max(long)
                        if lon < (amp*0.60 ):
                            continue

                        end2 = int(values2[3]) + lon

                    if start2 >= start and end2 <= end:
                        outfile.write(l + "\n")

        outfile.close()


        os.system("samtools view -Sb " + dir + "amplicons/" + name + "_" + values[3] + ".sam > " + dir + "amplicons/" + name + "_" + values[3] + ".bam")
        os.system("samtools index " + dir + "amplicons/" + name + "_" + values[3] + ".bam")
        os.system("samtools depth -d 9898989898 -a " + dir + "amplicons/" + name + "_" + values[3] + ".bam -r " +
                  values[0] + ":" + str(int(values[1]) + 1) + "-" + str(int(values[2])) +
                  " | awk '{print $1" + '"\\t"$2-1"\\t"$2"\\t"$3}' + "' > " + dir + "amplicons/" + name + "_" + values[3] + "_amplicons_depth")


        amplicons = {}
        amplicons[values[0] + "_" + values[1] + "_" + values[2] + "_" + values[3]] = [0, 9898989898, 0, 0]

        file3 = open(dir + "amplicons/" + name + "_" + values[3] + "_amplicons_depth", 'rb')
        f3 = file3.readlines()
        file3.close()


        if len(f3) == 0:
            outfile = open(dir + "amplicons/" + name + "_" + values[3] + "_amplicons_stats", 'wb')
            outfile.write(values[0] + "\t" + values[1] + "\t" + values[2] + "\t" + values[3] + "\t" +
                          "0" + "\t" + "0.0" + "\t" + "0" + "\t" + "0" + "\n")
            outfile.close()
            continue

        for line3 in f3:

            line3 = line3.rstrip()
            values3 = line3.split("\t")
            for amplicon in amplicons.keys():
                stats = amplicon.split("_")
                if values3[0] == stats[0]:
                    if int(values3[3]) < amplicons[amplicon][1] and values3[2] <= stats[2] and values3[2] >= stats[1]:
                        amplicons[amplicon][1] = int(values3[3])

                    if int(values3[3]) > amplicons[amplicon][2] and values3[2] <= stats[2] and values3[2] >= stats[1]:
                        amplicons[amplicon][2] = int(values3[3])

                    if int(values3[3]) >= int(args.cov):
                        if values3[2] <= stats[2] and values3[2] >= stats[1]:
                            amplicons[amplicon][0] += 1
                    if values3[2] <= stats[2] and values3[2] >= stats[1]:
                        amplicons[amplicon][3] += 1

        outfile = open(dir + "amplicons/" + name + "_" + values[3] + "_amplicons_stats", 'wb')
        for key, value in amplicons.items():
            stats = key.split("_")
            perc = str(round(float(value[0]) / value[3], 2))

            outfile.write(stats[0] + "\t" + stats[1] + "\t" + stats[2] + "\t" + stats[3] + "\t" + str(value[0]) + "\t" +
                          perc + "\t" + str(value[1]) + "\t" + str(value[2]) + "\n")

        outfile.close()

    sample = name
    type = "BRCA"

    os.system("cat " + dir + "amplicons/" + name + "*stats > " + dir + "amplicons/" + name + "_amplicons_stats")

    if name == args.name:
        os.system("echo '#Amplicon\\tChr\\tAmp_Start\\tAmp_End\\t#Amplicon_Nts_Covered (" + sample + ")\\t%Amplicon_Covered "
                      "(" + sample + ")\\tMin_Amp_Depth (" + sample + ")\\tMax_Amp_Depth (" + sample + ")' "
                      "> " + dir + "amplicons/" + "def/" + name + "_amplicons_stats_tmp")

        os.system("bedtools sort -i " + dir + "amplicons/" + name + "_amplicons_stats |"
                " awk '{print $4" + '"\\t"$1"\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t"$7"\\t"$8}' + "' >> " +
                  dir + "amplicons/" + "def/" + name + "_amplicons_stats_tmp")

    else:
        os.system(
            "echo '#Amplicon_Nts_Covered (" + sample + ")\\t%Amplicon_Covered "
                       "(" + sample + ")\\tMin_Amp_Depth (" + sample + ")\\tMax_Amp_Depth (" + sample + ")' "
                       "> " + dir + "amplicons/" + "def/" + name + "_amplicons_stats_tmp")

        os.system("bedtools sort -i " + dir + "amplicons/" + name + "_amplicons_stats |"
                                                     " awk '{print $5" + '"\\t"$6"\\t"$7"\\t"$8}' + "' >> " +
                  dir + "amplicons/" + "def/" + name + "_amplicons_stats_tmp")

    os.system("mv " + dir + "amplicons/" + "def/" + name + "_amplicons_stats_tmp " + dir + "amplicons/" + "def/" + name + "_amplicons_stats")

    i += 1


files = []
for file in os.listdir(dir + "amplicons/def/"):
    if file.endswith("_amplicons_stats"):
        files.append(os.path.join(dir + "amplicons/def/", file))

name = dir.split("/")[-3]
print name

os.system("paste -d '\\t' " + " ".join(
    sorted(files)) + " > " + dir + "amplicons/def/" + name + "_" + type + "_stats_by_amplicon.tsv")
os.system("rm " + dir + "amplicons/" + "*")

os.system("cp " + dir + "amplicons/def/" + name + "_" + type + "_stats_by_amplicon.tsv " + dir + "QC/mapping")

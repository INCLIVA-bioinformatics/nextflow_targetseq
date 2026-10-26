import os
import subprocess
import argparse

# Several arguments are needed to run the python script
parser = argparse.ArgumentParser()
# Arguments are added to the parser
parser.add_argument("file", help="VCF input file", type=str)
parser.add_argument("outfile", help="VCF conservation annotated ouput file", type=str)
parser.add_argument("conservation_file", help="Panel conservation file", type=str)
args = parser.parse_args()

file = open(args.file, 'rb')
f = file.readlines()
file.close()

outfile = open(args.outfile, 'wb')

for line in f:
    line = line.rstrip()
    if line.startswith("#"):

        if line.startswith("##INFO=<ID=ANN,Number"):
            #outfile.write("##CON=The conservation score for this site" + "\n")
            outfile.write("##INFO=<ID=Conservation,Number=A,Type=String,Description=" + '"##CON=The conservation score for this site">' + "\n")
            #line2 = line[:-2] + "|Conservation" + '">' + "\n"
            outfile.write(line + "\n")
            continue

        outfile.write(line + "\n")
        continue

    values = line.split("\t")
    info = values[7].split(";")
    for inf in info:
        if inf.startswith("ANN"):

            pos = values[1]
            ref = values[3]
            if len(ref) == 1:
                try:
                    output = subprocess.check_output("LC_ALL=C grep -F -w " + pos + " " + args.conservation_file, shell=True)
                except:
                    output = "0\n"

                out = output.split("\n")

                for l in out:
                    if l == "":
                        continue
                    stats = l.split("\t")

                    if stats[0] == values[0] and stats[2] == values[1]:
                        con = str(round(float(stats[3]), 3))
                        info = values[7].split(";")
                        for k in range(len(info)):
                            if info[k].startswith("ANN="):
                                ann = "Conservation=" + con + ";" + info[k]
                                if len(info) > k+1:
                                    info2 = ";".join(info[0:k]) + ";" + ann + ";" + ";".join(info[k+1:])
                                else:
                                    info2 = ";".join(info[0:k]) + ";" + ann
                                line2 = "\t".join(values[:7]) + "\t" + info2 + "\t" + "\t".join(values[8:])
                                outfile.write(line2 + "\n")

            else:
                con_t = ""
                dic_cont = {}
                for i in range(0, len(ref)):
                    posi = str(int(pos) + i)

                    try:
                        output = subprocess.check_output("LC_ALL=C grep -F -w " + pos + " " + args.conservation_file, shell=True)
                    except:
                        output = "0\n"


                    out = output.split("\n")

                    for l in out:
                        if l == "":
                            continue
                        stats = l.split("\t")

                        if stats[0] == values[0] and stats[2] == posi:

                            info = values[7].split(";")
                            for k in range(len(info)):
                                if info[k].startswith("ANN="):
                                    con = str(round(float(stats[3]), 3))
                                    con_t = con_t + "&" + con
                                    id = values[0] + "_" + values[1] + "_" + values[3]

                                    if id not in dic_cont:
                                        dic_cont[id] = [con]
                                    else:
                                        dic_cont[id].append(con)


                info = values[7].split(";")
                for k in range(len(info)):
                    if info[k].startswith("ANN="):
                        try:
                            con = "&".join(dic_cont[values[0] + "_" + values[1] + "_" + values[3]])
                        except:
                            con = "-"
                        ann = "Conservation=" + con + ";" + info[k]
                        if len(info) > k + 1:
                            info2 = ";".join(info[0:k]) + ";" + ann + ";" + ";".join(info[k + 1:])
                        else:
                            info2 = ";".join(info[0:k]) + ";" + ann
                        line2 = "\t".join(values[:7]) + "\t" + info2 + "\t" + "\t".join(values[8:])
                        outfile.write(line2 + "\n")

outfile.close()

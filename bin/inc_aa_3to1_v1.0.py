import argparse
import os

# Several arguments are needed to run the python script
parser = argparse.ArgumentParser()
# Arguments are added to the parser
parser.add_argument("input_csv", help="Input CSV file path", type=str)
args = parser.parse_args()

csv = args.input_csv


def replace(text, replacement):
    new = ''

    for c in text:
        try:
            new += replacement[c]
        except KeyError:
            new += c

    return new

aa_code = {
    "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F",
    "Gly": "G", "His": "H", "Ile": "I", "Lys": "K", "Leu": "L",
    "Met": "M", "Asn": "N", "Pro": "P", "Gln": "Q", "Arg": "R",
    "Ser": "S", "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y",
    "Ter": "*"
}

f = open(csv, 'rb')

outdir = csv.split("/")
outdir = "/".join(outdir[:-1]) + "/"

out = open(outdir + "out.csv", 'wb')

file_content = f.readlines()
f.close()

x = 1
for line in file_content:
    line = line.strip()
    values = line.decode('ascii').split(',')

    if x != 1:

        if (values[10] != "-"):

            aa_c = values[10]


            for i in range(len(aa_c)):
                if aa_c[i].isupper():
                    aa1 = aa_c[i] + aa_c[i+1] + aa_c[i+2]

                    aa2 = aa_code[aa1]
                    values[10] = (values[10].replace(aa1, aa2))
        new = values[10]

        line2 = ""
        for n in range(len(values)):
            if n == 10:
                line2 = line2 + new + ","
            else:
                line2 = line2 + values[n] + ","
        out.write(line2[:-1] + "\n")
    else:

        out.write(line + "\n")

    x = x + 1

out.close()

os.system("rm " + csv)
os.system("mv " + outdir + "out.csv " + csv)
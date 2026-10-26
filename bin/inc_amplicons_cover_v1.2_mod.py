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
        print ("ERROR: The entry directory path is wrong," \
            " please check the configuration file")

def list_samples(run_path_analysis):
    '''
    This function create a list of samples to process from a statistics file
    :return: List of samples to process
    '''
    for file_stats in os.listdir(outdir):
        if file_stats.endswith("_mapping_stats.csv"):
            file = open(outdir + file_stats, 'rb')

            f = file.readlines()
            file.close()
            f = f[1:]
            sample_list = []
            for line in f:
                line = line.strip()
                values = line.split(',')
                sample_list.append(values[0])
    return sample_list

def ampl_depth(outdir, sample, panel, type, bam):
    '''
    This function creates a depth statistics file for each amplicon from an specific panel
    '''
    os.system("samtools depth -d 9898989898 -a " + bam + " -b " +
              panel + " | awk '{print $1" + '"\\t"$2-1"\\t"$2"\\t"$3}' + "' > " + outdir + sample + "_" + type + "_amplicons_depth")

def create_dic(panel):
    '''
    This function creates an empty dictionary where each amplicon is a key of it
    '''
    file = open(panel, 'rb')
    f = file.readlines()
    file.close()

    amplicons = {}
    for line in f:
        line = line.rstrip()
        values = line.split("\t")
        amplicons[values[0] + "_" + values[1] + "_" + values[2] + "_" + values[3]] = [0, 9898989898, 0, 0]

    return amplicons

def fill_dic(amplicons, outdir, sample, type, cov):
    '''
    This function fills an empty dictionary with depth statistics (as values) for each amplicon from the panel (keys)
    '''
    file = open(outdir + sample + "_" + type + "_amplicons_depth", 'rb')
    f = file.readlines()
    file.close()
    
    for line in f:

        line = line.rstrip()
        values = line.split("\t")
        if type == "EPID":
            if values[0] == "chr11" and values[1] == '82398144' and values[2] == '82398145':
                continue
        for amplicon in amplicons.keys():
            stats = amplicon.split("_")
            if values[0] == stats[0]:
                if int(values[3]) < amplicons[amplicon][1] and values[2] <= stats[2] and values[2] >= stats[1]:
                    amplicons[amplicon][1] = int(values[3])

                if int(values[3]) > amplicons[amplicon][2] and values[2] <= stats[2] and values[2] >= stats[1]:
                    amplicons[amplicon][2] = int(values[3])

                if int(values[3]) >= int(cov):
                    if values[2] <= stats[2] and values[2] >= stats[1]:
                        amplicons[amplicon][0] += 1
                if values[2] <= stats[2] and values[2] >= stats[1]:
                    amplicons[amplicon][3] += 1
    
    return amplicons

def write_stats(outdir, sample, amplicons, type):
    '''
    This function creates a file for each sample with its amplicon depth statistics
    '''
    outfile = open(outdir + sample + "_" + type + "_amplicons_stats", 'wb')
    for key, value in amplicons.items():
        stats = key.split("_")
        try:
            perc = str(round(float(value[0]) / value[3], 2))
        except:
            perc = "0"

        outfile.write(stats[0] + "\t" + stats[1] + "\t" + stats[2] + "\t" + stats[3] + "\t" + str(value[0]) + "\t" +
                      perc + "\t" + str(value[1]) + "\t" + str(value[2]) + "\n")

    outfile.close()

def sort_stats(outdir, sample, type, i):
    '''
    This function sort and order an individual sample amplicon depth statistics
    '''
    if i == 0:
        os.system("echo '#Amplicon\\tChr\\tAmp_Start\\tAmp_End\\t#Amplicon_Nts_Covered (" + sample + ")\\t%Amplicon_Covered "
              "(" + sample + ")\\tMin_Amp_Depth (" + sample + ")\\tMax_Amp_Depth (" + sample + ")' "
              "> " + outdir + sample + "_" + type + "_tmp")

        os.system("bedtools  sort -i " + outdir + sample + "_" + type + "_amplicons_stats |"
              " awk '{print $4" + '"\\t"$1"\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t"$7"\\t"$8}' + "' >> " +
                  outdir + sample + "_" + type + "_tmp")
    else:
        os.system(
            "echo '#Amplicon_Nts_Covered (" + sample + ")\\t%Amplicon_Covered (" + sample + ")\\tMin_Amp_Depth (" +
            sample + ")\\tMax_Amp_Depth (" + sample + ")' > " + outdir + sample + "_" + type + "_tmp")
        os.system("bedtools sort -i " + outdir + sample + "_" + type + "_amplicons_stats |"
                                                                                     " awk '{print $5" + '"\\t"$6"\\t"$7"\\t"$8}' + "' >> " +
                  outdir + sample + "_" + type + "_tmp")

    os.system("mv " + outdir + sample + "_" + type + "_tmp " + outdir + sample + "_" + type + "_amplicons_stats")

def global_stats(outdir, type, name):
    '''
    This function creates a global file for all samples with its amplicon depth statistics from each individual
    sample file
    '''
    files = []
    for file in os.listdir(outdir):
        if file.endswith(type + "_amplicons_stats"):
            files.append(os.path.join(outdir, file))

    #name = run_path_analysis.split("/")[3]

    os.system("paste -d '\\t' " + " ".join(
        sorted(files)) + " > " + outdir + name + "_" + type + "_stats_by_amplicon.tsv")
    os.system("rm " + outdir + "*amplicons*")


def main():
    # Several arguments are needed to run the python script
    parser = argparse.ArgumentParser()
    # Arguments are added to the parser
    parser.add_argument("outdir", help="Out dir path", type=str)
    parser.add_argument("amplicons_panel_path", help="Panel path", type=str)
    parser.add_argument("cov", help="Minimum coverage value to consider a region covered", type=str)
    parser.add_argument("bam_list", help="BAM files names (comma separated)", type=str)
    parser.add_argument("name", help="Run name", type=str)
    # Nextflow: nombre del panel (en rutina se deduce de la ruta /nfs/home/panel_designs/<panel>/...)
    parser.add_argument("--panel_name", help="Panel name", type=str, default=None)
    args = parser.parse_args()


    #run_path_analysis = args.run_path_analysis
    #run_path_analysis = check_dir(run_path_analysis)
    #dir_exists(run_path_analysis)
    outdir = args.outdir
    outdir = check_dir(outdir)
    dir_exists(outdir)

    panel = args.amplicons_panel_path
    type = args.panel_name if args.panel_name else panel.split("/")[4]

    #sample_list = list_samples(run_path_analysis)
    bam_list = args.bam_list.split(",")

    i = 0
    for bam in bam_list:
        sample = bam.split("/")[-1]
        sample = sample.split(".")[0]

        ampl_depth(outdir, sample, panel, type, bam)

        amplicons = create_dic(panel)
        amplicons = fill_dic(amplicons, outdir, sample, type, args.cov)

        write_stats(outdir, sample, amplicons, type)

        sort_stats(outdir, sample, type,  i)

        i += 1

    global_stats(outdir, type, args.name)

if __name__ == '__main__':
    main()

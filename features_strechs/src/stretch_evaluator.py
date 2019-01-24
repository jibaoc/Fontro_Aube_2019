#!/user/bin/python3

"""
The goal of this script is to make a barplot indicating the proportion of low complexity sequenceq different \
kind of low complexity sequences for a nucleotide in the different sequence given in input.
It will allow to decide what type of low complexity sequence to use to better discriminate the input files
"""


# Set the environment :
from Bio import SeqIO
import pandas as pd
from math import sqrt
from ncephes import cprob
from matplotlib import pyplot as plt
import os
import sys
import argparse
import matplotlib.patches as mpatches


codon2aminoAcid = dict(TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L", ATT="I", ATC="I",
                       ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V", TCT="S", TCC="S", TCA="S", TCG="S",
                       CCT="P", CCC="P", CCA="P", CCG="P", ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A",
                       GCA="A", GCG="A", TAT="Y", TAC="Y", TAA="", TAG="", CAT="H", CAC="H", CAA="Q", CAG="Q",
                       AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E", TGT="C", TGC="C",
                       TGA="", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R", AGT="S", AGC="S", AGA="R", AGG="R",
                       GGT="G", GGC="G", GGA="G", GGG="G")

feature_dic = {
    "Very-small": ["A", "C", "G", "S"],
    "Small#2": ["A", "C", "D", "G", "N", "P", "S", "T"],
    "Large": ["F", "I", "K", "L", "M", "R", "W", "Y"],
    "Disorder-promoting#1": ["A", "E", "G", "K", "P", "Q", "R", "S"],
    "Order-promoting#1": ["C", "F", "I", "L", "N", "W", "V", "Y"],
    "Disorder-promoting#2": ["A", "E", "G", "K", "P", "Q", "S"],
    "Order-promoting#2": ["C", "F", "H", "I", "L", "M", "N", "W", "V", "Y"],
    "Polar-uncharged#1": ["C", "N", "Q", "S", "T", "Y"],
    "Polar-uncharged#2": ["N", "Q", "S", "T", "Y"],
    "Charged#1": ["R", "H", "K", "D", "E"],
    "Charged#2": ["R", "K", "D", "E"],
    "Hydrophilic#1": ["D", "E", "K", "N", "Q", "R"],
    "Hydrophobic#1": ["A", "C", "F", "I", "L", "M", "V"],
    "Neutral": ["G", "H", "P", "S", "T", "Y"],
    "Hydroxylic": ["S", "T", "Y"],
    "Negatively-charged": ["D", "E"],
    "Positively-charged#1": ["R", "H", "K"],
    "Positively-charged#2": ["R", "K"]
}

###
# Functions
###


def read_fasta_files(fasta_file, stretch_len):
    """
    :param fasta_file: (string) the name of a fasta file of interest
    :param stretch_len: (int) the len of the low complexity sequence of interest
    :return: the list of sequence in the fasta file, nt and amino
    """
    list_sequence = []
    nt_seq = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_aa = translator(record.seq)
        if len(seq_aa) > stretch_len-1:
            list_sequence.append(seq_aa)
        if len(record.seq) > stretch_len-1:
            nt_seq.append(record.seq)
    return nt_seq, list_sequence


def read_sequence(excel_file, stretch_len):
    """
    :param stretch_len: (int) the len of the low complexity sequence of interest
    :param excel_file: (string) path to a query file given by the tRNA program
    :return: (2 list of strings) the list of cds sequence and the list of amino acid sequence
    """
    # opening the excel file
    xl = pd.ExcelFile(excel_file)
    df = "NA"
    # opening the sheet of interest
    for sheet in xl.sheet_names:
        if "sequence" == sheet:
            df = xl.parse(sheet)
    # if the sheet "exon_skipping*" doesn't exist, the faRLine file cannot be used so we stop the program
    if str(df) == "NA":
        print("the sheet names sequence wasn't found")
        print("exiting...")
        exit(1)

    cds = []
    pep = []
    for row in df.itertuples():
        if not isinstance(row.CDS_genomic_sequence, float):
            if len(row.CDS_genomic_sequence) > stretch_len-1:
                cds.append(row.CDS_genomic_sequence)
        if not isinstance(row.CDS_peptide_sequence, float):
            if len(row.CDS_peptide_sequence) > stretch_len-1:
                pep.append(row.CDS_peptide_sequence)
    return cds, pep


def translator(seq):
    """
    :param seq: (string) a nucleotide sequence
    :return: the translated sequence
    """
    res = ""
    for i in range(0, len(seq)-2, 3):
        res += codon2aminoAcid[seq[i:i+3]]
    return res


def stretch_finder_feature(sequence, feature, stretch_len, stretch_content):
    """
    :param sequence: (string) an amino acid sequences
    :param feature: (string) the name of the feature to return
    :param stretch_len: (int) the length of the low complexity sequence of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there is a low complexity sequence in the sub-sequence
    :return: the number of low complexity sequences of the feature "feature" here
    """
    nb_stretch = 0
    for i in range(len(sequence) - stretch_len - 1):
        count = 0
        for letter in sequence[i:i+stretch_len]:
            if letter in feature_dic[feature]:
                count += 1
        if count >= stretch_content:
            nb_stretch += 1
    return nb_stretch


def stretch_finder_aa(sequence, aa, stretch_len, stretch_content):
    """
    :param sequence: (string) an amino acid sequences
    :param aa: (string) the name of the aa to return
    :param stretch_len: (int) the length of the low complexity sequence of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there is a low complexity sequence in the sub-sequence
    :return: the number of low complexity sequence of the aa "aa" here
    """
    nb_stretch = 0
    for i in range(len(sequence) - stretch_len - 1):
        count = 0
        for letter in sequence[i:i+stretch_len]:
            if letter == aa:
                count += 1
        if count >= stretch_content:
            nb_stretch += 1
    return nb_stretch


def stretch_finder_nt(sequence, nt, stretch_len, stretch_content):
    """
    :param sequence: (string) an amino acid sequences
    :param nt: (string) the name of the nt to return
    :param stretch_len: (int) the length of the low complexity sequence of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there is a low complexity sequence in the sub-sequence
    :return: the number of low complexity sequences of the nucleotide "nt" here
    """
    nb_stretch = 0
    if nt in ["A", "T", "G", "C"]:
        for i in range(len(sequence) - stretch_len - 1):
            count = 0
            for letter in sequence[i:i+stretch_len]:
                if letter == nt:
                    count += 1
            if count >= stretch_content:
                nb_stretch += 1
    else:
        iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
                 'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}
        for i in range(len(sequence) - stretch_len - 1):
            count = 0
            for letter in sequence[i:i+stretch_len]:
                if letter in iupac[nt]:
                    count += 1
            if count >= stretch_content:
                nb_stretch += 1
    return nb_stretch


def get_stretch(list_of_sequence, unit, unit_type, stretch_len, stretch_content):
    """

    :param list_of_sequence: (list of string) list of amino acid sequences
    :param unit: (string) he name of the unit of interest
    :param stretch_len: (int) the length of the low complexity sequence of interest
    :param unit_type: (string) the name of the unit of interest (aa, feature, nt)
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there is a low complexity sequence in the sub-sequence
    :return: (dictionary of int) the number of sequences having 0 to 10+ low complexity sequence\
     in the list of sequences
    """
    st_dic = {i: 0 for i in range(0, 11)}
    for sequence in list_of_sequence:
        if unit_type == "feature":
            nb_stretch = stretch_finder_feature(sequence, unit, stretch_len, stretch_content)
        elif unit_type == "aa":
            nb_stretch = stretch_finder_aa(sequence, unit, stretch_len, stretch_content)
        else:
            nb_stretch = stretch_finder_nt(sequence, unit, stretch_len, stretch_content)
        if nb_stretch > 10:
            nb_stretch = 10
        st_dic[nb_stretch] += 1
    counter = 0
    for key in st_dic.keys():
        counter += st_dic[key]
    st_dic["all"] = counter
    return st_dic


def frequency_test(obs1, tot1, obs2, tot2):
    """

    :param obs1: (int) the count number of an amino acid X in the set of protein 1.
    :param tot1: (int) the total number of amino acids in the set of protein 1.
    :param obs2: (int) the count number of an amino acid X in the set of protein 2.
    :param tot2: (int) the total number of amino acids in the set of protein 2.
    :return: proportion test p-value
    """
    if obs1 == 0 and obs2 == 0:
        return float('nan')
    mean1 = float(obs1) / tot1
    mean2 = float(obs2) / tot2

    var1 = float(obs1) * (1 - mean1) * (1 - mean1) + (tot1 - obs1) * mean1 * mean1
    var2 = float(obs2) * (1 - mean2) * (1 - mean2) + (tot2 - obs2) * mean2 * mean2

    df = tot1 + tot2 - 2
    svar = (var1 + var2) / df
    t = (mean1-mean2) / sqrt(svar*(1.0/tot1 + 1.0/tot2))

    return cprob.incbet(0.5*df, 0.5, df/(df+t*t))


def get_more(dic, num):
    """
    :param dic: dic of int) the number of exons having 0, 1, ..., 10 or more stretch in
    an exon set.
    :param num: (int) the number of exons having "num" low complexity sequence(s) or more
    :return: the prop of exon having "num" low complexity sequence(s) or more
    """
    count = 0
    for key in dic.keys():
        if key != "all":
            if key >= num:
                count += dic[key]
    return [float(count) / dic["all"], float(count)]


def barplot_maker(stretch_dics, name_stretch, name_file, unit_name, unit, output, project):
    """
    :param stretch_dics: (list of int dictionary) list of the number of exons in list of sequence
     having 1 to 10 stretches whose characteristics are described in "name_stretch"
    :param name_stretch: (list of string) the stretch of interest (i.e 6/5 : low complexity sequence\
     of length 6 having a content of unit_name (=feature, aa, nt) of unit(=A or C or G or T if unit_name=nt)
    :param name_file: (list of string) the name of the file where the list of exons analysed for low complexity \
    sequence content
    :param unit_name: (str) the name of the unit studied (nt, aa, feature
    :param unit: (str) the name of the unit studied
    :param output: (string) where the graph will be placed
    :param project: (string) the name of the project where the list of up and down exons are from
    """
    fig, ax = plt.subplots(figsize=(48. / 2.54 * 1.11, 27. * 0.91 / 2.54))

    # set the abscissa values for bars and tick_abscissa value for legend
    abscissa = []
    tick_abscissa = []
    i = 1
    for j in range(1, len(stretch_dics) + 1):
        abscissa.append(i)
        if j % 3 == 0:
            i += 2
        elif j % 3 == 2:
            tick_abscissa.append(i)
            i += 1
        else:
            i += 1

    if name_file[0] != "up exons":
        colors = ["blue", "red", "grey"] * (len(stretch_dics) / 3)
    else:
        colors = ["#FF9D09", "#8408B9", "grey"] * (len(stretch_dics) / 3)
    val = [get_more(my_dic, 1) for my_dic in stretch_dics]
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray', linestyle='dashed', alpha=0.3)
    ax.bar(abscissa,
           [get_more(my_dic, 1)[0] for my_dic in stretch_dics],
           width=0.8, color=colors)
    ax.set_xticks(tick_abscissa)
    ax.set_xticklabels(["stretches - " + str(cur_name) for cur_name in name_stretch])
    ax.set_title(u"Proportion of exons having more than 1\nstretches of " +
                 str(unit_name) + " " + str(unit) + " in CCE, UP and down set of exons\n"
                 + " Project : " + str(project))
    up_patch = mpatches.Patch(color=colors[0], label=name_file[0])
    down_patch = mpatches.Patch(color=colors[1], label=name_file[1])
    crtl_patch = mpatches.Patch(color=colors[2], label=name_file[2])
    for i in range(len(val)):
        ax.text(abscissa[i], val[i][0] + val[i][0] * 0.01,
                str(int(val[i][1])) + "/" + str(stretch_dics[i]["all"]),
                fontsize=10, horizontalalignment='center')

    plt.legend(handles=[up_patch, down_patch, crtl_patch])
    ax.set_ylabel(u"Proportion of exons")
    ax.set_xlabel(u"Number of stretches")
    fig.add_subplot(ax)
    plt.savefig(output + "stretches_exploration_" +
                str(unit_name) + "_" + str(unit) + "_" + project + "_figure.png", bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    return output + "stretches_exploration_" + str(unit_name) + "_" + str(unit) + "_" + project + "_figure.png"


def figure_creator(file_up, file_down, fasta, unit_type, output, project):
    """
    :param file_up: (string) the excel file containing the up exons sequences
    :param file_down: (string) the excel file containing the down exons sequences
    :param fasta: (string) a fasta file containing random sequences
    :param unit_type: (string) the name of the unit of interest (aa, feature, nt)
    :param output: (string) the path where the graphics will be created
    :param project: (string) the name of the project where the up and down or fasta file are from
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dir/")
    if unit_type == "feature":
        unit_list = feature_dic.keys()
    elif unit_type == "aa":
        unit_list = "ACDEFGHIKLMNPQRSTVWY"
    else:
        unit_list = "ACGTSWYRKM"
    if fasta:
        name_file = [file_up.replace(".fasta", "").split("/")[-1] + " sequences",
                     file_down.replace(".fasta", "").split("/")[-1] + "sequences", "CCE control exons"]
    else:
        name_file = ["up exons", "down exons", "CCE control exons"]
    for cur_unit in unit_list:
        st_dics = []
        name_stretch = []
        for stretch_content, stretch_len in [[4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12]]:
            name_stretch.append(str(stretch_content) + "/" + str(stretch_len))
            if fasta:
                up_nt_seq, up_aa_seq = read_fasta_files(file_up, stretch_len)
                down_nt_seq, down_aa_seq = read_fasta_files(file_down, stretch_len)
            else:
                up_nt_seq, up_aa_seq = read_sequence(file_up, stretch_len)
                down_nt_seq, down_aa_seq = read_sequence(file_down, stretch_len)
            if unit_type == "feature":
                up_seq = up_aa_seq
                down_seq = down_aa_seq
            elif unit_type == "aa":
                up_seq = up_aa_seq
                down_seq = down_aa_seq
            else:
                up_seq = up_nt_seq
                down_seq = down_nt_seq

            st_dics.append(get_stretch(up_seq, cur_unit, unit_type, stretch_len, stretch_content))
            st_dics.append(get_stretch(down_seq, cur_unit, unit_type, stretch_len, stretch_content))
            mod = __import__("CCE_" + str(unit_type) + "_" + str(stretch_len) + "_" + str(stretch_content) + "_dic")
            dep = cur_unit + "_dic"
            st_dics.append(eval("mod." + dep))
        barplot_maker(st_dics, name_stretch, name_file, unit_type, cur_unit, output, project)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""From 2 files (up or down file, fasta or excel file)
                                     give graphic for every unit of type unit
                                     displaying the number of sequence having more than one
                                     low complexity sequences in the two givens file for different kind of 
                                     low complexity sequences:
                                     i.e. different lenght and different content
    """)
    # Arguments for the parser

    parser.add_argument('--unit_type', dest='unit_type', help="the type of tunit",
                        default="nt")
    req_arg = parser.add_argument_group("required arguments")
    req_arg.add_argument('--output', dest='output', help="An output folder",
                         required=True)
    parser.add_argument('--up', dest='up', help="an up fasta or query_result excel file",
                        default=None)
    parser.add_argument('--down', dest='down', help="a down fasta or query_result excel file",
                        default=None)
    req_arg.add_argument('--fasta', dest='fasta', help="a fasta file",
                         required=True)
    req_arg.add_argument('--project', dest='project', help="The name of the project",
                         required=True)

    args = parser.parse_args()
    if args.output[-1] != "/":
        args.output += "/"
    if not os.path.isdir(args.output):
        print("ERROR : the given output folder does not exist")
        exit(1)

    if args.unit_type not in ["feature", "aa", "nt"]:
        print("ERROR : wrong value for the parameter unit_type...")
        print("It must be 'feature', 'nt' or 'aa'")
        exit(1)

    if args.up is not None or args.up is not None:
        if not os.path.isfile(args.up):
            print("Error : the up file does not exit")

        if not os.path.isfile(args.down):
            print("Error : the down file does not exit")

    if args.fasta == "True":
        args.fasta = True

    elif args.fasta == "False":
        args.fasta = False
    else:
        print("Wrong value for the fasta argument ! (True or False")
        print("Exiting...")
        exit(1)

    figure_creator(args.up, args.down, args.fasta, args.unit_type, args.output, args.project)


if __name__ == "__main__":
    launcher()

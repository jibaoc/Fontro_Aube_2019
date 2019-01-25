"""Description: script that can generate a fasta file with random sequences.
It will use a given probability to have a certain protein feature at each position of
an amino acid sequence."""


#############################
#          IMPORTS          #
#############################

import argparse
import random
import sys
import os
from dicitonary import *
import numpy as np

feature_dic = {
    "Very-small": ["A", "C", "G", "S"],
    "Small#2": ["A", "C", "D", "G", "N", "P", "S", "T"],
    "Large" : ["F", "I", "K", "L", "M", "R", "W", "Y"],
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


def ctrl_dic_adapter(dic):
    """
    Adaptation of a dictionary of counts.

    Description:
    1. Turn first the a dictionary of count to a dictionary
        of proportions.
        {A:1, B:7, C:2, all:10} => {A:0.1, B:0.7, C:0.2, all:1.0}
    2.Then give for each key the cumulative proportion giving their
        order.
        {A:0.1, B:0.7, C:0.2, all:1.0} => {A:0.1, B:0.8, C:1.0}
    3. Turn each value into an interval (list of value): the
    up border of the interval was the old value for each key except
     for the key with the proportion 1.0 (it will be 1 + 0.00000001).
     The down border will be 0 for the first key and the up border
     of the previous key for the other keys.
    {A:0.1, B:0.8, C:1.0} => {A:[0, 0.1], B:[0.1, 0.8], C:[0.8, 1.00000001]}
    :param dic: (dictionary of int)
    :return: the adapted dictionary
    """
    res_dic = {}
    for key in dic.keys():
        if key != "all":
            res_dic[key] = float(dic[key]) / dic["all"]
    list_key = res_dic.keys()
    tmp = 0.
    for key in list_key:
        res_dic[key] = res_dic[key] + tmp
        tmp = res_dic[key]
    sorted_res = sorted(res_dic.items(), key=lambda l: l[1], reverse=False)
    interval_dic = {}
    for i in range(len(sorted_res)):
        if i == 0:
            interval_dic[sorted_res[i][0]] = \
                [0, sorted_res[i][1]]
        elif i < len(sorted_res)-1:
            interval_dic[sorted_res[i][0]] = \
                [sorted_res[i-1][1], sorted_res[i][1]]
        else:
            interval_dic[sorted_res[i][0]] = \
                [sorted_res[i-1][1], sorted_res[i][1] + 0.00000001]
    return interval_dic


def generate_dic(ctrl_dic, unit_list):
    """

    :param ctrl_dic: (dictionary of int) for each unit (aa or codon) give it's count in the control set of exons
    :param unit_list: (list of string) list of aa or codon
    :return:
    """
    aa_dic = {}
    count = 0
    for aa in unit_list:
        aa_dic[aa] = ctrl_dic[aa]
    for aa in unit_list:
        count += aa_dic[aa]
    aa_dic['all'] = count
    res_aa = ctrl_dic_adapter(aa_dic)
    return res_aa


def get_cur_val(ctr_dic, value):
    """
    :param ctr_dic: (dictionary of list of 2 float)
    the keys of the dictionary are the codons, and their are link to
    an interval of 2 values.
    :param value: (int) a value
    :return: the key linked to an interval in (ctr_dic)
    containing the value "value"
    """
    for key in ctr_dic.keys():
        if ctr_dic[key][0] <= value < ctr_dic[key][1]:
            return key
    return None


def sequence_generator(length, prop, feature, ctrl):
    """

    :param length: (int) the length of the sequence to generate
    :param prop: (int) the percentage of the feature 'feature' in the sequence to create
    :param feature: (string) the name of the feature we want to analyse
    :param ctrl: (string) the control used for amino acid and codon determination
    :return:
    - fseq (string) a nucleotide sequence
    - ap : (float) the Adenine proportion in fseq
    - cp : (float) the Cytosine proportion in fseq
    - gp : (float) the Guanine proportion in fseq
    - tp : (float) the Thymine proportiob in  fseq
    - ftp : (float) the proportion of the feature 'feature' in fseq
    """
    if ctrl != "RD":
        file_dir = os.path.dirname(os.path.realpath(__file__))
        sys.path.insert(0, file_dir + "/control_dic/")
        mod = __import__(ctrl + "_dic")

    nbr_f = (prop + (np.random.randn() / 30)) * 100
    if nbr_f < 0:
        nbr_f = 0
    if nbr_f > 100:
        nbr_f = 100
    seq = ""
    ft_chooser = "F" * int(round(nbr_f)) + "N" * int(100 - (prop * 100))
    ft_count = 0
    for i in range(length/3):
        ft = ft_chooser[random.randint(0, len(ft_chooser)-1)]
        amino_acid_list = feature_dic[feature]
        if ft == "N":
            aa_list = "ACDEFGHIKLMNPQRSTVWY"
            new_seq = []
            for aa in aa_list:
                if aa not in amino_acid_list:
                    new_seq.append(aa)
            amino_acid_list = new_seq
        else:
            ft_count += 1
        if ctrl != "RD":
            aa_dic = generate_dic(mod.da, amino_acid_list)
            aa_chosen = get_cur_val(aa_dic, random.random())
            codon_list = amino_acid2codon[aa_chosen].split(",")
            codon_dic = generate_dic(mod.dc, codon_list)
            codon_chosen = get_cur_val(codon_dic, random.random())
        else:
            aa_chosen = amino_acid_list[random.randint(0, len(amino_acid_list)-1)]
            codon_list = amino_acid2codon[aa_chosen].split(",")
            codon_chosen = codon_list[random.randint(0, len(codon_list)-1)]
        seq += codon_chosen

    fseq = ""
    i = 0
    while i < len(seq):
        fseq += seq[i:i + 70] + "\n"
        i += 70

    ap = round(float(seq.count("A")) / len(seq), 2)
    cp = round(float(seq.count("C")) / len(seq), 2)
    gp = round(float(seq.count("G")) / len(seq), 2)
    tp = round(float(seq.count("T")) / len(seq), 2)
    ftp = round(float(ft_count) / (len(seq)/3), 2)
    return fseq, ap, cp, gp, tp, ftp


def adjusting_count_dic(dic_aa, prop, feature):
    """
    This function adjust the proportion of amino acid in groups given in the "feature" list so
    this group of amino acid reaches the proportion given in "prop" list.
    :param dic_aa: dictionary of count of amino acid
    :param prop: (list of int) list of proportion for the feature given in feature list
    :param feature: (list of string) list of feature
    :return: the dictionary count for amino acid adjusted
    """
    multplication_aa = {aa: [] for aa in "ACDEFGHIKLMNPQRSTVWY"}
    for i in range(len(feature)):
        count = 0
        for aa in feature_dic[feature[i]]:
            count += dic_aa[aa]
        for aa in feature_dic[feature[i]]:
            multplication_aa[aa].append(float(prop[i] * dic_aa["all"]) / count)
    for aa in multplication_aa.keys():
        if len(multplication_aa[aa]) > 0:
            multplication_aa[aa] = sum(multplication_aa[aa])/len(multplication_aa[aa])
        else:
            multplication_aa[aa] = 1
    for aa in multplication_aa.keys():
        dic_aa[aa] = (round(dic_aa[aa] * multplication_aa[aa], 0))
    count = 0
    count2 = 0
    for aa in multplication_aa.keys():
        if multplication_aa[aa] != 1:
            count += dic_aa[aa]
        else:
            count2 += dic_aa[aa]
    for aa in multplication_aa.keys():
        if multplication_aa[aa] == 1:
            if count2 != 0:
                multplication_aa[aa] = float(max(dic_aa["all"] - count, 0)) / count2
                dic_aa[aa] = (max(round(dic_aa[aa] * (dic_aa["all"] - count) / count2, 0), 0))
            else:
                multplication_aa[aa] = 0
                dic_aa[aa] = 0

    count = 0
    for aa in multplication_aa.keys():
        count += dic_aa[aa]
    dic_aa["all"] = count
    return dic_aa


def translator(seq):
    """

    :param seq: (list of string) list of codon
    :return: list of codon translated : i.e. amino acid list
    """
    aa = ""
    for i in range(0, len(seq)-2, 3):
        aa += codon2aminoAcid[seq[i:i+3]]
    return aa


def feature_frequency_calculator(seq, feature):
    """
    :param seq: (list of string) list of codons
    :param feature: (string) the feature of interest
    :return: the frequency of the feature in sequence 'seq'
    """
    aa_seq = translator(seq)
    count = 0
    for aa in feature_dic[feature]:
        count += aa_seq.count(aa)
    return float(count) / len(aa_seq)


def sequence_generator_with_multiple_feature(dic_aa, length, feature, ctrl):
    """
    :param dic_aa: dictionary of count of amino acid
    :param length: (int) the length of the sequence we will create
    :param feature: (list of string) list of feature
    :param ctrl: (string) the control exons used to generate the control sequence
    :return:
        - fseq : (string) the nucleotide sequence completed
        - ft_freq : (list of float) the proportion of feature (given in the same order as in "feature" list)
         in fseq
    """

    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")

    if ctrl == "RD":
        count = dic_aa["all"]
        dic_aa = {aa: int(count / 20) for aa in "ACDEFGHIKLMNPQRSTVWY"}
        count = 0
        for aa in dic_aa.keys():
            count += dic_aa[aa]
        dic_aa["all"] = count
    seq = ""
    for i in range(length/3):
        if ctrl != "RD":
            aa_dic = generate_dic(dic_aa, list("ACDEFGHIKLMNPQRSTVWY"))
            aa_chosen = get_cur_val(aa_dic, random.random())
            codon_list = amino_acid2codon[aa_chosen].split(",")
            codon_dic = generate_dic(mod.dc, codon_list)
            codon_chosen = get_cur_val(codon_dic, random.random())
        else:
            aa_dic = generate_dic(dic_aa, list("ACDEFGHIKLMNPQRSTVWY"))
            aa_chosen = get_cur_val(aa_dic, random.random())
            codon_list = amino_acid2codon[aa_chosen].split(",")
            codon_chosen = codon_list[random.randint(0, len(codon_list)-1)]
        seq += codon_chosen

    ft_freq = []
    for i in range(len(feature)):
        ft_freq.append(feature_frequency_calculator(seq, feature[i]))
    fseq = ""
    i = 0
    while i < len(seq):
        fseq += seq[i:i + 70] + "\n"
        i += 70

    return fseq, ft_freq


def header_generator(length, a_prop, c_prop, g_prop, t_prop, ft_prop, feature, num_seq):
    """
    :param length: (int) the length of the sequence we will create
    :param a_prop: (float) proportion of A in the sequence we will create
    :param t_prop: (float) proportion of t in the sequence we will create
    :param c_prop: (float) proportion of c in the sequence we will create
    :param g_prop: (float) proportion of g in the sequence we will create
    :param ft_prop: (float) proportion of the feature in the sequence
    :param feature: (string) the feature for which we want to control the frequency in the random sequences
    :param num_seq: (int) the number of sequence we will create
    :return: (string) the header of a sequence
    """
    header = ">seq" + str(num_seq) + " | " + feature + " : " + str(ft_prop) + " - A : " + str(a_prop) + \
             " - C : " + str(c_prop) + " - G : " + str(g_prop) + \
             " - T : " + str(t_prop) + " | length : " + str(length)
    return header


def header_generator_multiple_ft(length, ft_freq, features, num_seq):
    """

    :param length: (int) the length of the sequence we will create
    :param ft_freq: list of float) the proportion of feature (given in the same order as in "feature" list)
    :param features: (list of string) list of feature which we want to change the frequency
    :param num_seq: (int) the num of the current sequence we want to create
    :return: (string) an header for the sequence of interest
    """
    hd = ""
    for i in range(len(features)):
        hd += str(features[i]) + " : " + str(ft_freq[i]) + " - "

    header = ">seq" + str(num_seq) + " | " + hd + " | length : " + str(length)
    return header


def dictionnary_count_adaption(prop, feature, ctrl):
    """
    Adjustment of the CCE dictionary of count of amino acid : The iterative procedure is very useful
    when an amino acid is shared between to groups of feature which we want changed the frequency
    :param prop: (list of int) list of proportion for the feature given in feature list
    :param feature: (list of string) list of feature
    :param ctrl: (string) the control exons used to generate the control sequence
    :return: dictionary of count of amino acid
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")
    dic_aa = mod.da

    for i in range(300):
        dic_aa = adjusting_count_dic(dic_aa, prop, feature)
    return dic_aa


def fasta_generator(size_int, prop, feature, number_seq, output, out_name, ctrl):
    """
    Write a fasta file of number_seq sequences having a size in  'size_int' and proportion corresponding to a_prop,
    t_prop, c_prop, g_prop
    :param size_int: (list of 2 int) the min size possible and the max size possible of the sequences we want to create
    :param prop: (int) the proportion of the feature "feature" in the random sequence to generate
    :param feature: (string) the feature for which we want to control the frequency in the random sequences
    :param number_seq: (int) the number of sequence we will create
    :param output: (string) path where the fasta file will be created
    :param out_name: (string) the name of the fasta file to create
    :param ctrl: (string) the control exons used to generate the control sequence
    False if the proportion have to be equal or very close  to what the user specified
    """
    with open(output + out_name + ".fasta", "w") as outfile:
        list_ft_freq = []
        list_ft = [[] for i in range(len(feature))]
        if isinstance(feature, str):
            for i in range(1, number_seq+1):
                length = random.randint(size_int[0], size_int[1])
                seq, ap, cp, gp, tp, ftp = sequence_generator(length, prop, feature, ctrl)
                header = header_generator(len(seq), ap, cp, gp, tp, ftp, feature, i)
                list_ft_freq.append(ftp)
                outfile.write(header + "\n" + seq + "\n")
        else:
            dic_aa = dictionnary_count_adaption(prop, feature, ctrl)
            for i in range(1, number_seq + 1):
                length = random.randint(size_int[0], size_int[1])
                seq, ft_freq = sequence_generator_with_multiple_feature(dic_aa, length, feature, ctrl)
                header = header_generator_multiple_ft(length, ft_freq, feature, i)
                for i in range(len(ft_freq)):
                    list_ft[i].append(round(ft_freq[i], 3))
                outfile.write(header + "\n" + seq + "\n")

    if isinstance(feature, str):
        mean = 0
        for val in list_ft_freq:
            mean += val
        mean = float(mean)/len(list_ft_freq)
        print("frequence of " + feature + " amino acids in the file : " + str(mean))
    else:
        for i in range(len(list_ft)):
            print("frequence of " + feature[i] + " amino acids in the file : " + str(float(sum(list_ft[i]))/number_seq))


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""Given a number of sequence N and
                                       a feature and it's proportion,
                                       create a fasta file of N random sequences with proportions this feature
                                       specified above.

    Here are the accepted features :

    "Small#1": ["A", "C", "D", "G", "N", "P", "S", "T", "V"],
    "Small#2": ["A", "C", "D", "G", "N", "P", "S", "T"],
    "Large" : ["F", "I", "K", "L", "M", "R", "W", "Y"],
    "Disorder_promoting#1": ["A", "E", "G", "K", "P", "Q", "R", "S"],
    "Order_promoting#1": ["C", "F", "I", "L", "N", "W", "V", "Y"],
    "Disorder_promoting#2": ["A", "E", "G", "K", "P", "Q", "S"],
    "Order_promoting#2": ["C", "F", "H", "I", "L", "M", "N", "W", "V", "Y"],
    "Polar_uncharged#1": ["C", "N", "Q", "S", "T", "Y"],
    "Polar_uncharged#2": ["N", "Q", "S", "T", "Y"],
    "Charged": ["R", "H", "K", "D", "E"],
    "Hydrophilic#1": ["D", "E", "K", "N", "Q", "R"],
    "Hydrophobic#1": ["A", "C", "F", "I", "L", "M", "V"],
    "Hydrophilic#2": ["D", "E", "H", "K", "N", "Q", "R", "S", "T"],
    "Hydrophobic#2": ["A", "C", "F", "I", "L", "M", "P", "V", "W", "Y"],
    "Hydroxylic": ["S", "T", "Y"],
    "Negatively_charged": ["D", "E"],
    "Positively_charged": ["R", "H", "K"],

    """,
                                     usage='%(prog)s --feature a_feature_name --prob a_prob_name ')
    # Arguments for the parser

    parser.add_argument('--output', dest='output', help="An output folder",
                        default=".")
    parser.add_argument('--filename', dest='filename', help="the name of the fasta file, the program will produce ",
                        default="result")
    parser.add_argument('--size_inf', dest='size_inf', help="the smallest size possible of the sequences in the fasta",
                        default=50)
    parser.add_argument('--size_max', dest='size_max', help="the largest size possible of the sequences in the fasta",
                        default=300)
    parser.add_argument('--nbr_seq', dest='nbr_seq', help="the number of sequence in the fasta file",
                        default=300)
    parser.add_argument('--ctrl', dest='ctrl', help="the control used for fasta generation (CCE, ALL, ACE, RD)",
                        default="CCE")

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument('--feature', dest='feature', help="the feature tha will help for sequence generation",
                               required=True)
    required_args.add_argument('--prop', dest='prop', help="the prob to have the feature at each aa position (0-100)",
                               required=True)

    args = parser.parse_args()  # parsing arguments

    if args.filename == "result":
        args.filename = args.ctrl + "_" + args.feature + "_" + args.prop

    try:
        args.prop = float(args.prop)
        if 0 > args.prop or args.prop > 100:
            print("ERROR : wrong probability value")
            exit(1)
    except ValueError:
        try:
            args.prop = args.prop.split(",")
            for i in range(len(args.prop)):
                try:
                   args.prop[i] = int(args.prop[i])
                   if 0 < args.prop[i] > 100:
                       print("ERROR : wrong probability value")
                       exit(1)
                   args.prop[i] = float(args.prop[i]) / 100
                except ValueError:
                    print("ERROR : wrong probability value in the list")
                    exit(1)
        except ValueError:
            print("ERROR : wrong probability list")
            exit(1)


    if "," in args.feature:
        args.feature = args.feature.split(",")

    try:
        args.size_inf = int(args.size_inf)
        args.size_max = int(args.size_max)
        if args.size_inf > args.size_max:
            print("WARNING : maximum size value smaller than minimum size value")
            print("switching size value (min <=> max)")
            temp = args.size_inf
            args.size_inf = args.size_max
            args.size_max = temp
    except ValueError:
        print("ERROR : wrong size values.")
        print("Exiting...")
        exit(1)
    size_int = [args.size_inf, args.size_max]

    try:
        args.nbr_seq = int(args.nbr_seq)
        if args.nbr_seq < 0:
            int("a")
    except ValueError:
        print("ERROR : wrong 'nbr_seq' values")
        print("Exiting")
        exit(1)

    if not os.path.isdir(args.output):
        print("The given path in 'output' doesn't exist !")
        print("fasta file will be created in your current working directory")
        args.output = "./"

    if args.output[-1] != "/":
        args.output += "/"

    print(args.prop)
    fasta_generator(size_int, args.prop, args.feature, args.nbr_seq, args.output, args.filename, args.ctrl)

if __name__ == "__main__":
    launcher()

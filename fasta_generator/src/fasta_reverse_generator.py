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
from dicitonary import amino_acid2codon

feature_dic = {
    "Small": ["A", "C", "D", "G", "N", "P", "S", "T", "V"], "Tiny": ["A", "C", "G", "S", "T"],
    "Aliphatic": ["A", "G", "I", "L", "V"], "Aliphatic_s": ["I", "L", "V"],
    "side_chain_aliphatic_polar": ["C", "M", "S", "T"], "Aromatic": ["F", "W", "Y", "H"], "Aromatic_s": ["F", "W", "Y"],
    "Aromatic_NP": ["F", "W"], "Sulfuric": ["C", "M"], "Hydroxylic": ["S", "T", "Y"], "Amidic": ["N", "Q"],
    "Acidic_side_chain": ["D", "N", "E", "Q"], "Basic_amino_acid": ["H", "K", "R"],
    "Hydrophobic": ["A", "C", "I", "L", "M", "F", "P", "W", "Y", "V"],
    "Hydrophobic_NP": ["A", "G", "I", "L", "M", "F", "P", "W", "V"],
    "Hydrophobic_side_chain": ["A", "I", "L", "M", "F", "W", "Y", "V"],
    "Hydrophobic-Alkyl": ["A", "G", "I", "L", "M", "P", "V"],
    "Hydrophobic-aromatic": ["F", "W"],
    "Hydrophilic": ["E", "D", "H", "K", "N", "Q", "R", "S", "T"],
    "Hydrophilic_polar": ["N", "C", "Q", "S", "T", "Y", "E", "D", "R", "H", "K"],
    "Hydrophylic_side_chain_polar": ["N", "Q", "S", "T", "Y", "E", "D", "R", "H", "K"],
    "Hydrophilic_neutral": ["N", "C", "Q", "S", "T", "Y"],
    "Hydrophilic_side_chain_uncharged": ["N", "Q", "S", "T", "Y"],
    "Hydrophilic_charged": ["E", "D", "R", "H", "K"],
    "Hydrophilic_Acidic_negative_charged": ["D", "E"],
    "Hydrophilic_Basic_positive_charged": ["R", "H", "K"],
    "Hydrophilic_positively_charged": ["R", "K"],
    "Neutral": ["A", "C", "F", "G", "I", "L", "M", "N", "P", "Q", "S", "T", "V", "W"],
    "Neutral_s": ["A", "C", "N", "Q", "S", "T", "Y"],
    "Charged": ["R", "H", "K", "D", "E"],
    "Positively_charged": ["R", "H", "K"],
    "Positively_charged_s": ["R", "K"],
    "Negatively_charged": ["D", "E"],
    "Non_polar_1": ["G", "A", "V", "L", "I", "M", "P", "F", "W"],
    "Non_polar_2": ["A", "I", "L", "M", "P", "V", "F", "W"],
    "Non_polar_1s": ["G", "A", "V", "L", "I", "M"],
    "Non_polar_alkyl": ["G", "A", "V", "L", "I", "M", "P"],
    "Non_polar_aromatic": ["F", "W"],
    "Polar": ["Y", "S", "T", "C", "Q", "N", "E", "D", "K", "H", "R"],
    "Polar_uncharged1": ["G", "S", "T", "C", "Y", "N", "Q"],
    "Polar_uncharged2": ["S", "T", "Q", "N", "C", "P"],
    "Polar_uncharged3": ["Y", "S", "T", "C", "Q", "N"],
    "Polar_uncharged4": ["S", "T", "N", "Q"],
    "Polar_charged": ["E", "D", "R", "H", "K"],
    "Polar_positively_charged": ["R", "H", "K"],
    "Polar_positively_charged_s": ["R", "K"],
    "Polar_negatively_charged": ["D", "E"],
    "Low_complexity": ["S", "P", "G", "R", "K", "Y"],
    "Disorder_promoting": ["A", "R", "G", "Q", "S", "E", "K", "P"],
    "Disorder_promoting_s": ["S", "P", "G", "R"],
    "Order_promoting": ["W", "Y", "F", "I", "L", "V", "C", "N"],
    "Thiolation": ["K", "Q", "E"],
    "EPRS": ["P", "E"],
    "PEVK": ["P", "E", "V", "K"]
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
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")

    seq = ""
    ft_chooser = "F" * int(prop) + "N" * int(100 - prop)
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

        aa_dic = generate_dic(mod.da, amino_acid_list)
        aa_chosen = get_cur_val(aa_dic, random.random())
        codon_list = amino_acid2codon[aa_chosen].split(",")
        codon_dic = generate_dic(mod.dc, codon_list)
        codon_chosen = get_cur_val(codon_dic, random.random())
        seq += codon_chosen

    fseq = ""
    i = 0
    while i < len(seq):
        fseq += seq[i:i + 70] + "\r"
        i += 70

    ap = round(float(seq.count("A")) / len(seq), 2)
    cp = round(float(seq.count("C")) / len(seq), 2)
    gp = round(float(seq.count("G")) / len(seq), 2)
    tp = round(float(seq.count("T")) / len(seq), 2)
    ftp = round(float(ft_count) / (len(seq)/3), 2)
    return fseq, ap, cp, gp, tp, ftp


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
        for i in range(1, number_seq+1):
            length = random.randint(size_int[0], size_int[1])
            seq, ap, cp, gp, tp, ftp = sequence_generator(length, prop, feature, ctrl)
            header = header_generator(length, ap, cp, gp, tp, ftp, feature, i)
            outfile.write(header + "\n" + seq + "\n")
            list_ft_freq.append(ftp)
    mean = 0
    for val in list_ft_freq:
        mean += val
    mean = float(mean)/len(list_ft_freq)
    print("frequence of " + feature + " amino acids in the file : " + str(mean))


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

    Small: A, C, D, G, N, P, S, T, V, Tiny: A, C, G, S, T
    Aliphatic: A, G, I, L, V, Aliphatic_s: I, L, V
    side_chain_aliphatic_polar: C, M, S, T, Aromatic: F, W, Y, H, Aromatic_s: F, W, Y
    Aromatic_NP: F, W, Sulfuric: C, M, Hydroxylic: S, T, Y, Amidic: N, Q
    Acidic_side_chain: D, N, E, Q, Basic_amino_acid: H, K, R
    Hydrophobic: A, C, I, L, M, F, P, W, Y, V
    Hydrophobic_NP: A, G, I, L, M, F, P, W, V
    Hydrophobic_side_chain: A, I, L, M, F, W, Y, V
    Hydrophobic-Alkyl: A, G, I, L, M, P, V
    Hydrophobic-aromatic: F, W
    Hydrophilic: E, D, H, K, N, Q, R, S, T
    Hydrophilic_polar: N, C, Q, S, T, Y, E, D, R, H, K
    Hydrophylic_side_chain_polar: N, Q, S, T, Y, E, D, R, H, K
    Hydrophilic_neutral: N, C, Q, S, T, Y
    Hydrophilic_side_chain_uncharged: N, Q, S, T, Y
    Hydrophilic_charged: E, D, R, H, K
    Hydrophilic_Acidic_negative_charged: D, E
    Hydrophilic_Basic_positive_charged: R, H, K
    Hydrophilic_positively_charged: R, K
    Neutral: A, C, F, G, I, L, M, N, P, Q, S, T, V, W
    Neutral_s: A, C, N, Q, S, T, Y
    Charged: R, H, K, D, E
    Positively_charged: R, H, K
    Positively_charged_s: R, K
    Negatively_charged: D, E
    Non_polar_1: G, A, V, L, I, M, P, F, W
    Non_polar_2: A, I, L, M, P, V, F, W
    Non_polar_1s: G, A, V, L, I, M
    Non_polar_alkyl: G, A, V, L, I, M, P
    Non_polar_aromatic: F, W
    Polar: Y, S, T, C, Q, N, E, D, K, H, R
    Polar_uncharged1: G, S, T, C, Y, N, Q
    Polar_uncharged2: S, T, Q, N, C, P
    Polar_uncharged3: Y, S, T, C, Q, N
    Polar_uncharged4: S, T, N, Q
    Polar_charged: E, D, R, H, K
    Polar_positively_charged: R, H, K
    Polar_positively_charged_s: R, K
    Polar_negatively_charged: D, E
    Low_complexity: S, P, G, R, K, Y
    Disorder_promoting: A, R, G, Q, S, E, K, P
    Disorder_promoting_s: S, P, G, R
    Order_promoting: W, Y, F, I, L, V, C, N
    Thiolation: K, Q, E
    EPRS: P, E
    PEVK: P, E, V, K
    """,
                                     usage='%(prog)s --input input_file.txt [--output an output folder] ')
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
    parser.add_argument('--ctrl', dest='ctrl', help="the control used for fasta generation",
                        default="CCE")

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument('--feature', dest='feature', help="the feature tha will help for sequence generation",
                               required=True)
    required_args.add_argument('--prob', dest='prop', help="the prob to have the feature at each aa position (0-100)",
                               required=True)

    args = parser.parse_args()  # parsing arguments

    if args.filename == "result":
        args.filename = args.feature + "_" + args.prop

    try:
        args.prop = int(args.prop)
        if 0 < args.prop > 100:
            print("ERROR : wrong probability value")
            exit(1)
    except ValueError:
        print("ERROR : wrong probability value")
        exit(1)

    if args.feature not in feature_dic.keys():
        print("ERROR : unknown feature")
        exit(1)

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

    fasta_generator(size_int, args.prop, args.feature, args.nbr_seq, args.output, args.filename, args.ctrl)

if __name__ == "__main__":
    launcher()

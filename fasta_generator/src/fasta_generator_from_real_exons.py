import argparse
import random
import sys
import os
from dicitonary import *
from sets import Set
import numpy as np


feature_dic = {
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
}





def read_CCE_sequence(ctrl):
    """
    Read a file named ctrl + "_reading_frame.csv" in the folder of this script
    :param ctrl: (string) it can be CCE or ACE
    :return: a list of control sequence
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    list_seq = []
    with open(file_dir + "/" + ctrl + "_reading_frame.csv", "r") as seq_file:
        line = seq_file.readline()
        while line:
            list_seq.append(line.split("\t")[1].replace("\n", ""))
            line = seq_file.readline()
    return list_seq


def translator(seq):
    """

    :param seq: (list of string) list of codon
    :return: list of codon translated : i.e. amino acid list
    """
    aa = ""
    for codon in seq:
        aa += codon2aminoAcid[codon]
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
    return float(count) / len(seq)


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
    :return: a dictionary that only contains the unit in unit list not all the unit in ctrl_dic
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


def get_indices_of_feature(codon_list, feature):
    """

    :param codon_list: (list of string) list of codons
    :param feature: (string) the feature of interest
    :return: (list of int) the list of indices in codon_list where there are a codon encoding for
    the feature 'feature'
    """
    aa_seq = translator(codon_list)
    indices = []
    for i in range(len(aa_seq)):
        if aa_seq[i] in feature_dic[feature]:
            indices.append(i)
    return indices


def get_indices(codon_list, aa_list):
    """
    :param aa_list: (list of string) a list of amino acid
    :param codon_list: (list of string) list of codons
    :return: (list of int) the list of indices in codon_list where there are a codon encoding for
    the feature 'feature'
    """
    aa_seq = translator(codon_list)
    indices = []
    for i in range(len(aa_seq)):
        if aa_seq[i] in aa_list:
            indices.append(i)
    return indices


def second_feature_enrichment(seq, prop, feature, ctrl):
    """
    Calculation of the second enrichment/impoverishment of codons encoding amino acids from certain protein feature
    :param seq: (string) the sequence to modify
    :param prop: (list of float) list of proportion to reach for the features in "feature" list
    :param feature: (list of string) list of feature to enriched or impoverished
    :param ctrl: (string) CCE/RD
    :return: the sequence enriched/impoverished for features[1]
    """
    if ctrl != "RD":
        file_dir = os.path.dirname(os.path.realpath(__file__))
        sys.path.insert(0, file_dir + "/control_dic/")
        mod = __import__(ctrl + "_dic")

    my_seq = []
    for i in range(0, len(seq), 3):
        if i+3 <= len(seq):
            my_seq.append(seq[i:i+3])
    # enrichment of a nucleotide

    all_aa = list("ACDEFGHIKLMNPQRSTVWY")
    cur_prop = feature_frequency_calculator(my_seq, feature[1])
    inter = list(Set(feature_dic[feature[0]]).intersection(feature_dic[feature[1]]))
    if len(inter) == len(feature_dic[feature[0]]) or len(inter) == len(feature_dic[feature[1]]):
        inclusion = True
    else:
        inclusion = False
    if cur_prop > prop[1]:
        reg = "-"
        if inclusion:
            if len(feature_dic[feature[0]]) > len(feature_dic[feature[1]]):
                aa_modif = feature_dic[feature[1]]
                aa_to_choose = list(Set(feature_dic[feature[0]]).difference(feature_dic[feature[1]]))
            else:
                aa_modif = list(Set(feature_dic[feature[1]]).difference(feature_dic[feature[0]]))
                aa_to_choose = list(Set(all_aa).difference(Set(feature_dic[feature[1]])))
        else:
            aa_modif = list(Set(feature_dic[feature[1]]).difference(Set(inter)))
            union = Set(feature_dic[feature[1]]).union(feature_dic[feature[0]])
            aa_to_choose = list(Set(all_aa).difference(union))
    else:
        reg = "+"
        if inclusion:
            if len(feature_dic[feature[0]]) > len(feature_dic[feature[1]]):
                aa_modif = list(Set(feature_dic[feature[0]]).difference(feature_dic[feature[1]]))
                aa_to_choose = list(Set(feature_dic[feature[1]]))
            else:
                aa_modif = feature_dic[feature[1]]
                aa_to_choose = list(Set(feature_dic[feature[1]]).difference(feature_dic[feature[0]]))
        else:
            aa_modif = list(Set(all_aa).difference(Set(feature_dic[feature[0]]).union(feature_dic[feature[1]])))
            aa_to_choose = list(Set(feature_dic[feature[1]]).difference(inter))
    temp_reg = reg

    while temp_reg == reg:
        indice_list = get_indices(my_seq, aa_modif)
        if len(indice_list) == 0:
            return "".join(my_seq)
        val = indice_list[random.randint(0, len(indice_list) - 1)]
        aa_dic = generate_dic(mod.da, aa_to_choose)
        aa_chosen = get_cur_val(aa_dic, random.random())
        codon_list = amino_acid2codon[aa_chosen].split(",")
        codon_dic = generate_dic(mod.dc, codon_list)
        codon_chosen = get_cur_val(codon_dic, random.random())
        my_seq[val] = codon_chosen

        cur_prop = feature_frequency_calculator(my_seq, feature[1])
        if abs(prop[1] - cur_prop) > 0.02:
            if cur_prop > prop[1]:
                reg = "-"
            else:
                reg = "+"
        else:
            reg = "ok"

    rseq = "".join(my_seq)

    return rseq


def exon_sequence_generator_with_2_feature(size_int, list_seq, ctrl, feature_interest, prop_feature):
    """
    Generation of fasta sequences having the a feature frequency near of  prop_feature.
    The sequence will be first generated from existing CCE exons sequences
    ctrl variable.
    Those sequence can be enriched in one di-nucleotide if
    dnt_interest is not none.

    :param size_int: list of 2 int) the min size possible and the max size possible of the sequences we want to create
    :param list_seq: (list of string)  a list of control sequence
    :param ctrl: (string) CCE or ACE.
    :param feature_interest: (list of string) the name of 2 features of interest
    :param prop_feature: (list of float) list of proportion to reach for the features in "feature" list
    :return:
        - fseq : the sequence enriched in the 2 features given in  feature_interest
        - ft1 : proportion in feature_interest[0] in fseq
        - ft2, :  proportion in feature_interest[1] in fseq
        - len(rseq) : len of the sequence

    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")
    seq = ""
    while len(seq) < size_int[0] or len(seq) > size_int[1]:
        seq = list_seq[random.randint(0, len(list_seq)-1)]

    my_seq = []
    for i in range(0, len(seq), 3):
        if i+3 <= len(seq):
            my_seq.append(seq[i:i+3])
    # enrichment of a nucleotide

    # print(my_seq)
    cur_prop = feature_frequency_calculator(my_seq, feature_interest[0])
    if cur_prop > prop_feature[0]:
        reg = "-"
    else:
        reg = "+"
    temp_reg = reg
    # count = 0
    while temp_reg == reg:
        if reg == "+":
            val = random.randint(0, len(my_seq) - 1)
            if codon2aminoAcid[my_seq[val]] not in feature_dic[feature_interest[0]]:
                aa_dic = generate_dic(mod.da, feature_dic[feature_interest[0]])
                aa_chosen = get_cur_val(aa_dic, random.random())
                codon_list = amino_acid2codon[aa_chosen].split(",")
                codon_dic = generate_dic(mod.dc, codon_list)
                codon_chosen = get_cur_val(codon_dic, random.random())
                my_seq[val] = codon_chosen
        else:
            indice_list = get_indices_of_feature(my_seq, feature_interest[0])
            val = indice_list[random.randint(0, len(indice_list)-1)]
            aa_list = "ACDEFGHIKLMNPQRSTVWY"
            list_aa = []
            for aa in aa_list:
                if aa not in feature_dic[feature_interest[0]]:
                    list_aa.append(aa)
            aa_dic = generate_dic(mod.da, list_aa)
            aa_chosen = get_cur_val(aa_dic, random.random())
            codon_list = amino_acid2codon[aa_chosen].split(",")
            codon_dic = generate_dic(mod.dc, codon_list)
            codon_chosen = get_cur_val(codon_dic, random.random())
            my_seq[val] = codon_chosen

        cur_prop = feature_frequency_calculator(my_seq, feature_interest[0])
        if abs(prop_feature[0] - cur_prop) > 0.02:
            if cur_prop > prop_feature[0]:
                reg = "-"
            else:
                reg = "+"
        else:
            reg = "ok"

    rseq = "".join(my_seq)
    rseq = second_feature_enrichment(rseq, prop_feature, feature_interest, ctrl)

    fseq = ""
    i = 0
    while i < len(rseq):
        fseq += rseq[i:i + 70] + "\r"
        i += 70

    my_seq = []
    for i in range(0, len(rseq), 3):
        if i+3 <= len(rseq):
            my_seq.append(rseq[i:i+3])
    ft1 = feature_frequency_calculator(my_seq, feature_interest[0])
    ft2 = feature_frequency_calculator(my_seq, feature_interest[1])
    return fseq, ft1, ft2, len(rseq)


def header_generator_2_ft(length, ft1, ft2, feature_interest, num_seq):
    """
    :param length: (int) the length of the current sequence
    :param feature_interest: (list of string) the name of 2 features of interest
    :param ft1: (float) proportion in feature_interest[0] in the current sequence
    :param ft2: (float) proportion in feature_interest[1] in the current sequence
    :param num_seq: (int) the number of sequence we will create
    :return: (string) the header of a sequence
    """
    header = ">seq_" + str(num_seq) + " | " + str(feature_interest[0]) + " : " + \
             str(ft1) + " - " + str(feature_interest[1]) + \
             str(ft2) + " - A : " + \
             "| length : " + str(length)
    return header


def rounder(num, digits, up=True):
    import math
    mul = 10**digits
    if up:
        return math.ceil(num * mul)/mul
    else:
        return math.floor(num*mul)/mul


def exon_sequence_generator(size_int, list_seq, ctrl, feature_interest, prop_feature):
    """
    Generation of fasta sequences having the a feature frequency near of  prop_feature.
    The sequence will be first generated from existing CCE exons sequences
    ctrl variable.
    Those sequence can be enriched in one di-nucleotide if
    dnt_interest is not none.

    :param size_int: list of 2 int) the min size possible and the max size possible of the sequences we want to create
    :param list_seq: (list of string)  a list of control sequence
    :param ctrl: (string) CCE or ACE.
    :param feature_interest: (string) the name of a feature of interest
    :param prop_feature: (float) the wanted proportion of exons in the generated sequence
    :return: fseq, ap, cp, gp, tp, cur_prop, len(rseq)
        -fseq: (string) the nucleotide sequence having a proportion of prop_feature of the feature 'feature_interest'
        - ap : (float) Adenine proportion
        - cp : (float) Cytosine proportion
        - gp : (float) Guanine proportion
        - tp : (float) Thymine proportion
        - cur_prop : (float) feature_interest proportion
        - len(rseq) length of the generated sequence (without the \r)

    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")
    seq = ""
    while len(seq) < size_int[0] or len(seq) > size_int[1]:
        seq = list_seq[random.randint(0, len(list_seq)-1)]

    my_seq = []
    for i in range(0, len(seq), 3):
        if i+3 <= len(seq):
            my_seq.append(seq[i:i+3])
    # enrichment of a nucleotide

    # print(my_seq)
    cur_prop = feature_frequency_calculator(my_seq, feature_interest)
    lenseq = float(len(my_seq))
    if cur_prop > prop_feature:
        freq_interest = rounder(((prop_feature * 2 + (6. / lenseq )) / 2) + float(np.random.randn() / 30),
                                4, True)
    else:
        freq_interest = rounder(((prop_feature * 2 - (6. / lenseq)) / 2) + float(np.random.randn() / 30) ,
                                4, False)
    if freq_interest > 1.:
        freq_interest = 1
    if freq_interest < 0:
        freq_interest = 0
    if cur_prop > freq_interest:
        reg = "-"
    else:
        reg = "+"
    temp_reg = reg
    while temp_reg == reg:
        print(str(cur_prop) + " - " + str(feature_interest) + " - " + str(reg) + " - " + str(temp_reg))
        if reg == "+":
            val = random.randint(0, len(my_seq) - 1)
            if codon2aminoAcid[my_seq[val]] not in feature_dic[feature_interest]:
                aa_dic = generate_dic(mod.da, feature_dic[feature_interest])
                aa_chosen = get_cur_val(aa_dic, random.random())
                codon_list = amino_acid2codon[aa_chosen].split(",")
                codon_dic = generate_dic(mod.dc, codon_list)
                codon_chosen = get_cur_val(codon_dic, random.random())
                my_seq[val] = codon_chosen
        else:
            indice_list = get_indices_of_feature(my_seq, feature_interest)
            val = indice_list[random.randint(0, len(indice_list)-1)]
            aa_list = "ACDEFGHIKLMNPQRSTVWY"
            list_aa = []
            for aa in aa_list:
                if aa not in feature_dic[feature_interest]:
                    list_aa.append(aa)
            aa_dic = generate_dic(mod.da, list_aa)
            aa_chosen = get_cur_val(aa_dic, random.random())
            codon_list = amino_acid2codon[aa_chosen].split(",")
            codon_dic = generate_dic(mod.dc, codon_list)
            codon_chosen = get_cur_val(codon_dic, random.random())
            my_seq[val] = codon_chosen

        cur_prop = feature_frequency_calculator(my_seq, feature_interest)
        if cur_prop > prop_feature:
            temp_reg = "-"
        else:
            temp_reg = "+"
        if cur_prop == prop_feature:
            temp_reg = "ok"

    rseq = "".join(my_seq)
    fseq = ""
    i = 0
    while i < len(rseq):
        fseq += rseq[i:i + 70] + "\n"
        i += 70

    ap = round(float(rseq.count("A")) / len(rseq), 2)
    cp = round(float(rseq.count("C")) / len(rseq), 2)
    gp = round(float(rseq.count("G")) / len(rseq), 2)
    tp = round(float(rseq.count("T")) / len(rseq), 2)
    cur_prop = feature_frequency_calculator(my_seq, feature_interest)

    return fseq, ap, cp, gp, tp, cur_prop, len(rseq)


def header_generator(length, a_prop, c_prop, g_prop, t_prop, ft_prop, ft_name, num_seq):
    """
    :param length: (int) the length of the sequence we will create
    :param a_prop: (float) proportion of A in the sequence we will create
    :param t_prop: (float) proportion of t in the sequence we will create
    :param c_prop: (float) proportion of c in the sequence we will create
    :param g_prop: (float) proportion of g in the sequence we will create
    :param ft_prop: (float) the proportion of the feature ft_name in the sequence we will create
    :param ft_name: (string) the name of the feature of interest whose the frequency must be control
    in the generated sequence
    :param num_seq: (int) the number of sequence we will create
    :return: (string) the header of a sequence
    """
    header = ">seq_" + str(num_seq) + " | " + str(ft_name) + " : " + str(ft_prop) + " - A : " + \
             str(a_prop) + " - C : " + str(c_prop) + " - G : " + str(g_prop)
    header += " - T : " + str(t_prop) + " | length : " + str(length)
    return header


def fasta_generator(size_int, number_seq, output, out_name, feature, feature_prop, ctrl):
    """
    Write a fasta file of number_seq sequences having a size in  'size_int' and proportion corresponding to a_prop,
    t_prop, c_prop, g_prop
    :param size_int: (list of 2 int) the min size possible and the max size possible of the sequences we want to create
    :param number_seq: (int) the number of sequence we will create
    :param output: (string) path where the fasta file will be created
    :param out_name: (string) the name of the fasta file to create
    :param feature: (string) the name of a feature of interest
    :param ctrl: (string) CCE or ACE.
    :param feature_prop: (float) the wanted proportion of exons in the generated sequence
    False if the proportion have to be equal or very close  to what the user specified
    """
    list_seq = read_CCE_sequence(ctrl)
    with open(output + out_name + ".fasta", "w") as outfile:
        list_ft_freq = []
        list_features_prop = [[], []]
        if isinstance(feature, str):
            for i in range(1, number_seq+1):
                fseq, ap, cp, gp, tp, cur_prop, length = exon_sequence_generator(size_int,
                                                                                 list_seq, ctrl, feature, feature_prop)
                header = header_generator(length, ap, cp, gp, tp, cur_prop, feature, i)
                outfile.write(header + "\n" + fseq + "\n")
                list_ft_freq.append(cur_prop)
        else:
            for i in range(1, number_seq+1):
                fseq, ft1, ft2, length = exon_sequence_generator_with_2_feature(size_int,
                                                                                list_seq, ctrl, feature, feature_prop)
                header = header_generator_2_ft(length, ft1, ft2, feature, i)
                outfile.write(header + "\n" + fseq + "\n")
                list_features_prop[0].append(ft1)
                list_features_prop[1].append(ft2)
    if isinstance(feature, str):
        mean = 0
        for val in list_ft_freq:
            mean += val
        mean = float(mean)/len(list_ft_freq)
        print("frequence of " + feature + " amino acids in the file : " + str(mean))
    else:
        for i in range(len(list_features_prop)):
            print("frequence of " + feature[i] + "amino acids in the file : " +
                  str(float(sum(list_features_prop[i])) / number_seq))


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
    Serine : S
    Theronine: T
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
        if 0 > args.prop or args.prop > 1:
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

    if "," in args.feature:
        args.feature = args.feature.split(",")

    if args.output[-1] != "/":
        args.output += "/"

    fasta_generator(size_int, args.nbr_seq, args.output, args.filename, args.feature, args.prop, args.ctrl)

if __name__ == "__main__":
    launcher()

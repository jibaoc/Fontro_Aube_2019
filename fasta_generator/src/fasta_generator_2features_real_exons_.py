#!/usr/bin/python3

import argparse
import random
import sys
import os
from dicitonary import *
import numpy as np
from sets import Set


def debug_msg(msg):
    """
    Print teh enrichment info when debug is set to 1

    :param msg: (string) the message you want to print
    """
    if debug =="1":
        print(msg)


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
        freq_interest = rounder(((prop_feature * 2 + (1. / lenseq )) / 2) + float(np.random.normal(scale=1. / 30)),
                                4, True)
    else:
        freq_interest = rounder(((prop_feature * 2 - (1. / lenseq)) / 2) + float(np.random.normal(scale=1. / 30)),
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
        debug_msg(str(cur_prop) + " - " + str(feature_interest) + " - " + str(reg) + " - " + str(temp_reg))
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
        if cur_prop > freq_interest:
            temp_reg = "-"
        else:
            temp_reg = "+"
        if cur_prop == freq_interest:
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


def redefine_frequence_to_reach(my_seq, feature, interest_frequency):
    """
    This will redefine the frequency to reach for a feature based on the feature frequency and \
    on the given user frequency. And say if the sequence ``my_seq`` need to be enriched or \
    impoverished for this feature

    :param my_seq: (list of string) list of codon
    :param feature: (string) the feature of interest (i.e. to enriched or impoverished)
    :param interest_frequency: (float) the proportion to reach for the feature
    :return: (float)
        - freq_interest (float) : the frequency we want to reach
         - reg (string) "+" if we need to enriched the sequence, "-" if we want to impoverished the sequence
    """
    lenseq = float(len(my_seq))
    cur_prop = feature_frequency_calculator(my_seq, feature)
    if cur_prop > interest_frequency:
        freq_interest = rounder(((interest_frequency * 2 + (1. / lenseq)) / 2) + float(np.random.normal(scale=1. / 100)),
                                4, True)
    else:
        freq_interest = rounder(((interest_frequency * 2 - (1. / lenseq)) / 2) + float(np.random.normal(scale=1. / 100)),
                                4, False)
    if freq_interest > 1.:
        freq_interest = 1
    if freq_interest < 0:
        freq_interest = 0

    if cur_prop > freq_interest:
        reg = "-"
    else:
        reg = "+"

    return freq_interest, reg



def exon_sequence_generator_2feature(size_int, list_seq, ctrl, feature1, feature2, prop1, prop2):

    # Selection of sequences exons
    seq = ""
    while len(seq) < size_int[0] or len(seq) > size_int[1]:
        seq = list_seq[random.randint(0, len(list_seq)-1)]

    # turning the sequence in a list of codon
    my_seq = []
    for i in range(0, len(seq), 3):
        if i+3 <= len(seq):
            my_seq.append(seq[i:i+3])




    freq_interest1, reg1 = redefine_frequence_to_reach(my_seq, feature1, prop1)
    freq_interest2, reg2 = redefine_frequence_to_reach(my_seq, feature2, prop2)

    while freq_interest1 + freq_interest2 > 1:
        debug_msg("redefining freq_interest1 and freq_interest2")
        debug_msg(str(freq_interest1) + " - " + str(freq_interest2))
        debug_msg("done !")
        freq_interest1, reg1 = redefine_frequence_to_reach(my_seq, feature1, prop1)
        freq_interest2, reg2 = redefine_frequence_to_reach(my_seq, feature2, prop2)

    if reg1 != reg2:
        debug_msg("going into dependant enrichment")
        # enrichissement et appauvrissement commun
        my_seq = enrichment_n_impoverishment_2features(my_seq, feature1, feature2, freq_interest1, freq_interest2, reg1, reg2, ctrl)
    else:
        debug_msg("going into independant enrichment")
        my_seq = independant_enrichement_impoverishment_2feature(my_seq, feature1, feature2, freq_interest1, freq_interest2,
                                                        reg1, ctrl)



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

    cur_prop1 = feature_frequency_calculator(my_seq, feature1)
    cur_prop2 = feature_frequency_calculator(my_seq, feature2)

    return fseq, ap, cp, gp, tp, cur_prop1, cur_prop2, len(rseq)


def independant_enrichement_impoverishment_2feature(my_seq, feature1, feature2, ifreq1, ifreq2, reg, ctrl):
    """
    Function that will perform an independent impoverishment or enrichment of 2 features

    :param my_seq: (list of string) list of codons
    :param feature1: (string) the 1st feature of interest
    :param feature2: (string) the 2nd feature of interest
    :param ifreq1: (float) the frequency of interest for the first feature
    :param ifreq2: (float) the frequency of interest for the second feature
    :param reg: (string) "+" if we want to enrich or impoverish my_seq for the ``feature1`` and ``feature2``
    :param ctrl: (string) the control exons bias that we will use
    :return: (string) my_seq mutated having now the wanted frequency of feature1 and 2
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")
    features = [feature1, feature2]
    props = [ifreq1, ifreq2]

    cur_prop1 = feature_frequency_calculator(my_seq, feature1)
    cur_prop2 = feature_frequency_calculator(my_seq, feature2)

    debug_msg("frequency of " + feature1 + " to reach : " + str(ifreq1) + "  - observed : " + str(cur_prop1))
    debug_msg("frequency of " + feature2 + " to reach : " + str(ifreq2) + "  - observed : " + str(cur_prop2))

    for i in range(len(features)):
        cur_prop = feature_frequency_calculator(my_seq, features[i])
        temp_reg = reg
        while temp_reg == reg:
            debug_msg(str(cur_prop) + " - " + str(features[i]) + " - " + str(reg) + " - " + str(temp_reg))
            if reg == "+":
                val = random.randint(0, len(my_seq) - 1)
                j = 0
                while codon2aminoAcid[my_seq[val]] in feature_dic[features[0]] + feature_dic[features[1]]:
                    val = random.randint(0, len(my_seq) - 1)
                    j += 1
                    if j == len(my_seq) * 10:
                        print("WARNING : the enrichment was impossible to complete")
                        return my_seq
                aa_dic = generate_dic(mod.da, feature_dic[features[i]])
                aa_chosen = get_cur_val(aa_dic, random.random())
                codon_list = amino_acid2codon[aa_chosen].split(",")
                codon_dic = generate_dic(mod.dc, codon_list)
                codon_chosen = get_cur_val(codon_dic, random.random())
                my_seq[val] = codon_chosen
            else:
                indice_list = get_indices_of_feature(my_seq, features[i])
                val = indice_list[random.randint(0, len(indice_list)-1)]
                aa_list = "ACDEFGHIKLMNPQRSTVWY"
                list_aa = []
                for aa in aa_list:
                    if aa not in feature_dic[features[0]] + feature_dic[features[1]]:
                        list_aa.append(aa)
                if len(list_aa) == 0:
                    print("Warning : nothing to perform enrichment analysis")
                    return my_seq
                aa_dic = generate_dic(mod.da, list_aa)
                aa_chosen = get_cur_val(aa_dic, random.random())
                codon_list = amino_acid2codon[aa_chosen].split(",")
                codon_dic = generate_dic(mod.dc, codon_list)
                codon_chosen = get_cur_val(codon_dic, random.random())
                my_seq[val] = codon_chosen

            cur_prop = feature_frequency_calculator(my_seq, features[i])
            if cur_prop > props[i]:
                temp_reg = "-"
            else:
                temp_reg = "+"
            if cur_prop == props[i]:
                temp_reg = "ok"
    return my_seq


def enrichment_n_impoverishment_2features(my_seq, feature1, feature2, ifreq1, ifreq2, reg1, reg2, ctrl):
    """
    Function that perform an enrichment of a feature ``feature1`` while it impoverished the other feature ``feature2``

    :param my_seq: (list of string) list of codons
    :param feature1: (string) the 1st feature of interest
    :param feature2: (string) the 2nd feature of interest
    :param ifreq1: (float) the frequency of interest for the first feature
    :param ifreq2: (float) the frequency of interest for the second feature
    :param reg1: (string) "+" if we want to enrich or impoverish my_seq for the ``feature1``
    :param reg2: (string) "+" if we want to enrich or impoverish my_seq for the ``feature2``
    :param ctrl: (string) the control exons bias that we will use
    :return: (string) my_seq mutated having now the wanted frequency of feature1 and 2
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")

    if reg1 == "+":
        features = [feature1, feature2]
        props = [ifreq1, ifreq2]
        regs = [reg1, reg2]
    else:
        features = [feature2, feature1]
        props = [ifreq2, ifreq1]
        regs = [reg2, reg1]
    cur_props = [feature_frequency_calculator(my_seq, features[0]), feature_frequency_calculator(my_seq, features[1])]
    temp_regs = [regs[0], regs[1]]

    debug_msg("frequency of " + features[0] + " to reach : " + str(props[0]) + "  - observed : " + str(cur_props[0]))
    debug_msg("frequency of " + features[1] + " to reach : " + str(props[1]) + "  - observed : " + str(cur_props[1]))


    while temp_regs[0] == regs[0] and temp_regs[1] == regs[1]:
        debug_msg(str(cur_props[0]) + " - " + str(features[0]) + " - " + str(props[0]))
        debug_msg(str(cur_props[1]) + " - " + str(features[1]) + " - " + str(props[1]))
        val = random.randint(0, len(my_seq) - 1)
        while codon2aminoAcid[my_seq[val]] not in feature_dic[features[1]]:
            val = random.randint(0, len(my_seq) - 1)
        aa_dic = generate_dic(mod.da, feature_dic[features[0]])
        aa_chosen = get_cur_val(aa_dic, random.random())
        codon_list = amino_acid2codon[aa_chosen].split(",")
        codon_dic = generate_dic(mod.dc, codon_list)
        codon_chosen = get_cur_val(codon_dic, random.random())
        my_seq[val] = codon_chosen

        cur_props = [feature_frequency_calculator(my_seq, features[0]),
                     feature_frequency_calculator(my_seq, features[1])]
        if cur_props[0] < props[0]:
            temp_regs[0] = "+"
        else:
            temp_regs[0] = "ok"
        if cur_props[1] > props[1]:
            temp_regs[1] = "-"
        else:
            temp_regs[1] = "ok"

    if temp_regs[0] == "ok" and temp_regs[1] == "ok":
        return my_seq
    elif temp_regs[0] == "+" and temp_regs[1] == "ok":
        # enrichment of features[0] without touching features[1]
        debug_msg("enrichment_finisher")
        my_seq = enrichment_finisher(my_seq, features[0], props[0], features[1], ctrl)
        return my_seq
    elif temp_regs[0] == "ok" and temp_regs[1] == "-":
        # impoverishment of features[1] without touching features[0]
        debug_msg("impoverishment finisher")
        my_seq = impoverishment_finisher(my_seq, features[1], props[1], features[0], ctrl)
        return my_seq
    else:
        print("temp_regs[0] = " + str(temp_regs[0]) + " - temp_regs[1] = " + str(temp_regs[1]))
        print("Unexpected values !")
        print("Terminating...")
        exit(1)


def enrichment_finisher(my_seq, efeature, eprop, ifeature, ctrl):
    """
    This function terminates the enrichment of ``efeature`` in ``my_seq`` without touching the feature \
    ``ifeature``.

    :param my_seq: (list of string) list of codons
    :param efeature: (string) the feature for which we want to finish the enrichment
    :param eprop: (float) the proportion of ``efeature`` that we want to reach in ``my_seq``
    :param ifeature: (string) the feature we don't want to touch
    :param ctrl: (string) the control codons bias we want to use
    :return:(string) my_seq with the wanted enrichment
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")

    cur_prop = feature_frequency_calculator(my_seq, efeature)

    reg = "+"
    temp_reg = "+"

    while reg == temp_reg:
        debug_msg(str(cur_prop) + " - " + str(efeature) + " - " + str(eprop))
        val = random.randint(0, len(my_seq) - 1)
        j = 0
        while codon2aminoAcid[my_seq[val]] in feature_dic[efeature] + feature_dic[ifeature]:
            val = random.randint(0, len(my_seq) - 1)
            j += 1
            if j == len(my_seq) * 10:
                print("WARNING : the enrichment was impossible to complete")
                return my_seq
        aa_dic = generate_dic(mod.da, feature_dic[efeature])
        aa_chosen = get_cur_val(aa_dic, random.random())
        codon_list = amino_acid2codon[aa_chosen].split(",")
        codon_dic = generate_dic(mod.dc, codon_list)
        codon_chosen = get_cur_val(codon_dic, random.random())
        my_seq[val] = codon_chosen

        cur_prop = feature_frequency_calculator(my_seq, efeature)
        if cur_prop < eprop:
            temp_reg = "+"
        else:
            temp_reg = "-"
    return my_seq


def impoverishment_finisher(my_seq, ifeature, iprop, efeature, ctrl):
    """
    This function terminates the enrichment of ``efeature`` in ``my_seq`` without touching the feature \
    ``ifeature``.

    :param my_seq: (list of string) list of codons
    :param ifeature: (string) the feature for which we want to finish the impoverishment
    :param iprop: (float) the proportion of ``ifeature`` that we want to reach in ``my_seq``
    :param efeature: (string) the feature we don't want to touch
    :param ctrl: (string) the control codons bias we want to use
    :return:(string) my_seq with the wanted impoverishment
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")

    cur_prop = feature_frequency_calculator(my_seq, ifeature)

    reg = "-"
    temp_reg = "-"

    while reg == temp_reg:
        debug_msg(str(cur_prop) + " - " + str(ifeature) + " - " + str(iprop))
        indice_list = get_indices_of_feature(my_seq, ifeature)
        if len(indice_list) < 1:
            print("Warning, nothing to impoverish")
            return my_seq
        val = indice_list[random.randint(0, len(indice_list) - 1)]
        aa_list = "ACDEFGHIKLMNPQRSTVWY"
        list_aa = []
        for aa in aa_list:
            if aa not in feature_dic[ifeature] + feature_dic[efeature]:
                list_aa.append(aa)
        if len(list_aa) == 0:
            print("Warning : nothing to perform impoverishment finisher")
            return my_seq
        aa_dic = generate_dic(mod.da, list_aa)
        aa_chosen = get_cur_val(aa_dic, random.random())
        codon_list = amino_acid2codon[aa_chosen].split(",")
        codon_dic = generate_dic(mod.dc, codon_list)
        codon_chosen = get_cur_val(codon_dic, random.random())
        my_seq[val] = codon_chosen

        cur_prop = feature_frequency_calculator(my_seq, ifeature)
        if cur_prop > iprop:
            temp_reg = "-"
        else:
            temp_reg = "+"

    return my_seq






def header_generator(length, a_prop, c_prop, g_prop, t_prop, ft_prop1, ftprop2, ft_name1, ft_name2, num_seq):
    """
    :param length: (int) the length of the sequence we will create
    :param a_prop: (float) proportion of A in the sequence we will create
    :param t_prop: (float) proportion of t in the sequence we will create
    :param c_prop: (float) proportion of c in the sequence we will create
    :param g_prop: (float) proportion of g in the sequence we will create
    :param ft_prop1: (float) the proportion of the feature ft_name1 in the sequence we will create
    :param ft_prop1: (float) the proportion of the feature ft_name2 in the sequence we will create
    :param ft_name1: (string) the name of the feature1 of interest whose the frequency must be control
    in the generated sequence
    :param ft_name2: (string) the name of the feature2 of interest whose the frequency must be control
    in the generated sequence
    :param num_seq: (int) the number of sequence we will create
    :return: (string) the header of a sequence
    """
    header = ">seq_" + str(num_seq) + " | " + str(ft_name1) + " : " + str(ft_prop1) + \
             " -" + str(ft_name1) + " : " + str(ft_prop1) + " - A : " + \
             str(a_prop) + " - C : " + str(c_prop) + " - G : " + str(g_prop)
    header += " - T : " + str(t_prop) + " | length : " + str(length)
    return header


def fasta_generator(size_int, number_seq, output, out_name, ft1, ft2, prop1, prop2, ctrl):
    """
    Write a fasta file of number_seq sequences having a size in  'size_int' and proportion corresponding to a_prop,
    t_prop, c_prop, g_prop
    :param size_int: (list of 2 int) the min size possible and the max size possible of the sequences we want to create
    :param number_seq: (int) the number of sequence we will create
    :param output: (string) path where the fasta file will be created
    :param out_name: (string) the name of the fasta file to create
    :param ft1: (string) the name of a ft1 of interest
    :param ctrl: (string) CCE or ACE.
    :param prop1: (float) the wanted proportion of exons in the generated sequence
    False if the proportion have to be equal or very close  to what the user specified
    """
    list_seq = read_CCE_sequence(ctrl)
    with open(output + out_name + ".fasta", "w") as outfile:
        list_features_prop = [[], []]
        for i in range(1, number_seq+1):
            fseq, ap, cp, gp, tp, cur_prop1, cur_prop2, length = exon_sequence_generator_2feature(size_int, list_seq, ctrl,
                                                                                      ft1, ft2, prop1, prop2)
            header = header_generator(length, ap, cp, gp, tp, cur_prop1, cur_prop2, ft1, ft2, i)
            outfile.write(header + "\n" + fseq + "\n")
            list_features_prop[0].append(cur_prop1)
            list_features_prop[1].append(cur_prop2)

    features = [ft1, ft2]
    for i in range(len(features)):
        mean = 0
        for val in list_features_prop[i]:
            mean += val
        mean = float(mean)/len(list_features_prop[i])
        print("frequence of " + features[i] + " amino acids in the file : " + str(mean))


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

    Small#1: [A, C, D, G, N, P, S, T, V],
    Small#2: [A, C, D, G, N, P, S, T],
    Large : [F, I, K, L, M, R, W, Y],
    Disorder_promoting#1: [A, E, G, K, P, Q, R, S],
    Order_promoting#1: [C, F, I, L, N, W, V, Y],
    Disorder_promoting#2: [A, E, G, K, P, Q, S],
    Order_promoting#2: [C, F, H, I, L, M, N, W, V, Y],
    Polar_uncharged#1: [C, N, Q, S, T, Y],
    Polar_uncharged#2: [N, Q, S, T, Y],
    Charged: [R, H, K, D, E],
    Hydrophilic#1: [D, E, K, N, Q, R],
    Hydrophobic#1: [A, C, F, I, L, M, V],
    Hydrophilic#2: [D, E, H, K, N, Q, R, S, T],
    Hydrophobic#2: [A, C, F, I, L, M, P, V, W, Y],
    Hydroxylic: [S, T, Y],
    Negatively_charged: [D, E],
    Positively_charged: [R, H, K],
    """)
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
    parser.add_argument('--debug', dest="debug", help="1 if we want to debug 0 else", default="0")

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument('--feature1', dest='feature1', help="the feature that will be enriched in the sequence generation",
                               required=True)
    required_args.add_argument('--feature2', dest='feature2', help="the feature that will be impoverished in the sequence generation",
                               required=True)
    required_args.add_argument('--prop1', dest='prop1', help="the prop in the sequence in the feature1",
                               required=True)
    required_args.add_argument('--prop2', dest='prop2', help="the prop in the sequence in the feature 2",
                               required=True)

    args = parser.parse_args()  # parsing arguments

    if args.filename == "result":
        args.filename = args.ctrl + "_" + args.feature1 + "," + args.feature2  + "_" +  args.prop1 + "," + args.prop2

    try:
        args.prop1 = float(args.prop1)
        args.prop2 = float(args.prop2)
        if 0 > args.prop1 or args.prop1 > 1:
            print("ERROR : wrong probability value for prop 1")
            exit(1)
        if 0 > args.prop1 or args.prop1 > 1:
            print("ERROR : wrong probability value for prop2")
            exit(1)
    except ValueError:
        print("ERROR : wrong probability value")
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

    inter =  Set(feature_dic[args.feature1]).intersection(Set(feature_dic[args.feature2]))
    if len(inter) > 1:
        print("This program doesn't handle feature groups overlap")
        print("you must choose opposed features !")
        print("Exiting...")
        exit(1)

    global debug
    debug = args.debug



    fasta_generator(size_int, args.nbr_seq, args.output, args.filename, args.feature1, args.feature2,
                    args.prop1, args.prop2, args.ctrl)

if __name__ == "__main__":
    launcher()

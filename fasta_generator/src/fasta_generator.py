"""Description: script that can generate a fasta file with random sequences."""

################################
#          IMPORTS             #
################################
import random
import argparse
import os
from math import floor
import copy
from operator import itemgetter
import sys

################################
#      Global variables        #
################################

list_name = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT",
             "TA", "TC", "TT", "TG"]
# Link each IUPAC letter to their corresponding codons
iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
         'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}

################################
#         Functions            #
################################

#################################################################
# generation of fasta sequence based on nucleotide proportion
#################################################################


def sequence_generator(length, a_prop, t_prop, c_prop, g_prop):
    """
    :param length: (int) the length of the sequence we will create
    :param a_prop: (float) proportion of A in the sequence we will create
    :param t_prop: (float) proportion of t in the sequence we will create
    :param c_prop: (float) proportion of c in the sequence we will create
    :param g_prop: (float) proportion of g in the sequence we will create
    :return: (string) a sequence of length 'length' and with the proportion of nucleotides very near of the values
    given by the variables a_prop, c_prop, g_prop, t_prop
    """
    seq = "A" * int(floor(length * a_prop)) + "C" * int(floor(length * c_prop)) + "G" * int(floor(length * g_prop)) + \
          "T" * int(floor(length * t_prop))

    if len(seq) != length:
        new_seq = "A" * int(floor(100 * a_prop)) + "C" * int(floor(100 * c_prop)) + "G" * int(floor(100 * g_prop)) + \
          "T" * int(floor(100 * t_prop))
        while len(seq) != length:
            seq += new_seq[random.randint(0, len(new_seq)-1)]

    rseq = "".join(random.sample(seq, len(seq)))
    fseq = ""
    i = 0
    while i < len(rseq):
        fseq += rseq[i:i+70] + "\r"
        i += 70

    ap = round(float(rseq.count("A")) / len(rseq), 2)
    cp = round(float(rseq.count("C")) / len(rseq), 2)
    gp = round(float(rseq.count("G")) / len(rseq), 2)
    tp = round(float(rseq.count("T")) / len(rseq), 2)
    return fseq, ap, cp, gp, tp


def flexible_sequence_generator(length, a_prop, t_prop, c_prop, g_prop):
    """
    :param length: (int) the length of the sequence we will create
    :param a_prop: (float) proportion of A in the sequence we will create
    :param t_prop: (float) proportion of t in the sequence we will create
    :param c_prop: (float) proportion of c in the sequence we will create
    :param g_prop: (float) proportion of g in the sequence we will create
    :return: (string) a sequence of length 'length' and with the proportion of nucleotides quite near of the values
    given by the variables a_prop, c_prop, g_prop, t_prop
    """
    seq = ""
    new_seq = "A" * int(round(100 * a_prop)) + "C" * int(round(100 * c_prop)) + "G" * int(round(100 * g_prop)) + \
        "T" * int(round(100 * t_prop))
    while len(seq) != length:

        seq += new_seq[random.randint(0, len(new_seq)-1)]
    fseq = ""
    i = 0
    while i < len(seq):
        fseq += seq[i:i+70] + "\r"
        i += 70

    ap = round(float(seq.count("A")) / len(seq), 2)
    cp = round(float(seq.count("C")) / len(seq), 2)
    gp = round(float(seq.count("G")) / len(seq), 2)
    tp = round(float(seq.count("T")) / len(seq), 2)
    return fseq, ap, cp, gp, tp


def header_generator(length, a_prop, c_prop, g_prop, t_prop, num_seq):
    """
    :param length: (int) the length of the sequence we will create
    :param a_prop: (float) proportion of A in the sequence we will create
    :param t_prop: (float) proportion of t in the sequence we will create
    :param c_prop: (float) proportion of c in the sequence we will create
    :param g_prop: (float) proportion of g in the sequence we will create
    :param num_seq: (int) the number of sequence we will create
    :return: (string) the header of a sequence
    """
    header = ">seq" + str(num_seq) + " | A : " + str(a_prop) + " - C : " + str(c_prop) + " - G : " + str(g_prop)
    header += " - T : " + str(t_prop) + " | length : " + str(length)
    return header


def fasta_generator(size_int,  a_prop, t_prop, c_prop, g_prop, number_seq, output, out_name, flexible):
    """
    Write a fasta file of number_seq sequences having a size in  'size_int' and proportion corresponding to a_prop,
    t_prop, c_prop, g_prop
    :param size_int: (list of 2 int) the min size possible and the max size possible of the sequences we want to create
    :param a_prop: (float) proportion of A in the sequence we will create
    :param t_prop: (float) proportion of t in the sequence we will create
    :param c_prop: (float) proportion of c in the sequence we will create
    :param g_prop: (float) proportion of g in the sequence we will create
    :param number_seq: (int) the number of sequence we will create
    :param output: (string) path where the fasta file will be created
    :param out_name: (string) the name of the fasta file to create
    :param flexible: (boolean) True if the proportion of the sequence should only be close to what the user specified
    False if the proportion have to be equal or very close  to what the user specified
    """
    with open(output + out_name + ".fasta", "w") as outfile:
        for i in range(1, number_seq+1):
            length = random.randint(size_int[0], size_int[1])
            if flexible:
                seq, ap, cp, gp, tp = flexible_sequence_generator(length, a_prop, t_prop, c_prop, g_prop)
            else:
                seq, ap, cp, gp, tp = sequence_generator(length, a_prop, t_prop, c_prop, g_prop)
            header = header_generator(length, ap, cp, gp, tp, i)
            outfile.write(header + "\n" + seq + "\n")

#################################################################
# generation of fasta sequence based on di-nucleotide proportion
#################################################################


def next_dnt(my_dnt, last_nt):
    """
    Find the next dnt of a sequence being generated.

    Find the next dnt according to the last nucleotide of the sequence
    being generated. If we can choose a di-nucleotide beginning by this last
    nucleotide (because it is in my_dnt) we will only add the second nucleotide of the di-nucleotide chosen.
    If we cant we will add the entire di-nucleotide
    :param my_dnt: (dictionary of list of string of 2 letters), a dictionary that contains all the
    di-nucleotide that must be in the sequence being generated. The keys of the dictionary are A, T, G, C.
    Each key is associated with a list of di-nucleotides beginning by this key.
    :param last_nt: (character) the last nucleotide of the sequence being generated.
    :return: 1 or 2 nucleotides:
     1 nucleotide if we can choose a di-nucleotide beginning by last_nt
     2 nucleotides else.
    """
    len_dic = {"A": len(my_dnt["A"]), "C": len(my_dnt["C"]), "G": len(my_dnt["G"]), "T": len(my_dnt["T"])}
    sorted_ld = sorted(len_dic.items(), key=itemgetter(1), reverse=True)
    if len_dic[last_nt] > 0:
        # if we can find a di-nucleotide beginning by last_nt ...
        next_nt = my_dnt[last_nt][random.randint(0, len(my_dnt[last_nt])-1)]
        del (my_dnt[next_nt[0]][my_dnt[next_nt[0]].index(next_nt)])
        # ... we return only the last letter of the di-nucleotide chosen
        return next_nt[1]
    elif len_dic[last_nt] == 0:
        # if we can't find a di-nucleotide beginning by last_nt ...
        if sorted_ld[0][1] == 0:
            return "Done"
        else:
            list_of_remaning_nt = []
            for i in range(len(sorted_ld)):
                if sorted_ld[0][1] > 0:
                    list_of_remaning_nt.append(sorted_ld[0][0])
            last_nt = list_of_remaning_nt[random.randint(0, len(list_of_remaning_nt)-1)]
            next_nt = my_dnt[last_nt][random.randint(0, len(my_dnt[last_nt]) - 1)]
            del (my_dnt[next_nt[0]][my_dnt[next_nt[0]].index(next_nt)])
            # ... we return the all di-nucleotide chosen
            return next_nt


def dinucleotide_calculator(seq):
    """
    Calculate the di-nucleotide frequency in the sequence.
    Only unambiguous di-nucleotides.

    :param seq: (string) a nucleotide sequence
    :return: (dictionary) a dictionary containing the frequency of every possible di-nucleotides
    """
    dic = {"AA": 0., "AT": 0., "AG": 0., "AC": 0., "TA": 0., "TT": 0., "TG": 0., "TC": 0.,
           "GA": 0., "GT": 0., "GG": 0., "GC": 0., "CA": 0., "CT": 0., "CG": 0., "CC": 0.}

    cur = {"AA": 0., "AT": 0., "AG": 0., "AC": 0., "TA": 0., "TT": 0., "TG": 0., "TC": 0.,
           "GA": 0., "GT": 0., "GG": 0., "GC": 0., "CA": 0., "CT": 0., "CG": 0., "CC": 0.}
    if len(seq) > 1:
        for j in range(len(seq) - 1):
            cur[seq[j:j + 2]] += 1
        for key in dic.keys():
            dic[key] += float(cur[key]) / (len(seq) - 1)

    return dic


def flexible_dnt_sequence_generator(length, dnt_list):
    """
    Generation of completely random sequences.
    Generation based on di-nucleotide proportions

    :param length: (int) the approximate length of the
    sequence to generate.
    :param dnt_list: (list of float) the proportion of
    each dnt in the sequence - the proportions here
    corresponds to the di-nucleotides in global variable
    list_name
    :return: my_seq, dnt_prop_txt, dnt_prop
     - my_seq : (string) the random sequence generated
     - dnt_prop_txt : (string) the proportion of each
     di-nucleotides in my_seq
     - dnt_prop (dictionary of float) proportion (values)
     of each di-nucleotide (keys)
    """
    my_dnt_list = []
    for i in range(len(dnt_list)):
        my_dnt_list += [list_name[i]] * int(round(dnt_list[i] * 400))

    my_dnt = {"A": [], "C": [], "G": [], "T": []}

    # select only the wanted dnt
    for i in range(length-2):  # because of the init dnt
        dnt = my_dnt_list[random.randint(0, len(my_dnt_list)-1)]
        my_dnt[dnt[0]].append(dnt)

    my_seq = my_dnt_list[random.randint(0, len(my_dnt_list)-1)]

    # building the sequence
    res = ""
    while res != "Done":
        res = next_dnt(my_dnt, my_seq[-1])
        if res != "Done":
            my_seq += res

    # calculation of di-nucleotide proportion in my_seq
    dnt_prop_txt = ""
    dnt_prop = dinucleotide_calculator(my_seq)
    for key in dnt_prop.keys():
        dnt_prop_txt += str(key) + ": " + str(dnt_prop[key]) + " | "

    dnt_prop_txt = dnt_prop_txt[0:len(dnt_prop_txt)-3]

    return my_seq, dnt_prop_txt, dnt_prop


def header_dnt_generator(length, header_text, num_seq):
    """
    Generation of an header for the fasta sequence.

    :param length: (int) the length of the sequence
    :param header_text: (string) the text of the header : the dnt frequencies
    :param num_seq: (int) the number of the sequence generated.
    :return:
    """
    header = ">seq" + str(num_seq) + " | length : " + str(length) + " | " + header_text
    return header


def fasta_dnt_generator(size_int, dnt_list, number_seq, output, out_name):
    """
    Generate a fasta file containing random sequences
    generated with a list of dnt frequencies.

    :param size_int:(list of 2 float) first float : min length possible
    of the sequences in the fasta file. second float : max length possible of
    the sequences in the fasta file.
    :param dnt_list: (list of float)
    :param number_seq: (int) the number of sequences we want to generate.
    :param output: (string) the folder where the file will be created
    :param out_name: (string) the name of the fasta file.
    """
    res_stat = [0 for i in range(16)]
    with open(output + out_name + ".fasta", "w") as outfile:
        for i in range(1, number_seq+1):
            length = random.randint(size_int[0], size_int[1])
            seq, text_header, dnt_prop = flexible_dnt_sequence_generator(length, dnt_list)
            for j in range(len(list_name)):
                res_stat[j] += dnt_prop[list_name[j]]
            header = header_dnt_generator(len(seq), text_header, i)
            outfile.write(header + "\n" + seq + "\n")
    for j in range(len(res_stat)):
        res_stat[j] /= number_seq
    return res_stat

#################################################################
# generation of fasta sequence based on di-nucleotide proportion
# and codon frequency of ACCE/CCE/ALL fasterDB exons
#################################################################


def dinucleotide_calculator_bis(seq):
    """
    Calculate the di-nucleotide frequency in the sequence.
    For di-nucleotides expanded to iupac code.

    :param seq: (string) a nucleotide sequence
    :return: (dictionary) a dictionary containing the frequency of every possible di-nucleotides
    """
    dic = {"AA": 0., "AT": 0., "AG": 0., "AC": 0., "TA": 0., "TT": 0., "TG": 0., "TC": 0.,
           "GA": 0., "GT": 0., "GG": 0., "GC": 0., "CA": 0., "CT": 0., "CG": 0., "CC": 0.}

    cur = {"AA": 0., "AT": 0., "AG": 0., "AC": 0., "TA": 0., "TT": 0., "TG": 0., "TC": 0.,
           "GA": 0., "GT": 0., "GG": 0., "GC": 0., "CA": 0., "CT": 0., "CG": 0., "CC": 0.}
    if len(seq) > 1:
        for j in range(len(seq) - 1):
            cur[seq[j:j + 2]] += 1
        for key in dic.keys():
            dic[key] += float(cur[key]) / (len(seq) - 1)
    # Calculation of ambiguous di-nucleotide frequency
    for nt1 in iupac.keys():
        for nt2 in iupac.keys():
            dic[nt1 + nt2] = 0.
            for letter1 in iupac[nt1]:
                for letter2 in iupac[nt2]:
                    dic[nt1 + nt2] += dic[letter1 + letter2]
    return dic


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


def get_cur_codon(ctr_dic, value):
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


def exon_sequence_generator(length, ctrl, dnt_interest):
    """
    Generation of fasta sequences having the same codon frequency as
    the one in CCE/ACE/ALL exons in fasterDB according to the
    ctrl variable.
    Those sequence can be enriched in one di-nucleotide if
    dnt_interest is not none.

    :param length: (int) the length of the sequence to generate
    :param ctrl: (string) CCE or ACE or ALL.
    :param dnt_interest: (tuple of a string and a float) the first value is the dnt,
    the other is its proportion.
    :return: my_seq, dnt_prop_txt, dnt_prop
     - my_seq : (string) the random sequence generated
     - dnt_prop_txt : (string) the proportion of each
     di-nucleotides in my_seq
     - dnt_prop (dictionary of float) proportion (values)
     of each di-nucleotide (keys)
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")
    seq = ""
    ctrl_dic = ctrl_dic_adapter(mod.dc)

    # generation of the sequence
    for i in range(length/3):
        codon = get_cur_codon(ctrl_dic, random.random())
        if codon is None:
            print("Something went wrong ! ")
            exit(1)
        seq += str(codon)

    # enrichment of a di-nucleotide
    if dnt_interest is not None:
        seq = list(seq)
        dnt_prop = dinucleotide_calculator_bis("".join(seq))

        if dnt_prop[dnt_interest[0]] > dnt_interest[1]:
            reg = "-"
        else:
            reg = "+"
        temp_reg = reg
        while temp_reg == reg:
            if reg == "+":
                val = random.randint(0, len(seq)-2)
                if dnt_interest[0][0] in ["A", "T", "G", "C"]:
                    seq[val] = dnt_interest[0][0]
                else:
                    seq[val] = iupac[dnt_interest[0][0]][random.randint(0, len(iupac[dnt_interest[0][0]])-1)]
                if dnt_interest[0][1] in ["A", "T", "G", "C"]:
                    seq[val+1] = dnt_interest[0][1]
                else:
                    seq[val+1] = iupac[dnt_interest[0][1]][random.randint(0, len(iupac[dnt_interest[0][1]])-1)]
            else:
                break
            dnt_prop = dinucleotide_calculator_bis("".join(seq))
            if dnt_prop[dnt_interest[0]] >= dnt_interest[1]:
                temp_reg = "-"
            else:
                temp_reg = "+"
        seq = "".join(seq)
    dnt_prop_txt = ""
    dnt_prop = dinucleotide_calculator_bis(seq)
    for key in dnt_prop.keys():
        dnt_prop_txt += str(key) + ": " + str(dnt_prop[key]) + " | "

    dnt_prop_txt = dnt_prop_txt[0:len(dnt_prop_txt) - 3]

    return seq, dnt_prop_txt, dnt_prop


def nt_freq_calculator(seq, nt):
    """
    :param seq: (string) a nucleotide sequence
    :param nt: (string of one character) nt nt for which we want to calculate the frequency
    :return: the freq of nt in seq
    """
    if nt in ["A", "T", "G", "C"]:
        nt_prop = float(seq.count(nt)) / len(seq)
    else:
        count = 0
        for n in iupac[nt]:
            count += seq.count(n)
        nt_prop = float(count) / len(seq)
    return nt_prop


def exon_nt_sequence_generator(length, ctrl, nt_interest):
    """
    Generation of fasta sequences having the same codon frequency as
    the one in CCE/ACE/ALL exons in fasterDB according to the
    ctrl variable.
    Those sequence can be enriched in one di-nucleotide if
    dnt_interest is not none.

    :param length: (int) the length of the sequence to generate
    :param ctrl: (string) CCE or ACE or ALL.
    :param nt_interest: (tuple of a string and a float) the first value is the nt,
    the other is its proportion.
    :return: my_seq, dnt_prop_txt, nt_prop
     - my_seq : (string) the random sequence generated
     - nt_prop_txt : (string) the proportion of each
     nucleotides in my_seq
     - nt_prop (list of float) proportion
     of each nucleotide
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")
    seq = ""
    ctrl_dic = ctrl_dic_adapter(mod.dc)

    # generation of the sequence
    for i in range(length/3):
        codon = get_cur_codon(ctrl_dic, random.random())
        if codon is None:
            print("Something went wrong ! ")
            exit(1)
        seq += str(codon)

    # enrichment of a nucleotide
    if nt_interest is not None:
        seq = list(seq)
        nt_prop = nt_freq_calculator("".join(seq), nt_interest[0])

        if nt_prop > nt_interest[1]:
            reg = "-"
        else:
            reg = "+"
        temp_reg = reg
        while temp_reg == reg:
            if reg == "+":
                val = random.randint(0, len(seq)-1)
                if nt_interest[0] in ["A", "T", "G", "C"]:
                    seq[val] = nt_interest[0]
                else:
                    seq[val] = iupac[nt_interest[0]][random.randint(0, len(iupac[nt_interest[0]])-1)]
            else:
                break
            nt_prop = nt_freq_calculator("".join(seq), nt_interest[0])
            if nt_prop >= nt_interest[1]:
                temp_reg = "-"
            else:
                temp_reg = "+"
        seq = "".join(seq)
    nt_prop_txt = "A : " + str(float(seq.count("A")) / len(seq)) + " | C :" + str(float(seq.count("C")) / len(seq)) + " | "
    nt_prop_txt += "G : " + str(float(seq.count("G")) / len(seq)) + " | T :" + str(float(seq.count("G")) / len(seq))

    nt_prop = []
    nt_prop.append(float(seq.count("A")) / len(seq))
    nt_prop.append(float(seq.count("C")) / len(seq))
    nt_prop.append(float(seq.count("G")) / len(seq))
    nt_prop.append(float(seq.count("T")) / len(seq))

    return seq, nt_prop_txt, nt_prop


def ctrl_fasta_dnt_generator(size_int, dnt_interest, number_seq, output, out_name, ctrl):
    """

    :param size_int: (list of 2 float) first float : min length possible
    :param dnt_interest: (tuple of a string and a float) the first value is the dnt,
    the other is its proportion.
    :param number_seq: (int) the number of sequences we want to generate.
    :param output: (string) the folder where the file will be created
    :param out_name: (string) the name of the fasta file.
    :param ctrl: (string) CCE or ACE or ALL.
    """
    freq_nt = [0, 0, 0, 0]
    res_stat = [0 for i in range(16)]
    with open(output + out_name + ".fasta", "w") as outfile:
        for i in range(1, number_seq+1):
            length = random.randint(size_int[0], size_int[1])
            if len(dnt_interest[0]) > 1:
                seq, text_header, dnt_prop = exon_sequence_generator(length, ctrl, dnt_interest)
                for j in range(len(list_name)):
                    res_stat[j] += dnt_prop[list_name[j]]
            else:
                seq, text_header, nt_prop = exon_nt_sequence_generator(length, ctrl, dnt_interest)
                for i in range(len(freq_nt)):
                    freq_nt[i] += nt_prop[i]
            header = header_dnt_generator(len(seq), text_header, i)
            outfile.write(header + "\n" + seq + "\n")
    if len(dnt_interest[0]) > 1:
        for j in range(len(res_stat)):
            res_stat[j] /= number_seq
        return res_stat
    else:
        for i in range(len(freq_nt)):
            freq_nt[i] /=  number_seq
        seq = ""
        nt_list = ["A", "C", "G", "T"]
        for i in range(len(nt_list)):
            seq += str(nt_list[i]) + " : " + str(freq_nt[i]) + " - "
        return seq


######################################################
#             Manager functions
######################################################


def my_format(list_prop):
    """
    Turning string number into float and none value to none value if possible.
    If it's not possible, stop the program

    :param list_prop: (list of string of none) the list of a_prop, t_prop, c_prop, g_prop : even if those variable
    corresponds to a number, there are string when the are processing with this function or None value
    :return: (list of float or None)
    """
    percentage = False
    new_list = []
    tot_prop = 0.
    try:
        for val in list_prop:
            if val is not None:
                val = float(val)
                if val > 1:
                    percentage = True
        for val in list_prop:
            if val is not None:
                val = float(val)
                if val < 0:
                    print("ERROR : propotion value below 0")
                    exit(1)
                if percentage:
                    new_list.append(val / 100)
                    tot_prop += val / 100
                else:
                    new_list.append(val)
                    tot_prop += val
            else:
                new_list.append(None)
        if percentage and tot_prop > 100:
            print("ERROR : the sum of proportion given is greater than 100")
            exit(1)
        elif not percentage and tot_prop > 1:
            print("ERROR : the sum of proportion given is greater than 1")
            exit(1)
    except ValueError:
        print("ERROR : wrong proportion values.")
        print("Exiting...")
        exit(1)
    return new_list


def handling_nt_proportion(list_nt):
    """
    Turning all string value to int, if possible. The None value are estimated thanks to all the values in
    list_nt. The sum of every value in list_nt must equal 1
    """
    nbr_none = 0
    none_list = []
    prop = my_format(list_nt)
    temp = copy.deepcopy(prop)
    for i in range(len(prop)):
        if prop[i] is None:
            nbr_none += 1
            none_list.append(i)
    if nbr_none == 0:
        return prop
    else:
        for i in range(len(none_list)):
            temp[none_list[i]] = 0
        for i in range(len(none_list)):

            prop[none_list[i]] = (1. - sum(temp)) / len(none_list)
        return prop


def test_dnt_nt(nt_tuple, dnt_tuple):
    """
    Function to test if we launch a generation based on nucleotides
    of di-nucleotide frequencies.
    :param nt_tuple: (list of 4 float or None values)
    :param dnt_tuple: (list of 16 float or None values)
    :return: mix if both nt and dnt parameter was filled
    nt if none or one nt value was filled
    dnt if one dnt value was filled
    """
    nt_value = False
    dnt_value = False

    for val in nt_tuple:
        if val is not None:
            nt_value = True
    for val in dnt_tuple:
        if val is not None:
            dnt_value = True

    if nt_value and dnt_value:
        return "mix"

    elif nt_value:
        return "nt"

    elif dnt_value:
        return "dnt"

    else:
        return "nt"


def display_dnt_prop(list_dnt, message):
    """
    Display the di-nucleotide proportion.

    :param list_dnt: (list of float) list of di-nucleotides proportion
    :param message: (float) a string to display
    """
    print(message)
    res = ""
    for i in range(len(list_dnt)):
        res += list_name[i] + " : " + str(list_dnt[i]) + " | "
    print(res)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""Given a number of sequence N and 3 proportions values ( of A, C, G
                                     nucleotides)  create a fasta file of N random sequences with proportions of A, C, G
                                     as specified above.
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
    parser.add_argument("--prop_A", dest='prop_A', help="the proportion of alanine in the fasta file",
                        default=None)
    parser.add_argument('--prop_C', dest='prop_C', help="the proportion of cytosine in the fasta file",
                        default=None)
    parser.add_argument('--prop_G', dest='prop_G', help="the proportion of guanine in the fasta file",
                        default=None)
    parser.add_argument('--prop_T', dest='prop_T', help="the proportion of thymine in the fasta file",
                        default=None)
    parser.add_argument('--flexible', dest='flexible', help="true if you want to allow a litle distortion in your given"
                                                            " proportions, false else",
                        default=False)
    parser.add_argument('--ctrl', dest='ctrl', help="control dic we want to use",
                        default=None)
    parser.add_argument('--nt_dnt', dest='nt_dnt', help="the dnt or the nt you want to enriched in the ctrl (ACE/CCE/ALL) sequences",
                        default=None)
    parser.add_argument('--freq', dest='freq', help="the freq of the dnt you want to enriched in the ctrl "
                                                    "(ACE/CCE/ALL) sequences",
                        default=None)

    parser.add_argument('--AA', dest='AA', help="the proportion of AA in the fasta file",
                        default=None)
    parser.add_argument('--AC', dest='AC', help="the proportion of AC in the fasta file",
                        default=None)
    parser.add_argument('--AG', dest='AG', help="the proportion of AG in the fasta file",
                        default=None)
    parser.add_argument('--AT', dest='AT', help="the proportion of AT in the fasta file",
                        default=None)
    parser.add_argument('--CA', dest='CA', help="the proportion of CA in the fasta file",
                        default=None)
    parser.add_argument('--CC', dest='CC', help="the proportion of CC in the fasta file",
                        default=None)
    parser.add_argument('--CG', dest='CG', help="the proportion of CG in the fasta file",
                        default=None)
    parser.add_argument('--CT', dest='CT', help="the proportion of CT in the fasta file",
                        default=None)
    parser.add_argument('--GA', dest='GA', help="the proportion of GA in the fasta file",
                        default=None)
    parser.add_argument('--GC', dest='GC', help="the proportion of GC in the fasta file",
                        default=None)
    parser.add_argument('--GG', dest='GG', help="the proportion of GG in the fasta file",
                        default=None)
    parser.add_argument('--GT', dest='GT', help="the proportion of GT in the fasta file",
                        default=None)
    parser.add_argument('--TA', dest='TA', help="the proportion of TA in the fasta file",
                        default=None)
    parser.add_argument('--TC', dest='TC', help="the proportion of TC in the fasta file",
                        default=None)
    parser.add_argument('--TG', dest='TG', help="the proportion of TG in the fasta file",
                        default=None)
    parser.add_argument('--TT', dest='TT', help="the proportion of TT in the fasta file",
                        default=None)

    args = parser.parse_args()  # parsing arguments

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

    if args.flexible == "True":
        args.flexible = True
    elif args.flexible == "False" or args.flexible is False:
        args.flexible = False
    else:
        print("WARNING : unregonized boolean value for flexible argument")
        print("Setting it to False ! ")
        args.flexible = False

    nt_tuple = (args.prop_A, args.prop_T, args.prop_C, args.prop_G)
    dnt_tuple = (args.AA, args.AC, args.AG, args.AT, args.CA, args.CC, args.CG, args.CT, args.GA, args.GC, args.GG,
                 args.GT, args.TA, args.TC, args.TT, args.TG)

    res = test_dnt_nt(nt_tuple, dnt_tuple)
    size_int = [args.size_inf, args.size_max]
    if args.ctrl is None:
        if res == "nt":

            args.prop_A, args.prop_T, args.prop_C, args.prop_G = \
                handling_nt_proportion((args.prop_A, args.prop_T, args.prop_C, args.prop_G))

            print("Nucleotides proportion : ")
            print("A : " + str(args.prop_A) + " - C : " + str(args.prop_C) + " - G : " + str(args.prop_G) + " - T : " +
                  str(args.prop_T))

            fasta_generator(size_int,  args.prop_A, args.prop_T, args.prop_C, args.prop_G, args.nbr_seq, args.output,
                            args.filename, args.flexible)

        elif res == "dnt":
            dnt_tuple = handling_nt_proportion(dnt_tuple)
            display_dnt_prop(dnt_tuple, "di-nucleotides proportions")
            res_stat = fasta_dnt_generator(size_int, dnt_tuple, args.nbr_seq, args.output, args.filename)
            display_dnt_prop(res_stat, "proportion in the file : ")

        else:
            print("ouch will be hard")
    elif args.ctrl in ["CCE", "ACE"]:
        if args.nt_dnt is not None and args.freq is not None:
            interest_dnt = [args.nt_dnt, float(args.freq)]
        else:
            interest_dnt = None

        res_stat = ctrl_fasta_dnt_generator(size_int, interest_dnt, args.nbr_seq, args.output, args.filename, args.ctrl)

        if not isinstance(res_stat, str):
            display_dnt_prop(res_stat, "proportion in the file : ")
        else:
            print("proportion in the file : ")
            print(res_stat)
    else:
        print("Unrocognized control...")
        exit(1)


if __name__ == "__main__":
    launcher()

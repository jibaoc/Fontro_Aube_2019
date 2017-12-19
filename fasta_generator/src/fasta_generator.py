
#######
# IMPORTS
#######
import random
import argparse
import os
from math import floor
import copy
from operator import itemgetter
import sys

list_name = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT",
             "TA", "TC", "TT", "TG"]
# Link each IUPAC letter to their corresponding codons
iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
         'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}

def next_dnt(my_dnt, last_nt):

    len_dic = {"A" : len(my_dnt["A"]), "C" : len(my_dnt["C"]), "G": len(my_dnt["G"]), "T" : len(my_dnt["T"])}
    sorted_ld = sorted(len_dic.items(), key=itemgetter(1), reverse=True)
    if len_dic[last_nt] > 0:
        next_nt = None
        ct = 0
        while next_nt is None and ct <4:
            cut_nt = last_nt + sorted_ld[ct][0]
            if cut_nt in my_dnt[cut_nt[0]]:
                next_nt = cut_nt
                del (my_dnt[next_nt[0]][my_dnt[next_nt[0]].index(next_nt)])
                return next_nt
        if next_nt is None:
            print ("Error - the next_nt cannot be found")
    elif len_dic[last_nt] == 0:
        if sorted_ld[0][1] == 0:
            return "Done"
        else:
            cur_dnt = sorted_ld[0][0]
            next_nt = None
            ct = 0
            while next_nt is None and ct < 4:
                cut_nt = last_nt + sorted_ld[ct][0]
                if cut_nt in my_dnt[cut_nt[0]]:
                    next_nt = cut_nt
                    del (my_dnt[next_nt[0]][my_dnt[next_nt[0]].index(next_nt)])
                    return next_nt
            if next_nt is None:
                print ("Error - the next_nt cannot be found")

def next_dnt2(my_dnt, last_nt):
    len_dic = {"A" : len(my_dnt["A"]), "C" : len(my_dnt["C"]), "G": len(my_dnt["G"]), "T" : len(my_dnt["T"])}
    sorted_ld = sorted(len_dic.items(), key=itemgetter(1), reverse=True)
    if len_dic[last_nt] > 0:
        next_nt = my_dnt[last_nt][random.randint(0, len(my_dnt[last_nt])-1)]
        del (my_dnt[next_nt[0]][my_dnt[next_nt[0]].index(next_nt)])
        return next_nt[1]
    elif len_dic[last_nt] == 0:
        if sorted_ld[0][1] == 0:
            return "Done"
        else:
            list_of_remaning_nt = []
            for i in  range(len(sorted_ld)):
                if sorted_ld[0][1] > 0:
                    list_of_remaning_nt.append(sorted_ld[0][0])
            last_nt = list_of_remaning_nt[random.randint(0, len(list_of_remaning_nt)-1)]
            next_nt = my_dnt[last_nt][random.randint(0, len(my_dnt[last_nt]) - 1)]
            del (my_dnt[next_nt[0]][my_dnt[next_nt[0]].index(next_nt)])
            return next_nt


def dinucleotide_calculator(seq):
    """
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

def dinucleotide_calculator_bis(seq):
    """
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

    for nt1 in iupac.keys():
        for nt2 in iupac.keys():
            dic[nt1 + nt2] = 0.
            for letter1 in iupac[nt1]:
                for letter2 in iupac[nt2]:
                    dic[nt1 + nt2] += dic[letter1 + letter2]
    print(dic)
    return dic



def flexible_dnt_sequence_generator(length, dnt_list):
    seq = ""
    my_dnt_list = []
    for i in range(len(dnt_list)):
        my_dnt_list += [list_name[i]] * int(round(dnt_list[i] * 400))

    my_dnt =  {"A" : [], "C" : [], "G": [], "T" : []}

    # select only the wanted dnt
    for i in range(length-2): # because of the init dnt
        dnt = my_dnt_list[random.randint(0, len(my_dnt_list)-1)]
        my_dnt[dnt[0]].append(dnt)

    my_seq = my_dnt_list[random.randint(0, len(my_dnt_list)-1)]



    # building the sequence
    res = ""
    while res != "Done":
        res = next_dnt2(my_dnt, my_seq[-1])
        if res != "Done":
            my_seq += res


    dnt_prop_txt = ""
    dnt_prop = dinucleotide_calculator(my_seq)
    for key in dnt_prop.keys():
        dnt_prop_txt += key + ": " + str(dnt_prop[key]) + " | "

    dnt_prop_txt = dnt_prop_txt[0:len(dnt_prop_txt)-3]

    return my_seq, dnt_prop_txt, dnt_prop


def ctrl_dic_adapter(dic):
    res_dic = {}
    for key in dic.keys():
        if key != "all":
           res_dic[key] = float(dic[key])/dic["all"]
    list_key = res_dic.keys()
    tmp = 0.
    for key in list_key:
        res_dic[key] = res_dic[key] + tmp
        tmp = res_dic[key]
    sorted_res = sorted(res_dic.items(), key=lambda l: l[1], reverse=False)
    interval_dic = {}
    for i in range(len(sorted_res)):
        if i == 0:
            interval_dic[sorted_res[i][0]] = [0, sorted_res[i][1]]
        elif i < len(sorted_res)-1:
            interval_dic[sorted_res[i][0]] = [sorted_res[i-1][1], sorted_res[i][1]]
        else:
            interval_dic[sorted_res[i][0]] = [sorted_res[i-1][1], sorted_res[i][1] + 0.00000001]
    return interval_dic


def get_cur_codon(ctr_dic, value):
    for key in ctr_dic.keys():
        if ctr_dic[key][0] <= value < ctr_dic[key][1]:
            return key
    return None



def exon_sequence_generator(length, ctrl, dnt_interest):
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")
    seq = ""
    ctrl_dic = ctrl_dic_adapter(mod.dc)


    for i in range(length/3):
        codon = get_cur_codon(ctrl_dic, random.random())
        if codon is None:
            print("Something went wrong ! ")
            exit(1)
        seq += str(codon)

    if dnt_interest is not None:
        seq = list(seq)
        dnt_prop = dinucleotide_calculator_bis("".join(seq))

        if dnt_prop[dnt_interest[0]] > dnt_interest[1]:
            reg = "-"
        else:
            reg = "+"
        temp_reg = reg
        print(str(dnt_prop[dnt_interest[0]]) + " - " + str(dnt_interest[1]) + " - " + str(reg))
        g=0
        while temp_reg == reg and g==0:
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
                print("already enriched")
                break
            dnt_prop = dinucleotide_calculator_bis("".join(seq))
            print("".join(seq))
            print(str(dnt_prop[dnt_interest[0]]) + " - " + str(dnt_interest[1]) + " - " + str(reg) + " -" + str(temp_reg))
            if dnt_prop[dnt_interest[0]] >= dnt_interest[1]:
                print("ok")
                temp_reg = "-"
            else:
                print("ko")
                temp_reg = "+"
            g=0
        seq = "".join(seq)
    dnt_prop_txt = ""
    dnt_prop = dinucleotide_calculator_bis(seq)
    for key in dnt_prop.keys():
        dnt_prop_txt += key + ": " + str(dnt_prop[key]) + " | "

    dnt_prop_txt = dnt_prop_txt[0:len(dnt_prop_txt) - 3]

    return seq, dnt_prop_txt, dnt_prop

# with hexa nt control
"""
def dic_adapter(ctrl_dic, key):
    new_dic = {}
    list_key = [key + "A", key + "C", key + "G", key + "T"]
    count = 0.
    for key in list_key:
        val = ctrl_dic[key]
        new_dic[key] = val
        count += val

    new_dic["all"] = count
    res_dic = ctrl_dic_adapter(new_dic)
    return res_dic

def exon_sequence_generator2(length, ctrl, dnt_interest):
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")
    ctrl_dic = ctrl_dic_adapter(mod.d6)

    codon = get_cur_codon(ctrl_dic, random.random())
    seq = codon
    for i in range(length-6):
        cur_dic = dic_adapter(mod.d6, seq[-5:])
        codon = get_cur_codon(cur_dic, random.random())
        if codon is None:
            print("Something went wrong ! ")
            exit(1)
        seq += codon[-1]


    if dnt_interest is not None:
        seq = list(seq)
        dnt_prop = dinucleotide_calculator("".join(seq))

        if dnt_prop[dnt_interest[0]] > dnt_interest[1]:
            reg = "-"
        else:
            reg = "+"
        temp_reg = reg
        print(str(dnt_prop[dnt_interest[0]]) + " - " + str(dnt_interest[1]) + " - " + str(reg))
        g=0
        while temp_reg == reg and g==0:
            if reg == "+":
                val = random.randint(0, len(seq)-2)
                if dnt_interest[0][0] in ["A", "T", "G", "C"]:
                    seq[val] = dnt_interest[0][0]
                else:
                    seq[val] = iupac[dnt_interest[0][0]][random.randint(0, len(iupac[dnt_interest[0][0]]))]
                if dnt_interest[0][1] in ["A", "T", "G", "C"]:
                    seq[val] = dnt_interest[0][1]
                else:
                    seq[val] = iupac[dnt_interest[0][1]][random.randint(0, len(iupac[dnt_interest[0][1]]))]
            else:
                print("exiting")
                break
            dnt_prop = dinucleotide_calculator("".join(seq))
            print(str(dnt_prop[dnt_interest[0]]) + " - " + str(dnt_interest[1]) + " - " + str(reg) + " -" + str(temp_reg))
            if dnt_prop[dnt_interest[0]] >= dnt_interest[1]:
                print("ok")
                temp_reg = "-"
            else:
                print("ko")
                temp_reg = "+"
            g=0
        seq = "".join(seq)
    dnt_prop_txt = ""
    dnt_prop = dinucleotide_calculator(seq)
    for key in dnt_prop.keys():
        dnt_prop_txt += key + ": " + str(dnt_prop[key]) + " | "

    dnt_prop_txt = dnt_prop_txt[0:len(dnt_prop_txt) - 3]

    return seq, dnt_prop_txt, dnt_prop
"""

def header_dnt_generator(length, header_text, num_seq):
    header = ">seq" + str(num_seq) + " | length : " +  str(length) + " | " + header_text
    return header


def fasta_dnt_generator(size_int, dnt_list , number_seq, output, out_name):
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
        res_stat[j] = res_stat[j] / number_seq
    return res_stat

def ctrl_fasta_dnt_generator(size_int, dnt_interest , number_seq, output, out_name, ctrl):
    res_stat = [0 for i in range(16)]
    with open(output + out_name + ".fasta", "w") as outfile:
        for i in range(1, number_seq+1):
            length = random.randint(size_int[0], size_int[1])
            seq, text_header, dnt_prop = exon_sequence_generator(length, ctrl, dnt_interest)
            for j in range(len(list_name)):
                res_stat[j] += dnt_prop[list_name[j]]
            header = header_dnt_generator(len(seq), text_header, i)
            outfile.write(header + "\n" + seq + "\n")
    for j in range(len(res_stat)):
        res_stat[j] = res_stat[j] / number_seq
    return res_stat


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


def my_format(list_prop):
    """
    Turning string number into float and none value to none value if possible. If it's not possible, stop the program
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
    parser.add_argument('--dnt', dest='dnt', help="the dnt you want to enriched in the ctrl (ACE/CCE/ALL) sequences",
                        default=None)
    parser.add_argument('--freq', dest='freq', help="the  freq of the dnt you want to enriched in the ctrl (ACE/CCE/ALL) sequences",
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
    dnt_tuple = (args.AA, args.AC, args.AG, args.AT, args.CA, args.CC, args.CG, args.CT, args.GA, args.GC, args.GG, args.GT,
                 args.TA, args.TC, args.TT, args.TG)

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
    else:
        if args.dnt is not None and args.freq is not None:
            interest_dnt = [args.dnt, float(args.freq)]
        else:
            interest_dnt = None
        res_stat = ctrl_fasta_dnt_generator(size_int, interest_dnt, args.nbr_seq, args.output, args.filename, args.ctrl)
        display_dnt_prop(res_stat, "proportion in the file : ")



if __name__ == "__main__":
    launcher()

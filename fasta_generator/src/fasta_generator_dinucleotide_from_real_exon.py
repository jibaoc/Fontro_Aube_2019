import argparse
import os
import random
import sys
from sets import Set


list_name = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT",
             "TA", "TC", "TT", "TG"]
# Link each IUPAC letter to their corresponding codons
iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
         'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}


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


def get_indices(seq, dnt_interest):
    """
    :param seq: (string) a nucleotide sequence
    :param dnt_interest: (string) a di-nt of interest
    :return: list of in : the position where the dnt_of interest are found in seq
    """
    nt_list = ["A", "T", "G", "C"]
    if dnt_interest[0] in nt_list and dnt_interest[1] in nt_list:
        dnt_of_interest = [dnt_interest]
    elif dnt_interest[0] in nt_list and dnt_interest[1] not in nt_list:
        dnt_of_interest = [dnt_interest[0] + b for b in iupac[dnt_interest[1]]]
    elif dnt_interest[0] not in nt_list and dnt_interest[1] in nt_list:
        dnt_of_interest = [a + dnt_interest[1] for a in iupac[dnt_interest[0]]]
    else:
        dnt_of_interest = [a + b for a in iupac[dnt_interest[0]] for b in iupac[dnt_interest[1]]]
    indices = []
    for i in range(len(seq)-1):
        if "".join(seq[i:i+2]) in dnt_of_interest:
            indices.append(i)
    return indices


def generate_dic(ctrl_dic, unit_list):
    """

    :param ctrl_dic: (dictionary of int) for each unit (dnt or nt) give it's count in the control set of exons
    :param unit_list: (list of string) list of nt or dnt
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


def exon_sequence_generator(size_int, list_seq, ctrl, dnt_interest):
    """
    Generation of fasta sequences having the same codon frequency as
    the one in CCE/ACE/ALL exons in fasterDB according to the
    ctrl variable.
    Those sequence can be enriched in one di-nucleotide if
    dnt_interest is not none.
    :param list_seq: (list of string)  a list of control sequence
    :param size_int: (list of 2 float) first float : min length possible
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

    seq = ""
    while len(seq) < size_int[0] or len(seq) > size_int[1]:
        seq = list_seq[random.randint(0, len(list_seq)-1)]

    # enrichment of a nucleotide
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
                indices = get_indices(seq, dnt_interest[0])
                ind = indices[random.randint(0, len(indices)-1)]
                dnt_choosed = list_name[random.randint(0, len(list_name)-1)]
                seq[ind] = dnt_choosed[0]
                seq[ind+1] = dnt_choosed[1]

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


def get_nt_indices(seq, nt):
    """
    :param seq: (string) a nucleotide sequence
    :param nt: (string of one character) nt nt for which we want to find the indices in seq
    :return: (list of int) list of indices in
    """
    indices = []
    if nt in ["A", "T", "G", "C"]:
        for i in range(len(seq)):
            if seq[i] == nt:
                indices.append(i)
    else:
        for i in range(len(seq)):
            if seq[i] in iupac[nt]:
                indices.append(i)
    return indices

def exon_nt_sequence_generator(size_int, list_seq, ctrl, nt_interest):
    """
    Generation of fasta sequences having the same codon frequency as
    the one in CCE/ACE/ALL exons in fasterDB according to the
    ctrl variable.
    Those sequence can be enriched in one di-nucleotide if
    dnt_interest is not none.

    :param list_seq: (list of string)  a list of control sequence
    :param size_int: (list of 2 float) first float : min length possible
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
    while len(seq) < size_int[0] or len(seq) > size_int[1]:
        seq = list_seq[random.randint(0, len(list_seq)-1)]

    nt_list = []
    if nt_interest[0] in ["A", "T", "G", "C"]:
        for letter in ["A", "T", "G", "C"]:
            if letter != nt_interest[0]:
                nt_list.append(letter)
    else:
        for letter in ["A", "T", "G", "C"]:
            if letter not in iupac[nt_interest[0]]:
                nt_list.append(letter)
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
                indices = get_nt_indices(seq, nt_interest[0])
                ind = indices[random.randint(0, len(indices) - 1)]
                nt_dic = generate_dic(mod.nt, nt_list)
                nt_chosen = get_cur_val(nt_dic, random.random())
                seq[ind] = nt_chosen

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


def second_nt_enrichment(seq, prop_interest, nt_interest, ctrl):


    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__(ctrl + "_dic")


    nt_list = []
    for letter in ["A", "T", "G", "C"]:
        if letter not in nt_interest:
            nt_list.append(letter)
    # enrichment of a nucleotide
    if nt_interest is not None:
        seq = list(seq)
        nt_prop = nt_freq_calculator("".join(seq), nt_interest[1])

        if nt_prop > prop_interest[1]:
            reg = "-"
        else:
            reg = "+"
        temp_reg = reg
        while temp_reg == reg:
            if reg == "+":
                indices_to_keep = get_nt_indices(seq, nt_interest[0])
                all_indices = range(len(seq))
                indice = list(Set(all_indices).difference(indices_to_keep))
                val = indice[random.randint(0, len(indice) - 1)]
                seq[val] = nt_interest[1]
            else:
                indices = get_nt_indices(seq, nt_interest[1])
                ind = indices[random.randint(0, len(indices) - 1)]
                nt_dic = generate_dic(mod.nt, nt_list)
                nt_chosen = get_cur_val(nt_dic, random.random())
                seq[ind] = nt_chosen
            nt_prop = nt_freq_calculator("".join(seq), nt_interest[1])
            if nt_prop >= prop_interest[1]:
                temp_reg = "-"
            else:
                temp_reg = "+"

        seq = "".join(seq)

    return seq


def exon_2_nt_sequence_generator(size_int, list_seq, ctrl, nt_interest, prop_interest):
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

    nt_list = []

    for letter in ["A", "T", "G", "C"]:
        if letter != nt_interest[0]:
            nt_list.append(letter)

    # generation of the sequence
    while len(seq) < size_int[0] or len(seq) > size_int[1]:
        seq = list_seq[random.randint(0, len(list_seq)-1)]
    print(seq)
    # enrichment of a nucleotide
    if nt_interest is not None:
        seq = list(seq)
        nt_prop = nt_freq_calculator("".join(seq), nt_interest[0])

        if nt_prop > prop_interest[0]:
            reg = "-"
        else:
            reg = "+"
        temp_reg = reg
        while temp_reg == reg:
            if reg == "+":
                val = random.randint(0, len(seq)-1)
                seq[val] = nt_interest[0]
            else:
                indices = get_nt_indices(seq, nt_interest[0])
                ind = indices[random.randint(0, len(indices) - 1)]
                nt_dic = generate_dic(mod.nt, nt_list)
                nt_chosen = get_cur_val(nt_dic, random.random())
                seq[ind] = nt_chosen
            nt_prop = nt_freq_calculator("".join(seq), nt_interest[0])
            if nt_prop >= prop_interest[0]:
                temp_reg = "-"
            else:
                temp_reg = "+"
        seq = "".join(seq)

    seq = second_nt_enrichment(seq, prop_interest, nt_interest, ctrl)
    print(seq)
    nt_prop_txt = "A : " + str(float(seq.count("A")) / len(seq)) + " | C :" + str(float(seq.count("C")) / len(seq)) + " | "
    nt_prop_txt += "G : " + str(float(seq.count("G")) / len(seq)) + " | T :" + str(float(seq.count("G")) / len(seq))

    nt_prop = []
    nt_prop.append(float(seq.count("A")) / len(seq))
    nt_prop.append(float(seq.count("C")) / len(seq))
    nt_prop.append(float(seq.count("G")) / len(seq))
    nt_prop.append(float(seq.count("T")) / len(seq))
    return seq, nt_prop_txt, nt_prop


def header_dnt_generator(length, header_text, num_seq):
    """
    Generation of an header for the fasta sequence.

    :param length: (int) the length of the sequence
    :param header_text: (string) the text of the header : the dnt frequencies
    :param num_seq: (int) the number of the sequence generated.
    :return: (string) the header oft the sequence
    """
    header = ">seq" + str(num_seq) + " | length : " + str(length) + " | " + header_text
    return header


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


def ctrl_fasta_nt_dnt_generator(size_int, nt_dnt_interest, number_seq, output, out_name, ctrl):
    """
    :param size_int: (list of 2 float) first float : min length possible
    :param nt_dnt_interest: (tuple of a string and a float) the first value is the dnt or the nt,
    the other is its proportion.
    :param number_seq: (int) the number of sequences we want to generate.
    :param output: (string) the folder where the file will be created
    :param out_name: (string) the name of the fasta file.
    :param ctrl: (string) CCE or ACE or ALL.
    """
    freq_nt = [0, 0, 0, 0]
    list_seq = read_CCE_sequence(ctrl)
    res_stat = [0 for i in range(16)]
    with open(output + out_name + ".fasta", "w") as outfile:
        for i in range(1, number_seq+1):
            if len(nt_dnt_interest[0]) > 1:
                seq, text_header, dnt_prop = exon_sequence_generator(size_int, list_seq, ctrl, nt_dnt_interest)
                for j in range(len(list_name)):
                    res_stat[j] += dnt_prop[list_name[j]]
            else:
                seq, text_header, nt_prop = exon_nt_sequence_generator(size_int, list_seq, ctrl, nt_dnt_interest)
                for i in range(len(freq_nt)):
                    freq_nt[i] += nt_prop[i]
            header = header_dnt_generator(len(seq), text_header, i)
            outfile.write(header + "\n" + seq + "\n")
    if len(nt_dnt_interest[0]) > 1:
        for j in range(len(res_stat)):
            res_stat[j] /= number_seq
        return res_stat
    else:
        for i in range(len(freq_nt)):
            freq_nt[i] /= number_seq
        seq = ""
        nt_list = ["A", "C", "G", "T"]
        for i in range(len(nt_list)):
            seq += str(nt_list[i]) + " : " + str(freq_nt[i]) + " - "
        return seq


def ctrl_fasta_2_nt_generator(size_int, nt_interest, prop_interest, number_seq, output, out_name, ctrl):
    """

    :param size_int: (list of 2 float) first float : min length possible
    :param nt_interest: (tuple of a 2 string) the 2 nt we want to enriched impoverished
    :param prop_interest: (tuple of 2 float) the prop of the 2 nucleotide we want to enriched
    :param number_seq: (int) the number of sequences we want to generate.
    :param output: (string) the folder where the file will be created
    :param out_name: (string) the name of the fasta file.
    :param ctrl: (string) CCE or ACE or ALL.
    """
    freq_nt = [0, 0, 0, 0]
    list_seq = read_CCE_sequence(ctrl)
    with open(output + out_name + ".fasta", "w") as outfile:
        for i in range(1, number_seq + 1):
            seq, text_header, nt_prop = exon_2_nt_sequence_generator(size_int, list_seq, ctrl, nt_interest, prop_interest)
            for i in range(len(freq_nt)):
                freq_nt[i] += nt_prop[i]
            header = header_dnt_generator(len(seq), text_header, i)
            outfile.write(header + "\n" + seq + "\n")

    for i in range(len(freq_nt)):
        freq_nt[i] /= number_seq
    seq = ""
    nt_list = ["A", "C", "G", "T"]
    for i in range(len(nt_list)):
        seq += str(nt_list[i]) + " : " + str(freq_nt[i]) + " - "
    return seq


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
    parser.add_argument('--ctrl', dest='ctrl', help="control dic we want to use",
                        default="CCE")
    parser.add_argument('--nt_dnt', dest='nt_dnt', help="the dnt or nt you want to enriched in the ctrl (ACE/CCE/ALL) sequences",
                        default=None)
    parser.add_argument('--freq', dest='freq', help="the freq of the dnt you want to enriched in the ctrl "
                                                    "(ACE/CCE/ALL) sequences",
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

    size_int = [args.size_inf, args.size_max]
    interest_dnt = None
    if "," in args.nt_dnt:
        args.nt_dnt = args.nt_dnt.split(",")
        args.freq = args.freq.split(",")
        args.freq[0] = float(args.freq[0])
        args.freq[1] = float(args.freq[1])
    else:
        if args.nt_dnt is not None and args.freq is not None:
            interest_dnt = [args.nt_dnt, float(args.freq)]
        else:
            interest_dnt = None

    if not isinstance(args.nt_dnt, list):
        res_stat = ctrl_fasta_nt_dnt_generator(size_int, interest_dnt, args.nbr_seq, args.output, args.filename, args.ctrl)
        if len(interest_dnt[0]) > 1:
            display_dnt_prop(res_stat, "proportion in the file : ")
        else:
            print("proportion in the file : ")
            print(res_stat)
    else:
        if len(args.nt_dnt) > 2:
            print("List of nt too long...")
            exit(1)
        else:
            res_stat = ctrl_fasta_2_nt_generator(size_int, args.nt_dnt, args.freq, args.nbr_seq,
                                                 args.output, args.filename, args.ctrl)
            print("proportion in the file : ")
            print(res_stat)


if __name__ == "__main__":
    launcher()













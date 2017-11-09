
#######
# IMPORTS
#######
import random
import argparse
import os
from math import floor


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

    ap = round(float(rseq.count("A")) / len(rseq), 2)
    cp = round(float(rseq.count("C")) / len(rseq), 2)
    gp = round(float(rseq.count("G")) / len(rseq), 2)
    tp = round(float(rseq.count("T")) / len(rseq), 2)
    return rseq, ap, cp, gp, tp


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

    ap = round(float(seq.count("A")) / len(seq), 2)
    cp = round(float(seq.count("C")) / len(seq), 2)
    gp = round(float(seq.count("G")) / len(seq), 2)
    tp = round(float(seq.count("T")) / len(seq), 2)
    return seq, ap, cp, gp, tp


def header_generator(length, a_prop, t_prop, c_prop, g_prop, num_seq):
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
                    print "ERROR : propotion value below 0"
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
            print "ERROR : the sum of proportion given is greater than 100"
            exit(1)
        elif not percentage and tot_prop > 1:
            print "ERROR : the sum of proportion given is greater than 1"
            exit(1)
    except ValueError:
        print "ERROR : wrong proportion values."
        print "Exiting..."
        print exit(1)
    return new_list


def handling_nt_proportion(prop_A, prop_T, prop_C, prop_G):
    """
    Turning all string value to int, if possible. The None value are the estimated thanks to all the value in
    prop_A, prop_T, prop_C, prop_G. At the and the sum of prop_A, prop_T, prop_C, prop_G must be equal to 1
    :param prop_A: (string of a number or None value) proportion of A in the sequence we will create
    :param prop_T: (string of a number or None value) proportion of T in the sequence we will create
    :param prop_C: (string of a number or None value) proportion of C in the sequence we will create
    :param prop_G: (string of a number or None value) proportion of G in the sequence we will create
    :return:  prop[0], prop[1], prop[2], prop[3] float value
    """
    nbr_none = 0
    none_list = []
    prop = my_format([prop_A, prop_T, prop_C, prop_G])
    temp = my_format([prop_A, prop_T, prop_C, prop_G])
    for i in range(len(prop)):
        if prop[i] is None:
            nbr_none += 1
            none_list.append(i)
    if nbr_none == 0:
        return prop[0], prop[1], prop[2], prop[3]
    else:
        for i in range(len(none_list)):
            temp[none_list[i]] = 0
        for i in range(len(none_list)):
            prop[none_list[i]] = (1 - sum(temp)) / len(none_list)
        return prop[0], prop[1], prop[2], prop[3]


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
    parser.add_argument('--prop_A', dest='prop_A', help="the proportion of alanine in the fasta file",
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

    args = parser.parse_args()  # parsing arguments

    args.prop_A, args.prop_T, args.prop_C, args.prop_G = \
        handling_nt_proportion(args.prop_A, args.prop_T, args.prop_C, args.prop_G)

    print "Nucleotides proportion : "
    print "A : " + str(args.prop_A) + " - C : " + str(args.prop_C) + " - G : " + str(args.prop_G) + " - T : " + \
          str(args.prop_T)
    if args.prop_T == 0.25:
        if 0 <= args.prop_A + args.prop_C + args.prop_G <= 1:
            args.prop_T = 1. - (args.prop_A + args.prop_C + args.prop_G)
        else:
            print "negative value or not proportion value (i.e value between 0-1 or 0-100)"
            print "Exiting"
            exit(1)

    try:
        args.size_inf = int(args.size_inf)
        args.size_max = int(args.size_max)
        if args.size_inf > args.size_max:
            print "WARNING : maximum size value smaller than minimum size value"
            print "switching size value (min <=> max)"
            temp = args.size_inf
            args.size_inf = args.size_max
            args.size_max = temp
    except ValueError:
        print "ERROR : wrong size values."
        print "Exiting..."
        print exit(1)

    try:
        args.nbr_seq = int(args.nbr_seq)
        if args.nbr_seq < 0:
            int("a")
    except ValueError:
        print "ERROR : wrong 'nbr_seq' values"
        print "Exiting"
        exit(1)

    if not os.path.isdir(args.output):
        print "The given path in 'output' doesn't exist !"
        print "fasta file will be created in your current working directory"
        args.output = "./"

    if args.output[-1] != "/":
        args.output += "/"

    if args.flexible == "True":
        args.flexible = True
    elif args.flexible == "False":
        args.flexible = False
    else:
        print "WARNING : unrogonized boolean value for flexible argument"
        print "Setting it to False ! "
        args.flexible = False

    size_int = [args.size_inf, args.size_max]
    fasta_generator(size_int,  args.prop_A, args.prop_T, args.prop_C, args.prop_G, args.nbr_seq, args.output,
                    args.filename, args.flexible)

launcher()

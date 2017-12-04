"""Summary: this script wil create propensity chart of a set of exons."""

# imports
import pandas as pd
import numpy as np
import copy
from matplotlib import pyplot as plt
import argparse

# dictionary:
# hydrophobicity of amino acids - Source Composition profiler.
# Reference : Kyte J, and Doolittle RF. (1982)
# "A simple method for displaying the hydropathic character of a protein."
#  J. Mol. Biol. 157:105-132
aa2kyte_hydrophobicity = {
    "R": -4.5, "K": -3.9, "D": -3.5, "E": -3.5, "N": -3.5, "Q": -3.5,
    "H": -3.2, "P": -1.6, "Y": -1.3, "W": -0.9, "S": -0.8, "T": -0.7,
    "G": -0.4, "A": 1.8, "M": 1.9, "C": 2.5, "F": 2.8, "L": 3.8, "V": 4.2,
    "I": 4.5
}

# functions


def exons_reader(excel_file):
    """Description:read a query_result file produce by the tRNA program.

    :param excel_file: (string) path to an excel file.
        it must be produced by the tRNA program
    :return: (list of string) the list of each sequence in the mapping file.
    """
    sequences = []
    xl = pd.ExcelFile(excel_file)
    df = "NA"
    for sheet in xl.sheet_names:
        if "sequence" == sheet:
            df = xl.parse(sheet)
    # if the sheet nt info doesn't exist, we end the program
    if str(df) == "NA":
        print("The sheet sequence was not found in " + str(excel_file))
        print("Terminating...")
        exit(1)
    for row in df.itertuple():
        if isinstance(row.CDS_peptide_sequence, str):
            seq = row.CDS_peptide_sequence.replace("*", "")
            if len(seq) > 0:
                sequences.append(seq)
    return sequences


def cordinate_calculator(list_seq, dic):
    """Turn a list of sequence into cordinates given a dictionary dic.

    :param list_seq: (list of string) list of pepide sequences
    :param dic: (dict of float) associate each amino_acid to a value
    :return: 4 lists of floats : absicsa, ordinate, error, nb_seq
    """
    list_val = []
    max_length = 0
    nb_seq = []
    for i in range(len(list_seq)):
        if len(list_seq[i]) > max_length:
            max_length = len(list_seq[i])
    absissa = np.arange(max_length)
    for j in absissa:
        cur_val = []
        seq = 0
        for i in range(len(list_seq)):
            try:
                cur_val.append(dic[list_seq[i][j]])
                seq += 1
            except IndexError:
                seq += 0
        nb_seq.append(seq)
        list_val.append(copy.deep_copy(cur_val))
    ordinate = []
    error = []
    for i in range(len(list_val)):
        ordinate.append(np.mean(list_val[i]))
        error.append(np.sd(list_val[i]))
    return absissa, ordinate, error, nb_seq


def graphic_maker(exon_up_coord, exon_down_coord, name_scale, scale, output):
    """Create a graphic.

    :param exon_up_coord: tuple of 4 list of float/int
    :param exon_down_coord: tuple of 4 list of float/int
    """
    fig, ax = plt.subplots(figsize=(48. / 2.54, 27 / 2.54))
    absissa_up, ordinate_up, error_up, nb_seq_up = exon_up_coord
    absissa_down, ordinate_down, error_down, nb_seq_down = exon_down_coord
    if len(absissa_up) > len(absissa_down):
        absissa = absissa_up
        for i in range(len(absissa_up) - len(absissa_down)):
            ordinate_down.append(0)
            error_down.append(0)
            nb_seq_down.append(0)
    else:
        absissa = absissa_down
        for i in range(len(absissa_down) - len(absissa_up)):
            ordinate_up.append(0)
            error_up.append(0)
            nb_seq_up.append(0)
    label_up = "average " + str(name_scale) + " for up exons"
    label_down = "average " + str(name_scale) + " for down exons"
    area_up = "std of " + str(name_scale) + " for up exons"
    area_down = "std of " + str(name_scale) + " for down exons"
    ax.plot(absissa, ordinate_up, color="#59BADE", marker="-", label=label_up)
    ax.fill_between(absissa, ordinate_up - error_up, ordinate_up + error_up,
                    alpha=0.5, color="#59BADE", label=area_up)
    ax.plot(absissa, ordinate_down, color="#59BADE", marker="-",
            label=label_down)
    ax.fill_between(absissa, ordinate_down - error_down,
                    ordinate_down + error_down, alpha=0.5, color="#59BADE",
                    label=area_down)
    plt.title(name_scale + " by position of peptite coded by up/down exons")
    plt.xlabel("position in peptides")
    plt.ylabel(scale + " scale")
    plt.savefig(output + scale + "_figure.pdf", bbox_inches='tight')
    plt.savefig(output + scale + "_figure.pdf", bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


def wrap(excel_up, excel_down, output):
    """
    Wrap all the functions on top of this function.

    :param excel_up: (string) the name of the query result for the up exons
    :param excel_down: (string) the name of the query result for the down exons
    """
    dic = aa2kyte_hydrophobicity
    seq_up = exons_reader(excel_up)
    seq_down = exons_reader(excel_down)
    res_up = cordinate_calculator(seq_up, dic)
    res_down = cordinate_calculator(seq_down, dic)
    graphic_maker(res_up, res_down, "hydrophobicity",
                  "hydrophobicity(kyte-1982)", output)


def launcher():
    """Function that contains a parser to launch the program."""
    # description on how to use the program
    desc = """
    From 2 query result files produced by the  tRNA program produce a figure
    of kyte hydrophobicityfor those 2 sets of peptides within the files"""
    format_arg = argparse.RawDescriptionHelpFormatter
    usage = '%(prog)s --up query_up.xlsx --down query_down.xlsx'
    usage += ' [--output an outptut folder]'
    parser = argparse.ArgumentParser(formatter_class=format_arg,
                                     description=desc,
                                     usage=usage)
    # Arguments for the parser
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--up', dest='up', required=True,
                               help="file of the tRNA program that "
                               "contains the up exons")
    requiredNamed.add_argument('--down', dest='down', required=True,
                               help="file of the tRNA program that contains "
                               "the down exons")
    requiredNamed.add_argument('--output', dest='output', required=True,
                               help="the file where the graphic wille be "
                               "created")
    args = parser.parse_args()  # parsing arguments
    wrap(args.up, args.down, args.output)


if __name__ == "__main__":
    launcher()

"""Summary: this script wil create propensity chart of a set of exons.

Descritpion : Creation of a graphic separated in 3 part:
    - 1 : the value of a given propensity at each position of the
    up-regulated (red) and the down regulated set of sequences (blue).
    The sequences are mapped at their first base.
    - 2 :  a boxplot chart created with to sets of value:
        - sets up : containing for each up exons a value for a particular
        propensity
        - sets down : containing for each up exons a value for a particular
        propensity
    - 3 the number of exons (up and down) analyzed in function of their size

"""

# imports
import pandas as pd
import numpy as np
import copy
from matplotlib import pyplot as plt
import argparse
from scipy.stats import mannwhitneyu
from matplotlib import collections as mc

# dictionaries:
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
# Reference :
# Eisenberg D, Schwarz E, Komaromy M, and Wall R. (1984)
# "Analysis of membrane and surface protein sequences with
# the hydrophobic moment plot." J. Mol. Biol. 179:125-142.
aa2eisenberg_hydrophobicity = {
    "R": -2.53, "K": -1.5, "D": -0.9, "Q": -0.85, "N": -0.78, "E": -0.74,
    "H": -0.4, "S": -0.18, "T": -0.05, "P": 0.12, "Y": 0.26, "C": 0.29,
    "G": 0.48, "A": 0.62, "M": 0.64, "W": 0.81, "L": 1.06, "V": 1.08,
    "F": 1.19, "I": 1.38
}

# hydrophobicity of amino acids - Source Composition profiler
# Reference : Fauchere J-L, and Pliska VE. (1983)
# "Hydrophobic parameters pi of amino acid side chains from
# partitioning of N-acetyl-amino-acid amides." Eur. J. Med. Chem. 18:369-375.
aa2fauchere_hydrophobicity = {
    "R": -1.01, "K": -0.99, "D": -0.77, "E": -0.64, "N": -0.6, "Q": -0.22,
    "H": 0.13, "P": 0.72, "Y": 0.96, "W": 2.25, "S": -0.04, "T": 0.26, "G": 0.,
    "A": 0.31, "M": 1.23, "C": 1.54, "F": 1.79, "L": 1.7, "V": 1.22, "I": 1.8
}

# polarity of amino acids - Source Composition profiler
# Reference : Zimmerman JM, Eliezer N, and Simha R. (1968)
# "The characterization of amino acid sequences in
# proteins by statistical methods." J J. Theor. Biol. 21:170-201.
aa2zimmerman_polarity = {
    "A": 0, "G": 0, "I": 0.13, "L": 0.13, "V": 0.13, "F": 0.35, "M": 1.43,
    "C": 1.48, "P": 1.58, "Y": 1.61, "T": 1.66, "S": 1.67, "W": 2.1, "N": 3.38,
    "Q": 3.53, "K": 49.5, "D": 49.7, "E": 49.9, "H": 51.6, "R": 52
}

# polarity of amino acids - Source Composition profiler
# Reference : Grantham R. Science 185:862-864(1974).
aa2grantham_polarity = {
    "A": 8.100, "R": 10.500, "N": 11.600, "D": 13.000, "C": 5.500,
    "Q": 10.500, "E": 12.300, "G": 9.000, "H": 10.400, "I": 5.200, "L": 4.900,
    "K": 11.300, "M": 5.700, "F": 5.200, "P": 8.000, "S": 9.200, "T": 8.600,
    "W": 5.400, "Y": 6.200, "V": 5.900
}

# Conformational parameter for alpha helix. - source Protscale
# Reference - Deleage G., Roux B. Protein Engineering 1:289-294(1987).
aa2deleage_alpha = {
    "A": 1.489, "R": 1.224, "N": 0.772, "D": 0.924, "C": 0.966,
    "Q": 1.164, "E": 1.504, "G": 0.510, "H": 1.003, "I": 1.003, "L": 1.236,
    "K": 1.172, "M": 1.363, "F": 1.195, "P": 0.492, "S": 0.739, "T": 0.785,
    "W": 1.090, "Y": 0.787, "V": 0.990
}

# Normalized frequency for alpha helix.  - source protscale
# Reference :  Levitt M. -  Biochemistry 17:4277-4285(1978).
aa2levitt_alpha = {
    "A": 1.290, "R": 0.960, "N": 0.900, "D": 1.040, "C": 1.110, "Q": 1.270,
    "E": 1.440, "G": 0.560, "H": 1.220, "I": 0.970, "L": 1.300, "K": 1.230,
    "M": 1.470, "F": 1.070, "P": 0.520, "S": 0.820, "T": 0.820, "W": 0.990,
    "Y": 0.720, "V": 0.910
}

# Conformational parameter for alpha helix (computed from 29 proteins).
# source protscale
# Reference Chou P.Y., Fasman G.D. - Adv. Enzym. 47:45-148(1978).
aa2chou_alpha = {
    "A": 1.420, "R": 0.980, "N": 0.670, "D": 1.010, "C": 0.700, "Q": 1.110,
    "E": 1.510, "G": 0.570, "H": 1.000, "I": 1.080, "L": 1.210, "K": 1.160,
    "M": 1.450, "F": 1.130, "P": 0.570, "S": 0.770, "T": 0.830, "W": 1.080,
    "Y": 0.690, "V": 1.060
}

# Beta structure frequency (Nagano, 1973) - source composition profiler
# Reference : Nagano K. (1973)
# "Local analysis of the mechanism of protein folding. I.
# Prediction of helices, loops, and beta-structures from primary structure."
# J. Mol. Biol. 75:401-420.
aa2nagano_beta = {
    "E": 0.33, "R": 0.67, "N": 0.72, "P": 0.75, "S": 0.77, "K": 0.81,
    "H": 0.87, "D": 0.9, "G": 0.9, "A": 0.96, "Y": 1.07, "C": 1.13, "W": 1.13,
    "Q": 1.18, "T": 1.23, "L": 1.26, "M": 1.29, "F": 1.37, "V": 1.41, "I": 1.54
}


# Conformational parameter for beta-sheet. - source protscale
# Deleage G., Roux B. -  Protein Engineering 1:289-294(1987).
aa2deleage_beta = {
    "A": 0.709, "R": 0.920, "N": 0.604, "D": 0.541, "C": 1.191, "Q": 0.840,
    "E": 0.567, "G": 0.657, "H": 0.863, "I": 1.799, "L": 1.261, "K": 0.721,
    "M": 1.210, "F": 1.393, "P": 0.354, "S": 0.928, "T": 1.221, "W": 1.306,
    "Y": 1.266, "V": 1.965
}

# Conformational parameter for beta-sheet (computed from 29 proteins).
# source protparam
# Chou P.Y., Fasman G.D. - Adv. Enzym. 47:45-148(1978).
aa2chou_beta = {
    "A": 0.830, "R": 0.930, "N": 0.890, "D": 0.540, "C": 1.190, "Q": 1.100,
    "E": 0.370, "G": 0.750, "H": 0.870, "I": 1.600, "L": 1.300, "K": 0.740,
    "M": 1.050, "F": 1.380, "P": 0.550, "S": 0.750, "T": 1.190, "W": 1.370,
    "Y": 1.470, "V": 1.700
}

# Conformational parameter for beta-turn. source - protscale
# Deleage G., Roux B. Protein Engineering 1:289-294(1987).
aa2deleage_bturn = {
    "A": 0.788, "R": 0.912, "N": 1.572, "D": 1.197, "C": 0.965, "Q": 0.997,
    "E": 1.149, "G": 1.860, "H": 0.970, "I": 0.240, "L": 0.670, "K": 1.302,
    "M": 0.436, "F": 0.624, "P": 1.415, "S": 1.316, "T": 0.739, "W": 0.546,
    "Y": 0.795, "V": 0.387
}

# Normalized frequency for beta-turn.  - source protscale
# Levitt M. Biochemistry 17:4277-4285(1978).
aa2levitt_bturn = {
    "A": 0.770, "R": 0.880, "N": 1.280, "D": 1.410, "C": 0.810, "Q": 0.980,
    "E": 0.990, "G": 1.640, "H": 0.680, "I": 0.510, "L": 0.580, "K": 0.960,
    "M": 0.410, "F": 0.590, "P": 1.910, "S": 1.320, "T": 1.040, "W": 0.760,
    "Y": 1.050, "V": 0.470
}

# Conformational parameter for beta-turn (computed from 29 proteins).
#  Chou P.Y., Fasman G.D.  -  Adv. Enzym. 47:45-148(1978).
aa2chou_bturn = {
    "A": 0.660, "R": 0.950, "N": 1.560, "D": 1.460, "C": 1.190, "Q": 0.980,
    "E": 0.740, "G": 1.560, "H": 0.950, "I": 0.470, "L": 0.590, "K": 1.010,
    "M": 0.600, "F": 0.600, "P": 1.520, "S": 1.430, "T": 0.960, "W": 0.960,
    "Y": 1.140, "V": 0.500
}

# Coil propensity (Nagano, 1973)
# Nagano K. (1973) "Local analysis of the mechanism of protein folding.
# I. Prediction of helices, loops,
# and beta-structures from primary structure." J. Mol. Biol. 75:401-420.
aa2nagano_coil = {
    "F": 0.58, "M": 0.62, "L": 0.63, "A": 0.72, "E": 0.75, "H": 0.76, "I": 0.8,
    "Q": 0.81, "V": 0.83, "K": 0.84, "W": 0.87, "C": 1.01, "T": 1.03,
    "D": 1.04, "R": 1.33, "S": 1.34, "G": 1.35, "Y": 1.35, "N": 1.38, "P": 1.43
}


# Conformational parameter for coil.
# Deleage G., Roux B. - Protein Engineering 1:289-294(1987).
aa2deleage_coil = {
    "A": 0.824, "R": 0.893, "N": 1.167, "D": 1.197, "C": 0.953, "Q": 0.947,
    "E": 0.761, "G": 1.251, "H": 1.068, "I": 0.886, "L": 0.810, "K": 0.897,
    "M": 0.810, "F": 0.797, "P": 1.540, "S": 1.130, "T": 1.148, "W": 0.941,
    "Y": 1.109, "V": 0.772
}

list_dic = [aa2kyte_hydrophobicity, aa2eisenberg_hydrophobicity,
            aa2fauchere_hydrophobicity, aa2zimmerman_polarity,
            aa2grantham_polarity, aa2deleage_alpha, aa2levitt_alpha,
            aa2chou_alpha, aa2nagano_beta, aa2deleage_beta, aa2chou_beta,
            aa2deleage_bturn, aa2levitt_bturn, aa2chou_bturn, aa2nagano_coil,
            aa2deleage_coil]
scale_name = ["hydrophobicity", "hydrophobicity", "hydrophobicity",
              "polarity", "polarity", "alpha helix prediction",
              "alpha helix prediction", "alpha helix prediction",
              "beta helix prediction", "beta helix prediction",
              "beta helix prediction", "beta turn prediction",
              "beta turn prediction", "beta turn prediction",
              "coil prediction", "coil prediction"]

scale = ["hydrophobicity(kyte)", "hydrophobicity(eisenberg)",
         "hydrophobicity(fauchere)", "polarity(zimmerman)",
         "polarity(grantham)", "alpha helix prediction(deleage)",
         "alpha helix prediction(levitt)", "alpha helix prediction(chou)",
         "beta helix prediction(nagano)", "beta helix prediction(deleage)",
         "beta helix prediction(chou)", "beta turn prediction(deleage)",
         "beta turn prediction(levitt)", "beta turn prediction(chou)",
         "coil prediction(nagano)", "coil prediction(deleage)"]

# functions


def exons_reader(excel_file, max_size):
    """Description:read a query_result file produce by the tRNA program.

    :param excel_file: (string) path to an excel file.
        it must be produced by the tRNA program
    :param max_size: (int) the max size allowed for the peptide
    :return: (list of string) the list of each sequence in the mapping file.
    the sequences cannot have a size greater than max_size.
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
    for row in df.itertuples():
        if isinstance(row.CDS_peptide_sequence, unicode):
            seq = row.CDS_peptide_sequence.replace("*", "")
            if 0 < len(seq) < max_size + 1:
                sequences.append(seq)
    return sequences


def get_exons_value(list_seq, dic):
    """
    Calculate the propensity scale given in dic for the sequences in list_seq.

    Give for each sequence in list_seq, its value according to each amino acid
    values in dic.
    :param list_seq: (list of string), list of peptide sequences
    :param dic: (dictionary) each amino acid (key) is associated with a float
    value (value)
    :return: (list of float) the list of propensity value for each sequences
    in list_seq
    """
    list_val = []
    correction = True
    for key in dic.keys():
        if dic[key] < 0:
            correction = False
    for i in range(len(list_seq)):
        val = 0.
        for j in range(len(list_seq[i])):
            val += dic[list_seq[i][j]]
        if correction:
            list_val.append(val / len(list_seq[i]))
        else:
            list_val.append(val)
    return list_val


def cordinate_calculator(list_seq, dic):
    """Turn a list of sequence into coordinates given a dictionary dic.

    :param list_seq: (list of string) list of peptide sequences
    :param dic: (dict of float) associate each amino_acid to a value
    :return: 4 lists of floats : abscissa, ordinate, error, nb_seq
    and a list of list of float.
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
        list_val.append(copy.deepcopy(cur_val))
    ordinate = []
    error = []
    for i in range(len(list_val)):
        ordinate.append(np.mean(list_val[i]))
        error.append(np.std(list_val[i]))
    return absissa, ordinate, error, nb_seq, list_val


def make_man_withney_test(list_up, list_down):
    """
    Make a mann-whitney test.

    :param list_up: (list of float)
    :param list_down: (list of float)
    :return: (float) the p-value
    """
    try:
        p_val = mannwhitneyu(list_up, list_down)
    except ValueError:
        p_val = ["NA", "NA"]
    return p_val[1]


def line_maker(list_pval, up_mean, down_mean, up_value):
    """
    Create a list of lines, with their colour.

    Create a line if the a p-value in list_pval[i] is below 0.05.
    If up_mean[i] > down_mean[i] the line will be green, purple else.
    :param list_pval: (list of float), list of pvalue get by the comparison
    of a propensity scale in a particular sequence position in an
    up-regulated and down_regulated set of sequences.
    :param up_mean: (list of float) the mean propensity scale in all position
    of a regulated set of sequence.
    :param down_mean: (list of float) the mean propensity scale in all position
    of a regulated set of sequence.
    :param up_value: (int) the ordinate coordinates where the line will be
    placed.
    :return: lines - (list of list of 2 tuple), the list of 2 tuple corresponds
    to a lines with the coordinates [(x1, y1), (x2, y2)]
    """
    lcolor = []
    lines = []
    for i in range(len(list_pval)):
        if list_pval[i] < 0.05:
            lines.append([(i - 0.5, up_value), (i + 0.5, up_value)])
            if up_mean[i] > down_mean[i]:
                lcolor.append("#66FF66")  # green
            else:
                lcolor.append("#B266FF")  # purple

    return lines, lcolor


def graphic_maker(exon_up_coord, exon_down_coord, box_up, box_down, list_pval,
                  name_scale, scale_n, output):
    """Create the recap graphic.

    :param exon_up_coord: tuple of 4 list of float/int
    :param exon_down_coord: tuple of 4 list of float/int
    :param box_up: (list of floats) the propensity value of each up exons
    :param box_down: (list of floats) the propensity value of each down exons
    :param list_pval: (list of float), list of pvalue get by the comparison
    of a propensity scale in a particular sequence position in an
    up-regulated and down_regulated set of sequences.
    :param name_scale: (string) the short name of a scale
    :param scale_n: (string) the full name of a scale
    :param output: (string) the path where the figures will be created
    """
    fig = plt.figure(figsize=(48. / 2.54, 27 / 2.54))
    ax = fig.add_subplot(2, 1, 1)
    abscissa_up, ordinate_up, error_up, nb_seq_up = exon_up_coord
    abscissa_down, ordinate_down, error_down, nb_seq_down = exon_down_coord
    if len(abscissa_up) > len(abscissa_down):
        abscissa = abscissa_up
        for i in range(len(abscissa_up) - len(abscissa_down)):
            ordinate_down.append(0)
            error_down.append(0)
            nb_seq_down.append(0)
    else:
        abscissa = abscissa_down
        for i in range(len(abscissa_down) - len(abscissa_up)):
            ordinate_up.append(0)
            error_up.append(0)
            nb_seq_up.append(0)
    label_up = "average " + str(name_scale) + " for up exons"
    label_down = "average " + str(name_scale) + " for down exons"
    area_up = "std of " + str(name_scale) + " for up exons"
    area_down = "std of " + str(name_scale) + " for down exons"
    ax.plot(abscissa, ordinate_up, color="#EB4C4C", label=label_up)
    ord_1 = [x - y for x, y in zip(ordinate_up, error_up)]
    ord_2 = [x + y for x, y in zip(ordinate_up, error_up)]
    ax.fill_between(abscissa, ord_1, ord_2,
                    alpha=0.5, color="#EB4C4C", label=area_up)
    ax.plot(abscissa, ordinate_down, color="#59BADE",
            label=label_down)
    ord2_1 = [x - y for x, y in zip(ordinate_down, error_down)]
    ord2_2 = [x + y for x, y in zip(ordinate_down, error_down)]
    ax.fill_between(abscissa, ord2_1, ord2_2, alpha=0.5, color="#59BADE",
                    label=area_down)
    up_value = max(ord_2 + ord2_2)
    lines, lcolor = line_maker(list_pval, ordinate_up, ordinate_down, up_value)
    lc = mc.LineCollection(lines, colors=lcolor, linewidths=2)
    ax.add_collection(lc)
    ax.set_title(name_scale + " by position of peptide coded by up/down exons")
    ax.set_xlabel("position in peptides")
    ax.set_ylabel(scale_n + " scale by position")
    ax.plot([], [], color="#66FF66", label="up greater than down (p<0.05)")
    ax.plot([], [], color="#B266FF", label="down greater than up (p<0.05)")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, frameon=True, loc='lower center', ncol=3)

    ax2 = fig.add_subplot(2, 2, 3)
    lab = ["up", "down"]
    pval = make_man_withney_test(box_up, box_down)
    box = ax2.boxplot([box_up, box_down], labels=lab, showmeans=True,
                      notch=False, patch_artist=True)
    box['boxes'][0].set_facecolor("#EB4C4C")
    box['boxes'][1].set_facecolor("#59BADE")
    title = name_scale + " by peptides coded by up/down exons"
    ytitle = scale_n + "scale by peptide"
    ax2.set_title(title)
    ax2.set_xlabel("man_whitney test : p=" + str(pval))
    ax2.set_ylabel(ytitle)

    ax3 = fig.add_subplot(2, 2, 4)
    ax3.plot(abscissa, nb_seq_up, color="#EB4C4C", label="up")
    ax3.plot(abscissa, nb_seq_down, color="#59BADE", label="down")
    title = "number of amino acid studied at each peptide position"
    ytitle = "number of amino acid studied"
    xtitle = "position"
    ax3.set_title(title)
    ax3.set_ylabel(ytitle)
    ax3.set_xlabel(xtitle)
    handles, labels = ax3.get_legend_handles_labels()
    ax3.legend(handles, labels)

    plt.savefig(output + scale_n + "_figure.pdf", bbox_inches='tight')
    plt.savefig(output + scale_n + "_figure.pdf", bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()


def make_list_comparison(list_val_up, list_val_down):
    """
    Make multiple mann whitney tests.

    :param list_val_up: (list of list of floats)
    :param list_val_down: (list of lists of floats)
    :return: list of pvalue at each position of up and down regulated sequences
    """
    list_pval = []
    min_size = min(len(list_val_up), len(list_val_down))
    for i in range(min_size):
        if len(list_val_up) > 5 and len(list_val_down) > 5:
            list_pval.append(make_man_withney_test(list_val_up[i],
                                                   list_val_down[i]))
        else:
            break
    return list_pval


def wrap(excel_up, excel_down, output, max_size):
    """
    Wrap all the functions on top of this function.

    :param excel_up: (string) the name of the query result for the up exons
    :param excel_down: (string) the name of the query result for the down exons
    :param output: (string) the path  where the figure will be created
    :param max_size: (int) the maximum size of peptide sequence allowed.
    if there are longer than max_size, they will not be kept for the graphics
    """
    for i in range(len(list_dic)):
        seq_up = exons_reader(excel_up, max_size)
        seq_down = exons_reader(excel_down, max_size)
        exon_value_up = get_exons_value(seq_up, list_dic[i])
        exon_value_down = get_exons_value(seq_down, list_dic[i])
        res_up = cordinate_calculator(seq_up, list_dic[i])
        res_down = cordinate_calculator(seq_down, list_dic[i])
        list_pval = make_list_comparison(res_up[-1], res_down[-1])
        res_up = res_up[:-1]
        res_down = res_down[:-1]
        graphic_maker(res_up, res_down, exon_value_up, exon_value_down,
                      list_pval, scale_name[i], scale[i], output)


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
    required_named = parser.add_argument_group('required arguments')
    required_named.add_argument('--up', dest='up', required=True,
                                help="file of the tRNA program that "
                                "contains the up exons")
    required_named.add_argument('--down', dest='down', required=True,
                                help="file of the tRNA program that contains "
                                "the down exons")
    required_named.add_argument('--output', dest='output', required=True,
                                help="the file where the graphic wille be "
                                "created")
    parser.add_argument("--max_size", dest="max_size", default=100000,
                        help="The max size of sequences in the file")
    args = parser.parse_args()  # parsing arguments

    try:
        args.max_size = int(args.max_size)
    except ValueError:
        print("Wrong size value")
        print("Exiting...")
        exit(1)

    wrap(args.up, args.down, args.output, args.max_size)


if __name__ == "__main__":
    launcher()

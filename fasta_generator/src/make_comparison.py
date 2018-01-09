"""Description: Create an enrichment file that will compare
a fasta file to a desired control"""

######################################
#                IMPORT              #
######################################

import argparse
from Bio import SeqIO
from dicitonary import *
from scipy.stats import hypergeom
import sys
import os
import re
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import xlsxwriter
import math

##################################################
#                 Function
##################################################


def fasta_reader(fasta_file):
    """
    Read a fasta file.

    :param fasta_file: (string) a fasta file
    :return: (list of string) list of all sequence in the fasta file
    """
    list_seq = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        list_seq.append(str(record.seq))
    return list_seq


def calcul_dic(dic, seq):
    """
    Add the count of hexanucleotide in the sequence 'seq'
    to a dictionary 'dic' that already contains the hexanucleotid count of mainy previous
    nucleotide sequence

    :param dic: (a dictionary of float) freq of the 6-length word of interest
    :param seq: (string) the nt sequence
    :return: 'dic' completed
    """
    for i in range(len(seq)-(6-1)):
        if seq[i:i+6] not in dic.keys():
            dic[seq[i:i+6]] = 1
        else:
            dic[seq[i:i+6]] += 1
    return dic


def create_an_hexanucleotid_dic(list_seq):
    """
    Count the total number of hexa-nucleotide in a list
    if sequence list_seq

    :param list_seq: (list of string) list of DNA sequence
    :return: dic, count
     -dic (dictionary of int) the number of each hexa-nucleotide in
     the list of sequence
     -count : the total number of hexa-nucleotides
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    mod = __import__("CCE_dic")

    dic = {}
    for i in range(len(list_seq)):
        print(i)
        dic = calcul_dic(dic, list_seq[i])
    count = 0
    for key in dic.keys():
        count += dic[key]
    for key in mod.d6.keys():
        if key not in dic.keys() and key != "all":
            dic[key] = 0
    return dic, count


def nucleotide_calculator(dic, seq):
    """
    Add the count of di-nucleotide in the sequence 'seq'
    to a dictionary 'dic' that already contains the di_nucleotide count of many previous
    nucleotide sequences

    :param dic: (dictionary of int) the number of amino acid for a
    given set of sequence
    :param seq: (string) the current exons (cds) sequence studied
    :return: (dictionary) a dictionary containing the frequency of every possible di-nucleotides
    """
    if len(seq) > 1:

        for j in range(len(seq) - 1):
            dic[seq[j]] += 1
    return dic


def create_a_nt_dic(list_seq):
    """
    Count the total number of di-nucleotide in a list
    if sequence list_seq

    :param list_seq:(list of string) List of sequence in the fasta file
    :return: a dictionary that contains the number of dnt
    in all exons in tuple list
    """
    dic = {"A": 0, "T": 0, "G": 0, "C": 0}
    for i in range(len(list_seq)):
        dic = nucleotide_calculator(dic, "".join(list_seq[i]))
    count = 0
    for key in dic.keys():
        count += dic[key]
    dic["Y"] = dic["C"] + dic["T"]
    dic["R"] = dic["A"] + dic["G"]
    dic["W"] = dic["A"] + dic["T"]
    dic["S"] = dic["C"] + dic["G"]
    dic["K"] = dic["T"] + dic["G"]
    dic["M"] = dic["C"] + dic["A"]
    return dic, count


def dinucleotide_calculator(dic, seq):
    """
    Add the count of di-nucleotide in the sequence 'seq'
    to a dictionary 'dic' that already contains the di_nucleotide count of many previous
    nucleotide sequences

    :param dic: (dictionary of int) the number of amino acid for a
    given set of sequence
    :param seq: (string) the current exons (cds) sequence studied
    :return: (dictionary) a dictionary containing the frequency of every possible di-nucleotides
    """
    if len(seq) > 1:

        for j in range(len(seq) - 1):
            dic[seq[j:j + 2]] += 1
    return dic


def create_a_dnt_dic(list_seq):
    """
    Count the total number of di-nucleotide in a list
    if sequence list_seq

    :param list_seq:(list of string) List of sequence in the fasta file
    :return: a dictionary that contains the number of dnt
    in all exons in tuple list
    """
    dic = {"AA": 0, "AT": 0, "AG": 0, "AC": 0, "TA": 0, "TT": 0, "TG": 0, "TC": 0,
           "GA": 0, "GT": 0, "GG": 0, "GC": 0, "CA": 0, "CT": 0, "CG": 0, "CC": 0}
    for i in range(len(list_seq)):
        dic = dinucleotide_calculator(dic, "".join(list_seq[i]))
    count = 0
    for key in dic.keys():
        count += dic[key]
    return dic, count


def calcul_dic_codon(dic, seq):
    """
    Add the count of codon in the sequence 'seq'
    to a dictionary 'dic' that already contains the codon count of many previous
    nucleotide sequences

    :param dic: (dictionary of int) the number of codon for a
    given set of sequences
    :param seq: (string) nucleotides sequence
    :return: dic with the updated codon content
    """
    for i in range(0, len(seq)-2, 3):
        if seq[i:i+3] not in dic.keys():
            dic[str(seq[i:i+3])] = 1
        else:
            dic[str(seq[i:i+3])] += 1
    return dic


def create_a_codon_dic(list_seq):
    """
    Count the total number of codon in a list
    if sequence list_seq

    :param list_seq: List of sequence in the fasta file
    :return:dic, count
        - dic : (dictionary of int) the number of each codon in the
            list of sequence given by list_seq
        - count : the total number of codon in the list of sequence
        in list_seq
    """
    dic = {}
    for i in range(len(list_seq)):
        print(i)
        dic = calcul_dic_codon(dic, list_seq[i])
    count = 0
    for key in dic.keys():
        count += dic[key]
    return dic, count


def calcul_dic_aa(dic, seq):
    """
    Add the count of amino acid encoded in a nucleotide sequence 'seq'
    to a dictionary 'dic' that already contains the encoded amino acid count of many previous
    nucleotide sequences

    :param dic: (dictionary of int) the number of amino acids for a
    given set of sequences
    :param seq: (string) a nucleotides sequence
    :return: dic with the updated amino acids content
    """
    for i in range(0, len(seq)-2, 3):
        if codon2aminoAcid[seq[i:i+3]] != "":
            if codon2aminoAcid[seq[i:i+3]] not in dic.keys():
                dic[codon2aminoAcid[seq[i:i+3]]] = 1
            else:
                dic[codon2aminoAcid[seq[i:i+3]]] += 1
    return dic


def create_an_aa_dic(list_seq):
    """
    Count the total number of amino acid in a list
    if sequence list_seq

    :param list_seq: List of sequence in the fasta file
    :return:dic, count
        - dic : (dictionary of int) the number of each amino acid in the
            list of sequence given by list_seq
        - count : the total number of amino acid in the list of sequence
        in list_seq
    """
    dic = {}
    for i in range(len(list_seq)):
        if len(list_seq[i]) > 0:
            dic = calcul_dic_aa(dic, list_seq[i])
    count = 0
    for key in dic.keys():
        count += dic[key]
    return dic, count


def create_a_custom_dic(aa_dic, feature_dic):
    """
    Count the total number of amino acid that participate to a particular
    protein feature (from a dictionary of count encoded amino acids in a list of sequence)

    :param aa_dic:  a dictionary that contains the number of every amino acids
    in all sequence in the fasta file
    :param feature_dic: (dictionary of list of character) link each feature
    to their corresponding amino acids
    :return: a dictionary that link for each amino acid feature their number in the
    given fasta file
    """
    res = {}
    for key in feature_dic.keys():
        res[key] = 0
        for aa in feature_dic[key]:
            if aa in aa_dic.keys():
                res[key] += aa_dic[aa]
    return res


# def get_exons_value(list_seq, dic):
#     """
#     Calculate the propensity scale given in dic for the sequences in list_seq.
#
#     Give for each sequence in list_seq, its value according to each amino acid
#     values in dic.
#     :param list_seq: (list of string), list of peptide sequences
#     :param dic: (dictionary) each amino acid (key) is associated with a float
#     value (value)
#     :return: (list of float) the list of propensity value for each sequences
#     in list_seq
#     """
#     list_val = []
#     correction = True
#     for key in dic.keys():
#         if dic[key] < 0:
#             correction = False
#     for i in range(len(list_seq)):
#         val = 0.
#         if len(list_seq[i]) > 0:
#             for j in range(len(list_seq[i])):
#                 val += dic[list_seq[i][j]]
#             if correction:
#                 list_val.append(val / len(list_seq[i]))
#             else:
#                 list_val.append(val)
#     return list_val


def create_a_codon_pos_dic(codon_dic):
    """
    Count the nucleotide frequency at each codon position from a dic
    of codon.

    :param codon_dic: (dictionary of int) contains the codon count of a list of sequence
    :return: (dictionary of int) the nucleotide frequency at each codon position
    """
    nt_pos_dic = {"A1": 0, "A2": 0, "A3": 0, "C1": 0, "C2": 0, "C3": 0,
                  "G1": 0, "G2": 0, "G3": 0, "T1": 0, "T2": 0, "T3": 0}
    count = 0
    for i in range(0, 3, 1):
        for codon in codon_dic.keys():
            if codon != "all":
                count += codon_dic[codon]
                nt_pos_dic[codon[i] + str(i+1)] += codon_dic[codon]

    return nt_pos_dic, count / 3


def create_dic(fasta_file):
    """
    Create the count dictionary for the hexanucleotide, di-nucleotide, codon,
    codon position nucleotides, amino acids, and other (i.e. group of amino acids
    that encode for a particular type of protein feature)

    :param fasta_file:  (string) a fasta file
    :return: res_dic, [all_hexa, all_dnt, all_codon, all_codon_pos, all_aa]
         - res_dic : (list of 10 dictionaries of int)
         this list of dictionaries contains : the count dictionary for the hexanucleotide,
         di-nucleotide, codon,codon position nucleotides, amino acids, and other
         (i.e. group of amino acids that encode for a particular type of protein feature)
         - [all_hexa, all_dnt, all_codon, all_codon_pos, all_aa] : list of int:
         the total number of hexa-nucleotide, di-nucleotide, codons, nucleotide for the
         3 different codon position, and encoded amino acids of every sequence in the fasta file

    """
    res_dic = []
    print("reading fasta")
    list_seq = fasta_reader(fasta_file)
    print("generation of hex dics")
    dic, all_hexa = create_an_hexanucleotid_dic(list_seq)
    res_dic.append(dic)
    print("generation of di-nucleotide dictionary")
    dic, all_dnt = create_a_dnt_dic(list_seq)
    res_dic.append(dic)
    print("generation of nucleotides dictionary")
    dic, all_nt = create_a_nt_dic(list_seq)
    res_dic.append(dic)
    print("generation of codon dics")
    dic, all_codon = create_a_codon_dic(list_seq)
    res_dic.append(dic)
    print("generation of codon pos dics")
    dic, all_codon_pos = create_a_codon_pos_dic(dic)
    res_dic.append(dic)
    print("generation of aa dics")
    dic, all_aa = create_an_aa_dic(list_seq)
    res_dic.append(dic)
    print("generating group dics")
    res_dic.append(create_a_custom_dic(res_dic[4], schain2aa))
    res_dic.append(create_a_custom_dic(res_dic[4], hydro_info2aa))
    res_dic.append(create_a_custom_dic(res_dic[4], charge_info2aa))
    res_dic.append(create_a_custom_dic(res_dic[4], polarity_info2aa))
    res_dic.append(create_a_custom_dic(res_dic[4], misc2aa))
    """
    print("generating propensity dics...")
    list_dic = [aa2kyte_hydrophobicity, aa2eisenberg_hydrophobicity,
                aa2fauchere_hydrophobicity, aa2zimmerman_polarity,
                aa2grantham_polarity, aa2deleage_alpha, aa2levitt_alpha,
                aa2chou_alpha, aa2nagano_beta, aa2deleage_beta, aa2chou_beta,
                aa2deleage_bturn, aa2levitt_bturn, aa2chou_bturn, aa2nagano_coil,
                aa2deleage_coil]
    scale = ["hydrophobicity_kyte", "hydrophobicity_eisenberg",
             "hydrophobicity_fauchere", "polarity_zimmerman",
             "polarity_grantham", "alpha_helix_prediction_deleage",
             "alpha_helix_prediction_levitt", "alpha_helix_prediction_chou",
             "beta_helix_prediction_nagano", "beta_helix_prediction_deleage",
             "beta_helix_prediction_chou", "beta_turn_prediction_deleage",
             "beta_turn_prediction_levitt", "beta_turn_prediction_chou",
             "coil_prediction_nagano", "coil_prediction_deleage"]
    for i in range(len(list_dic)):
        print(scale[i])
        res_dic.append(get_exons_value(list_seq, list_dic[i]))
    """
    return res_dic, [all_hexa, all_dnt, all_nt, all_codon, all_codon_pos, all_aa]


def get_dic_pvalue(dic_fasta, dic_control, ctrl_all, fasta_all):
    """
    Calculate the p-value with an hypergeomtric test of the count
    of the hexanucleotide, di-nucleotide, codon, codon position nucleotides,
    amino acids, or others in a control set of exons and the exons in the fasta file

    :param dic_fasta:(dictionary of int) the count of the hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, or others in the fasta file given by the user
    :param dic_control:(dictionary of int) the count of the hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, or others in the control exons (exons CCE/ALL/ACE of fasterDB)
    :param ctrl_all: (int) the total number of the hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, or others in the control set of exons
    :param fasta_all:(int) the total number of the hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, or others in the fasta file
    :return: (dictionary of float) for each element give it's calculated p-value
    """
    p_val = {}
    for key in dic_fasta.keys():
        p1 = hypergeom.cdf(dic_fasta[key], ctrl_all, dic_control[key], fasta_all)
        p2 = 1 - hypergeom.cdf(dic_fasta[key], ctrl_all, dic_control[key], fasta_all)
        if p1 < p2:
            p_val[key] = p1
        else:
            p_val[key] = p2
    return p_val


# def get_propensity_pvalue(dic_fasta, dic_control):
#     p_val = {}
#     for key in dic_fasta.keys():
#         p_val[key] = mannwhitneyu(dic_fasta[key], dic_control[key])
#     return  p_val


def calculate_random_dic():
    """
    Create random dictionary of count of hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, or others elements.

    :return:
    d6_r : dictionary of count of hexanuleotide if they where randomly distributed in exons.
    ddnt_r : dictionary of count of hdi-nuleotide if they where randomly distributed in exons.
    dc_r: dictionary of count of di-nuleotide if they where randomly distributed in exons.
    dcp_r: dictionary of count of codon nucleotide position if they where randomly distributed in exons.
    da_r: dictionary of count of encoded amino acid if they where randomly distributed in exons
    sh, hy,ch, po, mi: dictionaries of count of amino acid encoding particular feature if they where randomly
     distributed in exons
    """
    mod = __import__("CCE_dic")
    codon_number = mod.dc["all"]
    dnt_number = mod.ddnt["all"]
    hexa_number = mod.d6["all"]
    nt_number = mod.nt["all"]
    d6_r = {}
    ddnt_r = {}
    nt_r = {}
    dc_r = {}
    da_r = {}
    dcp_r = {}
    for key in mod.d6.keys():
        if key != "all":
            d6_r[key] = int(math.pow(0.25, 6) * hexa_number)
    for key in mod.ddnt.keys():
        if key != "all":
            if key in ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]:
                ddnt_r[key] = int(math.pow(0.25, 2) * dnt_number)
            else:
                ddnt_r[key] = int(math.pow(0.5, 2) * dnt_number)
    for key in mod.nt.keys():
        if key in ["A", "T", "G", "C"]:
            nt_r[key] = int(0.25 * nt_number)
        elif key in ["S", "W", "K", "M", "Y", "R"]:
            nt_r[key] = int(0.5 * nt_number)
    for key in mod.dc.keys():
        if key != "all":
            dc_r[key] = int(math.pow(0.25, 3) * codon_number)
    for key in mod.dcp.keys():
        if key != "all":
            dcp_r[key] = int(math.pow(0.25, 3) * codon_number)
    print(int(math.pow(0.25, 3) * codon_number * 3))
    print(codon_number)
    for key in dc_r.keys():
        if codon2aminoAcid[key] != "":
            if codon2aminoAcid[key] not in da_r.keys():
                da_r[codon2aminoAcid[key]] = dc_r[key]
            else:
                da_r[codon2aminoAcid[key]] += dc_r[key]
    sh = create_a_custom_dic(da_r, schain2aa)
    hy = create_a_custom_dic(da_r, hydro_info2aa)
    ch = create_a_custom_dic(da_r, charge_info2aa)
    po = create_a_custom_dic(da_r, polarity_info2aa)
    mi = create_a_custom_dic(da_r, misc2aa)
    count = 0
    for key in d6_r.keys():
        count += d6_r[key]
    d6_r['all'] = count
    count = 0
    for key in ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]:
        count += ddnt_r[key]
    ddnt_r['all'] = count
    nt_r['all'] = nt_r["A"] + nt_r["T"] + nt_r["G"] + nt_r["C"]
    count = 0
    for key in dc_r.keys():
        count += dc_r[key]
    dc_r['all'] = count
    count = 0
    for key in ["A1", "C1", "G1", "T1", "A2", "C2", "G2", "T2", "A3", "C3", "G3", "T3"]:
        count += dcp_r[key]
    dcp_r['all'] = count/3
    count = 0
    for key in da_r.keys():
        count += da_r[key]
    da_r['all'] = count
    return [d6_r, ddnt_r, nt_r, dc_r, dcp_r, da_r, sh, hy, ch, po, mi]


def calculate_all_p_value_dic(fasta_dics, exon_type, all_list):
    """
    Calculate the p-value between the count of hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, and others elements in a control set and in the
    fasta set of exons

    :param fasta_dics: (list of dictionaries of in) dictionaries of count of hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, and others elements in the sequence of the fasta file.
    :param exon_type: (string) name of the control set of exons (CCE/ACE/ALL)
    :param all_list: (list of int) the total number of hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, and others elements.
    :return: (list of dictionaries of float) : p-value for hexanucleotide, di-nucleotide,
    codon, codon position nucleotides, amino acids, and others
    """
    print("calculation of pvalues")
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    p_val_dics = []
    if exon_type == "RD":
        ctrl_dics = calculate_random_dic()
    else:
        mod = __import__(exon_type + "_dic")
        ctrl_dics = [mod.d6, mod.ddnt, mod.nt, mod.dc, mod.dcp, mod.da, mod.sh, mod.hy, mod.ch, mod.po, mod.mi]
    """
    ctrl_dics = [mod.d6, mod.dc, mod.da, mod.sh, mod.hy, mod.ch, mod.po, mod.mi,
                mod.hydrophobicity_kyte, mod.hydrophobicity_eisenberg,
                mod.hydrophobicity_fauchere, mod.polarity_zimmerman,
                mod.polarity_grantham, mod.alpha_helix_prediction_deleage,
                mod.alpha_helix_prediction_levitt, mod.alpha_helix_prediction_chou,
                mod.beta_helix_prediction_nagano, mod.beta_helix_prediction_deleage,
                mod.beta_helix_prediction_chou, mod.beta_turn_prediction_deleage,
                mod.beta_turn_prediction_levitt, mod.beta_turn_prediction_chou,
                mod.coil_prediction_nagano, mod.coil_prediction_deleage]
    """
    for i in range(len(ctrl_dics)):
        print(str(i) + "dics out of " + str(11))
        if i < 6:
            p_val_dics.append(get_dic_pvalue(fasta_dics[i], ctrl_dics[i], ctrl_dics[i]["all"], all_list[i]))
        elif i < 11:
            p_val_dics.append(get_dic_pvalue(fasta_dics[i], ctrl_dics[i], ctrl_dics[5]["all"], all_list[5]))
        """
        elif i >= 8:
            p_val_dics.append(get_propensity_pvalue(fasta_dics[i], ctrl_dics[i]))
        """
    return p_val_dics, ctrl_dics


def find_hexa_motif(motif, hexa_fasta):
    """
    Give a dictionary that links each hexa nucleotide to the value
    'motif' if 'motif' is in the heaxnucleotide

    :param motif: (string) a regular expression or a simple nucleotid motif
    :param hexa_fasta: (dictionary of int) the dictionary that count each hexanucleotide
    in the fasta file.
    :return: motif_dic - (dictionary of string) inks each hexa nucleotide to the value
    'motif' if 'motif' is in the heaxnucleotide
    """
    motif_dic = {}
    for hexant in hexa_fasta.keys():
        res = re.findall(str(motif), str(hexant))
        if len(res) > 0:
            motif_dic[hexant] = "MOTIF"
        else:
            motif_dic[hexant] = "-"
    return motif_dic


def get_p_adjust_dic(pval_dic):
    """
    Correct the p-values by the benjamini hotchberg method
    :param pval_dic: (dictionary of float) the float values given in the dictionary are the p-value to correct
    :return: (dictionary of float) the corrected p-values
    """
    p_adj = {}
    list_key = pval_dic.keys()
    list_pval = []
    for key in list_key:
        list_pval.append(pval_dic[key])
    rstats = importr('stats')
    p_adjust_codon = rstats.p_adjust(FloatVector(list_pval), method="BH")
    for i in range(len(list_key)):
        p_adj[list_key[i]] = p_adjust_codon[i]
    return p_adj


def get_hexanucleotide_content(hexa_ctrl, hexa_fasta, pval_hexa, exon_type, motif, ctrl_all, fasta_all):
    """
    Create a list of list of string variable that represent the content of then sheet hexanucleotide
    in the enrichment report file.

    :param hexa_ctrl: (dictionary of int) dictionary of count of hexanucleotide in the control set of exons
    :param hexa_fasta: (dictionary of int) dictionary of count of hexanucleotide in the fasta file
    :param pval_hexa:(dictionary of float) dictionary of p-value associated to each hexanucleotide
    :param exon_type: (string) the name of the control set of exons used
    :param motif: (string) a regular expression or a simple nucleotid motif
    :param ctrl_all: (int) the total number of hexanucleotid in the control set of exons
    :param fasta_all: (int) the total number of hexanucleotid in the fasta file
    :return: list of list of string variable that represent the content of then sheet hexanucleotide
    in the enrichment report file.
    """
    content = [["Hexanucleotides", "Motif", "fasta_content", exon_type + "_content", "pvalues", "p_corr", "reg_pval",
                "reg_p_cor"]]
    sorted_res = sorted(pval_hexa.items(), key=lambda l: l[1], reverse=False)
    motif_dic = find_hexa_motif(motif, hexa_fasta)
    p_adj = get_p_adjust_dic(pval_hexa)
    for i in range(len(sorted_res)):
        key = sorted_res[i][0]
        reg = " = "
        reg_cor = " = "
        if pval_hexa[key] < 0.05:
            if float(hexa_fasta[key])/fasta_all > float(hexa_ctrl[key])/ctrl_all:
                reg = " + "
            else:
                reg = " - "
        if p_adj[key] < 0.05:
            reg_cor = reg
        content.append([key, motif_dic[key], float(hexa_fasta[key])/fasta_all,
                        float(hexa_ctrl[key])/ctrl_all, pval_hexa[key], p_adj[key],
                        reg, reg_cor])
    return content


def get_codon_content(codon_ctrl, codon_fasta, pval_codon, exon_type, ctrl_all, fasta_all):
    """

    Create a list of list of string variable that represent the content of then sheet codon
    in the enrichment report file.

    :param codon_ctrl: (dictionary of int) dictionary of count of codon in the control set of exons
    :param codon_fasta: (dictionary of int) dictionary of count of codon in the fasta file
    :param pval_codon: (dictionary of float) dictionary of p-value associated to each codon
    :param exon_type: (string) the name of the control set of exons used
    :param ctrl_all:  (int) the total number of codon in the control set of exons
    :param fasta_all:  (int) the total number of codons in the fasta file
    :return: list of list of string variable that represent the content of then sheet codon
    in the enrichment report file.
    """
    content = [["codon", "aa", "fasta_content", exon_type + "_content", "pvalues", "p_corr", "reg_pval",
                "reg_p_cor"]]
    sorted_res = sorted(pval_codon.items(), key=lambda l: l[1], reverse=False)
    p_adj = get_p_adjust_dic(pval_codon)
    for i in range(len(sorted_res)):
        key = sorted_res[i][0]
        reg = " = "
        reg_cor = " = "
        if pval_codon[key] < 0.05:
            if float(codon_fasta[key]) / fasta_all > float(codon_ctrl[key]) / ctrl_all:
                reg = " + "
            else:
                reg = " - "
        if p_adj[key] < 0.05:
            reg_cor = reg
        content.append([key, codon2aminoAcid[key], float(codon_fasta[key]) /
                        fasta_all, float(codon_ctrl[key]) / ctrl_all, pval_codon[key], p_adj[key],
                        reg, reg_cor])
    return content


def get_content(ctrl, fasta, pval, exon_type, type_ft, ctrl_all, fasta_all):
    """
    Create a list of list of string variable that represent the content of the sheet our_group
    in the enrichment report file.

    :param ctrl:(dictionary of int) dictionary of count of a particular feature in the control set of exons
    :param fasta:(dictionary of int) dictionary of count of a particular feature in the fasta file
    :param pval: (dictionary of float) dictionary of p-value associated to each a particular feature
    :param exon_type:(string) the name of the control set of exons used
    :param type_ft: (string) the name of the particular feature
    :param ctrl_all:(int) the total number of a particular feature in the control set of exons
    :param fasta_all:(int) the total number of a particular feature in the fasta file
    :return:
    """
    content = [[type_ft, "fasta_content", exon_type + "_content", "pvalues", "p_corr", "reg_pval",
                "reg_p_cor"]]
    sorted_res = sorted(pval.items(), key=lambda l: l[1], reverse=False)
    p_adj = get_p_adjust_dic(pval)
    for i in range(len(sorted_res)):
        key = sorted_res[i][0]
        reg = " = "
        reg_cor = " = "
        if pval[key] < 0.05:
            if float(fasta[key]) / fasta_all > float(ctrl[key]) / ctrl_all:
                reg = " + "
            else:
                reg = " - "
        if p_adj[key] < 0.05:
            reg_cor = reg
        content.append([key, float(fasta[key]) / fasta_all,
                        float(ctrl[key]) / ctrl_all, pval[key], p_adj[key], reg, reg_cor])
    return content


def size_adaptater(content):
    """
    :param content: (list of list of string) the content of a sheet
    :return: (list of int), list of column size : it will allow to set the appropriate size of each column in all
    the sheets of the query_results.xlsx
    """
    val = 0
    for row in content:
        if len(row) > val:
            val = len(row)
    res = [0] * val
    for row in content:
        for i in range(len(row)):
            if res[i] < len(str(row[i])):
                res[i] = len(str(row[i]))
    return res


def sheet_filler(content, a_sheet, header_format, normal_format):
    """
    Fills the sheet 'a_sheet' with the content 'content'
    :param content: (list of list of string) the content of a sheet
    :param a_sheet: (instance of xlsxwriter.worksheet.Worksheet) the sheet you want to fill
    :param header_format: (instance of xlsxwriter.format.Format) the format of the header line
    :param normal_format: (instance of xlsxwriter.format.Format) the format of normal lines
    """
    i = 0
    for row in content:
        if i == 0 or 'Name : ' in row[0]:
            a_sheet.write_row(i, 0, row, header_format)
        else:
            a_sheet.write_row(i, 0, row, normal_format)
        i += 1
    cell_size = size_adaptater(content)
    for i in range(len(cell_size)):
        a_sheet.set_column(i, i, cell_size[i]+2)


def creating_report(hexa_content, dnt_content, nt_content, codon_content, codon_pos_content, aa_content, group_content, output):
    """
    Create the enrichment report file

    :param hexa_content: list of list of string variable that represent the content of then sheet hexanucleotide
    in the enrichment report file.
    :param dnt_content: list of list of string variable that represent the content of then sheet di-nucleotides
    in the enrichment report file.
    :param codon_content: list of list of string variable that represent the content of then sheet codon
    in the enrichment report file.
    :param codon_pos_content: list of list of string variable that represent the content of then sheet codon_pos
    in the enrichment report file.
    :param aa_content: list of list of string variable that represent the content of then sheet amino acid
    in the enrichment report file.
    :param group_content: list of list of string variable that represent the content of then sheet our_group
    in the enrichment report file.
    :param output:(string) the path where the result will be created
    """
    workbook = xlsxwriter.Workbook(output + "enrichment_report.xlsx")
    header_format = workbook.add_format({'bg_color': '#00DCFF', 'align': 'center', 'valign': 'vcenter', 'border': True})
    normal_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': True})

    nt_sheet = workbook.add_worksheet("nt")
    dnt_sheet = workbook.add_worksheet("dnt")
    hexa_sheet = workbook.add_worksheet("hexanucleotide")
    codon_sheet = workbook.add_worksheet("codon")
    codon_pos_sheet = workbook.add_worksheet("codon position")
    aa_sheet = workbook.add_worksheet("amino_acid")
    grp_sheet = workbook.add_worksheet("our_group")
    # prop_sheet = workbook.add_worksheet("propensity")

    sheet_filler(nt_content, nt_sheet, header_format, normal_format)
    sheet_filler(dnt_content, dnt_sheet, header_format, normal_format)
    sheet_filler(hexa_content, hexa_sheet, header_format, normal_format)
    sheet_filler(codon_content, codon_sheet, header_format, normal_format)
    sheet_filler(codon_pos_content, codon_pos_sheet, header_format, normal_format)
    sheet_filler(aa_content, aa_sheet, header_format, normal_format)
    sheet_filler(group_content, grp_sheet, header_format, normal_format)
    # sheet_filler(propensity_content, prop_sheet, header_format, normal_format)

    workbook.close()


def create_fusions_dic(list_dics):
    """
    :param list_dics: (list of dic float values)
    :return: a dictionary that contains all the keys and their associated value given by the list_dics
    """
    grp_dics = {}
    prop_dics = {}
    for i in range(6, 11):
        for key in list_dics[i].keys():
            grp_dics[key] = list_dics[i][key]
    """
    for i in range(8, 24):
        for key in list_dics[i].keys():
            prop_dics[key] = list_dics[i][key]
    """
    return grp_dics, prop_dics


def main(fasta_file, output, motif, exon_type):
    """
    Execute the core program

    :param fasta_file:(string) a fasta file
    :param output: (string) the path where the result will be created
    :param motif: (string) a nucleotide motif (it can be a regular expression)
    :param exon_type: (string) the name of the control set of exons
    """
    fasta_dics, all_list = create_dic(fasta_file)
    p_val_dics, ctrl_dics = calculate_all_p_value_dic(fasta_dics, exon_type, all_list)
    print("Creating hexa content")
    hexa_content = get_hexanucleotide_content(ctrl_dics[0], fasta_dics[0], p_val_dics[0],
                                              exon_type, motif, ctrl_dics[0]["all"], all_list[0])
    print("Creating dnt content")
    dnt_content = get_content(ctrl_dics[1], fasta_dics[1], p_val_dics[1],
                                             exon_type, "dnt", ctrl_dics[1]["all"], all_list[1])
    print("Creating nt content")
    nt_content = get_content(ctrl_dics[2], fasta_dics[2], p_val_dics[2],
                                             exon_type, "nt", ctrl_dics[2]["all"], all_list[2])
    print("Creating codon content")
    codon_content = get_codon_content(ctrl_dics[3], fasta_dics[3], p_val_dics[3],
                                      exon_type, ctrl_dics[3]["all"], all_list[3])
    print("Creating codon position content")
    codon_pos_content = get_content(ctrl_dics[4], fasta_dics[4], p_val_dics[4],
                                    exon_type, "codon_pos", ctrl_dics[4]["all"],
                                    all_list[4])
    print("Creating aa content")
    aa_content = get_content(ctrl_dics[5], fasta_dics[5], p_val_dics[5],
                             exon_type, "aa", ctrl_dics[5]["all"], all_list[5])
    print("Creating grp and prop content")
    grp_fasta_dics, prop_fasta_dics = create_fusions_dic(fasta_dics)
    grp_ctrl_dics, prop_ctrl_dics = create_fusions_dic(ctrl_dics)
    grp_pval_dics, prop_pval_dics = create_fusions_dic(p_val_dics)
    grp_content = get_content(grp_ctrl_dics, grp_fasta_dics, grp_pval_dics,
                              exon_type, "feature_groups", ctrl_dics[5]["all"], all_list[5])
    """
    prop_content = get_content(prop_ctrl_dics, prop_fasta_dics, prop_pval_dics, exon_type, "scaled_groups")
    """
    creating_report(hexa_content, dnt_content, nt_content, codon_content, codon_pos_content, aa_content, grp_content, output)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""Create an enrichment file that will compare
a fasta file to a desired control set of exons (taken from fasterDB)
    """,
                                     usage='%(prog)s --output output_folder --fasta a_fasta_file --ctrl ctrl_set'
                                           'of exons --motif nucleotid_motif')
    # Arguments for the parser

    required_group = parser.add_argument_group("required arguments")
    required_group.add_argument('--output', dest='output', help="An output folder", required=True)
    required_group.add_argument('--fasta', dest='fasta', help="Your fasta file", required=True)
    required_group.add_argument('--ctrl', dest='ctrl', help='your control either (ACE, CCE, ALL)', required=True)
    required_group.add_argument('--motif', dest='motif', help="a motif it can be a regexp motif", required=True)

    args = parser.parse_args()

    main(args.fasta, args.output, args.motif, args.ctrl)


if __name__ == "__main__":
    launcher()

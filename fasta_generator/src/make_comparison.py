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
from scipy.stats import mannwhitneyu
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
    :param dic: (a dictionary of float) freq of the word af interest
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


def calcul_dic_codon(dic, seq):
    """
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
    return  dic


def create_a_codon_dic(list_seq):
    """
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
    :param dic: (dictionary of int) the number of amino acids for a
    given set of sequences
    :param seq: (string) a nucleotides sequence
    :return: dic with the updated amino acids content
    """
    for i in range(0, len(seq)-2, 3):
        if codon2aminoAcid[seq[i:i+3]]!="":
            if codon2aminoAcid[seq[i:i+3]] not in dic.keys():
                dic[codon2aminoAcid[seq[i:i+3]]] = 1
            else:
                dic[codon2aminoAcid[seq[i:i+3]]] += 1
    return  dic


def create_an_aa_dic(list_seq):
    """
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
            if aa in  aa_dic.keys():
                res[key] += aa_dic[aa]
    return res


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
        if len(list_seq[i]) > 0:
            for j in range(len(list_seq[i])):
                val += dic[list_seq[i][j]]
            if correction:
                list_val.append(val / len(list_seq[i]))
            else:
                list_val.append(val)
    return list_val


def create_dic(fasta_file):
    """

    :param fasta_file:  (string) a fasta file
    :return: res_dic, [all_hexa, all_codon, all_aa]
         - res_dic : list of dictionary of
    """
    res_dic = []
    print("reading fasta")
    list_seq = fasta_reader(fasta_file)
    print("generation of hex dics")
    dic, all_hexa = create_an_hexanucleotid_dic(list_seq)
    res_dic.append(dic)
    print("generation of codon dics")
    dic, all_codon = create_a_codon_dic(list_seq)
    print("***")
    print(dic["TTT"])
    res_dic.append(dic)
    print(res_dic[-1])
    print("generation of aa dics")
    dic, all_aa = create_an_aa_dic(list_seq)
    res_dic.append(dic)
    print("generating group dics")
    res_dic.append(create_a_custom_dic(res_dic[2], schain2aa))
    res_dic.append(create_a_custom_dic(res_dic[2], hydro_info2aa))
    res_dic.append(create_a_custom_dic(res_dic[2], charge_info2aa))
    res_dic.append(create_a_custom_dic(res_dic[2], polarity_info2aa))
    res_dic.append(create_a_custom_dic(res_dic[2], misc2aa))
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
    return res_dic, [all_hexa, all_codon, all_aa]


def get_dic_pvalue(dic_fasta, dic_control, ctrl_all, fasta_all):
    p_val = {}
    for key in dic_fasta.keys():
        p1 = hypergeom.cdf(dic_fasta[key], ctrl_all, dic_control[key], fasta_all)
        p2 = 1 - hypergeom.cdf(dic_fasta[key], ctrl_all, dic_control[key], fasta_all)
        if p1 < p2:
            p_val[key] = p1
        else:
            p_val[key] = p2
    return p_val

def get_propensity_pvalue(dic_fasta, dic_control):
    p_val = {}
    for key in dic_fasta.keys():
        p_val[key] = mannwhitneyu(dic_fasta[key], dic_control[key])
    return  p_val

def calculate_random_dic():
    mod = __import__("CCE_dic")
    codon_number = mod.dc["all"]
    hexa_number = mod.d6["all"]
    d6_r = {}
    dc_r = {}
    da_r = {}
    for key in mod.d6.keys():
        if key != "all":
            d6_r[key] = int(math.pow(0.25, 6) * hexa_number)
    for key in mod.dc.keys():
        if key != "all":
            dc_r[key] = int(math.pow(0.25, 3) * codon_number)
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
    for key in dc_r.keys():
        count += dc_r[key]
    dc_r['all'] = count
    count = 0
    for key in da_r.keys():
        count += da_r[key]
    da_r['all'] = count
    return [d6_r, dc_r, da_r, sh, hy, ch, po, mi]



def calculate_all_p_value_dic(fasta_dics, exon_type, all_list):
    print("calculation of pvalues")
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, file_dir + "/control_dic/")
    p_val_dics = []
    if exon_type == "RD":
        ctrl_dics = calculate_random_dic()
    else:
        mod = __import__(exon_type + "_dic")
        ctrl_dics = [mod.d6, mod.dc, mod.da, mod.sh, mod.hy, mod.ch, mod.po, mod.mi]
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
        print(str(i) + "dics out of " + str(8))
        if i < 3:
            p_val_dics.append(get_dic_pvalue(fasta_dics[i], ctrl_dics[i], ctrl_dics[i]["all"], all_list[i]))
        elif i < 8:
            p_val_dics.append(get_dic_pvalue(fasta_dics[i], ctrl_dics[i], ctrl_dics[2]["all"], all_list[2]))
        """
        elif i >= 8:
            p_val_dics.append(get_propensity_pvalue(fasta_dics[i], ctrl_dics[i]))
        """
    return p_val_dics, ctrl_dics


def find_hexa_motif(motif, hexa_fasta):
    motif_dic = {}
    for hexant  in hexa_fasta.keys():
        res = re.findall(str(motif), str(hexant))
        if len(res) > 0:
            motif_dic[hexant] = "MOTIF"
        else:
            motif_dic[hexant] = "-"
    return motif_dic

def get_p_adjust_dic(pval_dic):
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
        content.append([key, motif_dic[key], float(hexa_fasta[key])/fasta_all, float(hexa_ctrl[key])/ctrl_all, pval_hexa[key], p_adj[key],
                        reg, reg_cor])
    return content

def get_codon_content(codon_ctrl, codon_fasta, pval_codon, exon_type, ctrl_all, fasta_all):
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
        content.append([key, codon2aminoAcid[key], float(codon_fasta[key]) / fasta_all, float(codon_ctrl[key]) / ctrl_all, pval_codon[key], p_adj[key],
                        reg, reg_cor])
    return content

def get_content(ctrl, fasta, pval, exon_type, type_ft, ctrl_all, fasta_all):
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
        content.append([key, float(fasta[key]) / fasta_all, float(ctrl[key]) / ctrl_all, pval[key], p_adj[key], reg, reg_cor])
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


def creating_report(hexa_content, codon_content, aa_content, group_content, output):
    workbook = xlsxwriter.Workbook(output + "enrichment_report.xlsx")
    header_format = workbook.add_format({'bg_color': '#00DCFF', 'align': 'center', 'valign': 'vcenter', 'border': True})
    normal_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': True})

    hexa_sheet = workbook.add_worksheet("haxnucleotide")
    codon_sheet = workbook.add_worksheet("codon")
    aa_sheet = workbook.add_worksheet("amino_acid")
    grp_sheet = workbook.add_worksheet("our_group")
    #prop_sheet = workbook.add_worksheet("propensity")

    sheet_filler(hexa_content, hexa_sheet, header_format, normal_format)
    sheet_filler(codon_content, codon_sheet, header_format, normal_format)
    sheet_filler(aa_content, aa_sheet, header_format, normal_format)
    sheet_filler(group_content, grp_sheet, header_format, normal_format)
    #sheet_filler(propensity_content, prop_sheet, header_format, normal_format)

    workbook.close()


def create_fusions_dic(list_dics):
    grp_dics = {}
    prop_dics = {}
    for i in range(3, 8):
        for key in list_dics[i].keys():
            grp_dics[key] = list_dics[i][key]
    """
    for i in range(8, 24):
        for key in list_dics[i].keys():
            prop_dics[key] = list_dics[i][key]
    """
    return grp_dics, prop_dics

def main(fasta_file, output, motif, exon_type):
    fasta_dics, all_list = create_dic(fasta_file)
    p_val_dics, ctrl_dics = calculate_all_p_value_dic(fasta_dics, exon_type, all_list)
    print("Creating hexa content")
    hexa_content = get_hexanucleotide_content(ctrl_dics[0], fasta_dics[0], p_val_dics[0], exon_type, motif, ctrl_dics[0]["all"], all_list[0])
    print("Creating codon content")
    codon_content = get_codon_content(ctrl_dics[1], fasta_dics[1], p_val_dics[1], exon_type, ctrl_dics[1]["all"], all_list[1])
    print("Creating aa content")
    aa_content = get_content(ctrl_dics[2], fasta_dics[2], p_val_dics[2], exon_type, "aa", ctrl_dics[2]["all"], all_list[2])
    print("Creating grp and prop content")
    grp_fasta_dics, prop_fasta_dics = create_fusions_dic(fasta_dics)
    grp_ctrl_dics, prop_ctrl_dics = create_fusions_dic(ctrl_dics)
    grp_pval_dics, prop_pval_dics = create_fusions_dic(p_val_dics)
    grp_content = get_content(grp_ctrl_dics, grp_fasta_dics, grp_pval_dics, exon_type, "feature_groups", ctrl_dics[2]["all"], all_list[2])
    """
    prop_content = get_content(prop_ctrl_dics, prop_fasta_dics, prop_pval_dics, exon_type, "scaled_groups")
    """
    creating_report(hexa_content, codon_content, aa_content, grp_content, output)


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

    required_group = parser.add_argument_group("required arguments")
    required_group.add_argument('--output', dest='output', help="An output folder", required=True)
    required_group.add_argument('--fasta', dest='fasta', help="Your fasta file", required=True)
    required_group.add_argument('--ctrl', dest='ctrl', help='your control either (ACE, CCE, ALL)', required=True)
    parser.add_argument('--motif', dest='motif', help="a motif it can be a regexp motif", required=True)

    args = parser.parse_args()

    main(args.fasta, args.output, args.motif, args.ctrl)



if __name__ == "__main__":
    launcher()

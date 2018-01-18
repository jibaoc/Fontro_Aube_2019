##############################################################################
#                                  Description                               #
##############################################################################
# This script allows  the creation of the query_summary.xlsx file and the
# creation of the enrichment report.xlsx file

##############################################################################
#                                 File                                       #
##############################################################################

import xlsxwriter
from dictionnary import *

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import numpy as np


group_nature2complete_name = { "NP": "Non polar", "NP-Alkyl": "Non polar alkyl",
                               "NP-aromatic": "Non polar and aromatic", "NPR": "Non polar restrained", "P": "Polar",
                               "PN": "Polar neutral", "PNC": "Polar uncharged", "PNC1": "Polar uncharged 1",
                               "PNC2" : "Polar uncharged 2", "PC": "Polar charged", "P-NC": "Polar negatively charged",
                               "P-PC": "Polar positively charged", "HC": "Hydrophobic chain", "H": "Hydrophobic",
                               "Aliphatic": "Aliphatic", "HS": "hydroxylated/sulfured", "Aromatic": "Aromatic"

}

##############################################################################
#                             Functions                                      #
##############################################################################

# --------------------------- query_summary.xlsx -----------------------------


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


def writing_query_results_file(outpath, input_content, exon_list):
    """
    Creates the file query_results.xlsx
    :param outpath: (string) the path where the query_results.xlsx file will be created
    :param input_content: (list of list of string) the content of the input sheet
    :param exon_list: (instance of ListExon) a list of exon
    """
    workbook = xlsxwriter.Workbook(outpath + "query_results.xlsx")
    # setting the formats
    header_format = workbook.add_format({'bg_color': '#00DCFF', 'align': 'center', 'valign': 'vcenter', 'border': True})
    normal_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': True})
    # getting the sheet content
    mapping_content = exon_list.get_content_mapping_sheet()
    sequence_content = exon_list.get_content_sequence_sheet()
    feature_content = exon_list.get_content_feature_sheet()
    # creating the sheets...
    input_sheet = workbook.add_worksheet("input")
    mapping_sheet = workbook.add_worksheet("mapping")
    sequence_sheet = workbook.add_worksheet("sequence")
    feature_sheet = workbook.add_worksheet("feature")
    # filling the sheets...
    sheet_filler(input_content, input_sheet, header_format, normal_format)
    sheet_filler(mapping_content, mapping_sheet, header_format, normal_format)
    sheet_filler(sequence_content, sequence_sheet, header_format, normal_format)
    sheet_filler(feature_content, feature_sheet, header_format, normal_format)
    workbook.close()

# --------------------------- enrichment report.xlsx -----------------------------


def check_regulation(user_frequency, ic_95, p_value, p_value_corrected):
    """
    :param user_frequency: (float) a frequency of a codon/amino acid/nature from the exon set given by the user
    :param ic_95: (list of 2 floats) an interval containing 95% of the frequency of the codon/amino acid/nature of
    interest
    :param p_value: (float) a p_value indicating if the frequency of the codon/amino acid/nature of interest is
    different from the frequency of the same codon/amino acid/nature in the control set
    :param p_value_corrected: (float) the corrected p_value (fdr method)
    :return: (2 strings) : regulation that can be equal to "=" if the pvalue is <0.05 or "+" or "-" if the pvalue is
    <0.05. The choice of the regulation if done thanks to the confidence interval that contains 95% of the frequency of
    the codon/amino acid/nature of interest 'ic_95'
    """
    regulation = " = "
    regulation_fdr = " = "
    if p_value <= 0.05:
        if user_frequency < ic_95[0]:
            regulation = "-"
        else:
            regulation = "+"
    if p_value_corrected <= 0.05:
        if user_frequency < ic_95[0]:
            regulation_fdr = "-"
        else:
            regulation_fdr = "+"
    return regulation, regulation_fdr


def calculate_ic_95(control_frequencies):
    """
    :param control_frequencies: frequencies of the control sets
    :return: the interval containing 95% of the frequencies for a given codon (there are as much as frequencies for a
    codon than the number of control sets)
    """
    ic_95 = dict()
    for codon_or_amino_acid_or_nature in control_frequencies:
        control_frequencies[codon_or_amino_acid_or_nature].sort()
        ic_95[codon_or_amino_acid_or_nature] = list()
        ic_95[codon_or_amino_acid_or_nature].append(control_frequencies[codon_or_amino_acid_or_nature][int(
            len(control_frequencies[codon_or_amino_acid_or_nature]) * 0.025)])
        ic_95[codon_or_amino_acid_or_nature].append(control_frequencies[codon_or_amino_acid_or_nature][int(
            len(control_frequencies[codon_or_amino_acid_or_nature]) * 0.975)])
    return ic_95


def get_content_codon_enrichment(control_frequencies, interest_frequencies, interest_frequencies_5p, interest_frequencies_3p, dic_p_val, set_number):
    """
    :param control_frequencies: (dictionary of float) a dictionary containing the codon frequencies of the control sets
    :param interest_frequencies: (dictionary of float) a dictionary frequency of each codon in the user set of exons
    :param dic_p_val: (dic of floats) a dictionary containing the p_values
    :param set_number: (int) the number of set to create
    :return: (list of list of strings) the content of the codon sheet
    """
    codon_list = ["TAG" , "TAA" , "TGA", "GCA" , "GCC" , "GCG" , "GCT" , "TGT" , "TGC" , "GAT" , "GAC" , "GAA" ,
                  "GAG" , "TTC" , "TTT" , "GGT" , "GGG" , "GGA" , "GGC" , "CAT" , "CAC" , "ATC" , "ATA" , "ATT" ,
                  "AAG" , "AAA" , "CTG" , "CTA" , "CTC" , "CTT" , "TTG" , "TTA" , "ATG" , "AAC" , "AAT" , "CCT" ,
                  "CCG" , "CCA" , "CCC" , "CAA" , "CAG" , "AGG" , "AGA" , "CGA" , "CGC" , "CGG" , "CGT" , "AGC" ,
                  "AGT" , "TCG" , "TCC" , "TCA" , "TCT" , "ACC" , "ACA" , "ACG" , "ACT" , "GTA" , "GTC" , "GTG" ,
                  "GTT" , "TGG" , "TAT" , "TAC"]
    dic_padjust = {}
    content = [["codon", "tRNA", "amino_acid", "frequencies_of_the_interest_set", "frequencies_interest_set_5p", "frequencies_interest_set_3p",
                "average_frequencies_of_the_"+str(set_number)+"_sets", "IC_95_of_the_"+str(set_number)+"_sets",
                "p_values_like", "FDR", "regulation_(p<=0.05)", "regulation(fdr<=0.05)", "codon_info"]]
    ic_95 = calculate_ic_95(control_frequencies)
    p_vals = list()
    for codon in codon_list:
        p_vals.append(dic_p_val[codon])
    rstats = importr('stats')
    p_adjust = rstats.p_adjust(FloatVector(p_vals), method="BH")
    i = 0
    for codon in codon_list:
        if codon2rare[codon] == "+":
            info_codon = "Frequent"
        elif codon2rare[codon] == "-":
            info_codon = "Rare"
        else:
            info_codon = ""
        regulation, regulation_fdr = check_regulation(interest_frequencies[codon], ic_95[codon], dic_p_val[codon],
                                                      p_adjust[i])
        content.append([str(codon), str(codon2anticodon[codon]), str(codon2aminoAcid[codon]),
                        str(interest_frequencies[codon]), str(interest_frequencies_5p[codon]), str(interest_frequencies_3p[codon]), str(np.mean(control_frequencies[codon])),
                        str(ic_95[codon]), str(dic_p_val[codon]), str(p_adjust[i]), str(regulation),
                        str(regulation_fdr), str(info_codon)])
        dic_padjust[codon] = p_adjust[i]
        i += 1
    return content, dic_padjust


def get_content_amino_acid_enrichment(control_frequencies, interest_frequencies, interest_frequencies_5p, interest_frequencies_3p, dic_p_val, set_number):
    """
    :param control_frequencies: (dictionary of floats) a dictionary containing the amino acid frequencies of the control
    sets
    :param interest_frequencies: (dictionary of floats) a dictionary frequency of each amino acid in the user set of
    exons
    :param dic_p_val: (dictionary of floats) a dictionary containing the p_values
    :param set_number: (int) the number of set to create
    :return: (list of list of strings) the content of the amino_acid sheet ! Each sublist correspond to a row in the
    amino_acid sheet of the enrichment_report.xlsx file
    """
    dic_padjust = {}
    content = [["amino_acid", "frequencies_of_the_interest_set", "frequencies_interest_set_5p", "frequencies_interest_set_3p",
                "average_frequencies_of_the_"+str(set_number)+"_sets",
                "IC_95_of_the_"+str(set_number)+"_sets", "p_values_like", "FDR", "regulation_(p<=0.05)",
                "regulation(fdr<=0.05)"]]
    ic_95 = calculate_ic_95(control_frequencies)
    p_vals = list()
    for amino_acid in sorted(dic_p_val.keys()):
        p_vals.append(dic_p_val[amino_acid])
    rstats = importr('stats')
    p_adjust = rstats.p_adjust(FloatVector(p_vals), method="BH")
    i = 0
    for amino_acid in sorted(dic_p_val.keys()):
        regulation, regulation_fdr = check_regulation(interest_frequencies[amino_acid], ic_95[amino_acid],
                                                      dic_p_val[amino_acid], p_adjust[i])
        content.append([str(amino_acid), str(interest_frequencies[amino_acid]), str(interest_frequencies_5p[amino_acid]), str(interest_frequencies_3p[amino_acid]),
                        str(np.mean(control_frequencies[amino_acid])), str(ic_95[amino_acid]),
                        str(dic_p_val[amino_acid]), str(p_adjust[i]), str(regulation), str(regulation_fdr)])
        dic_padjust[amino_acid] = p_adjust[i]
        i += 1
    return content, dic_padjust


def get_content_group_enrichment(control_frequencies, interest_frequencies, interest_frequencies_5p, interest_frequencies_3p, dic_p_val, set_number, name, list_key=None):
    """
    :param control_frequencies: (dictionary of floats) a dictionary containing the amino acid nature frequencies of the
    control sets
    :param interest_frequencies: (dictionary of floats) a dictionary frequency of each amino acid nature in the user set
    of exons
    :param dic_p_val: (dictionary of floats a dictionary containing the p_values
    :param set_number: (int) the number of set to create
    :param name: (string) : the name of the first column
    :param list_key: (list of string) list of key contained in control_frequencies, interest_frequencies, dic_p_val
    :return: (list of list of strings) the content of the nature sheet ! Each sublist correspond to a row in the
    nature sheet of the enrichment_report.xlsx file
    """
    dic_padjust = {}
    content = [[name, "frequencies_of_the_interest_set", "frequencies_interest_set_5p", "frequencies_interest_set_3p",
                "average_frequencies_of_the_"+str(set_number)+"_sets", "IC_95_of_the_"+str(set_number)+"_sets",
                "p_values_like", "FDR", "regulation_(p<=0.05)", "regulation(fdr<=0.05)"]]
    ic_95 = calculate_ic_95(control_frequencies)
    p_vals = list()
    for nature in dic_p_val.keys():
        p_vals.append(dic_p_val[nature])
    rstats = importr('stats')
    p_adjust = list(rstats.p_adjust(FloatVector(p_vals), method="BH"))
    padjust = {}
    i=0
    for nature in dic_p_val.keys():
        padjust[nature] = p_adjust[i]
        i += 1
    if list_key is None:
        for nature in dic_p_val.keys():
            regulation, regulation_fdr = check_regulation(interest_frequencies[nature], ic_95[nature],
                                                      dic_p_val[nature], padjust[nature])
            content.append([str(nature), str(interest_frequencies[nature]),  str(interest_frequencies_5p[nature]),  str(interest_frequencies_3p[nature]),
                        str(np.mean(control_frequencies[nature])), str(ic_95[nature]), str(dic_p_val[nature]),
                        str(padjust[nature]), str(regulation), str(regulation_fdr)])


    else:
        for nature in list_key:
            if nature == " ":
                content.append([" " for i in range(7)])
            else:
                regulation, regulation_fdr = check_regulation(interest_frequencies[nature], ic_95[nature], dic_p_val[nature], padjust[nature])
                content.append([str(nature), str(interest_frequencies[nature]), str(interest_frequencies_5p[nature]),  str(interest_frequencies_3p[nature]),
                        str(np.mean(control_frequencies[nature])), str(ic_95[nature]), str(dic_p_val[nature]),
                        str(padjust[nature]), str(regulation), str(regulation_fdr)])

    return content, dic_padjust

def create_iupac_dic(dic):
    """

    :param dic: (dictionary of float) must have the following keys : A, T, G, C
    :return: a dictionary with ambiguous iupac nt
    """
    res_iupac = {}
    res_iupac["Y"] = dic["C"] + dic["T"]
    res_iupac["R"] = dic["A"] + dic["G"]
    res_iupac["W"] = dic["T"] + dic["A"]
    res_iupac["S"] = dic["C"] + dic["G"]
    res_iupac["K"] = dic["T"] + dic["G"]
    res_iupac["M"] = dic["A"] + dic["C"]
    return res_iupac

def sorted_string(a_dic, nt_string=None):
    """
    :param a_dic: (dictionary of float) - the key of the dictionary are nucleotides and their associated value their
    frequencies in a sequence of interest
    :param nt_string: (string) sequence of nucleotides that corresponds to frequencies in "a_dic" dictionary
    :return: (list of tuple) ordered by the float value in the dictionary
    """
    res_str = ""
    res_prop = ""
    list_tuple = sorted(a_dic.items(), key=lambda x: x[1], reverse=True)
    for my_tuple in list_tuple:
        res_str += str(my_tuple[0]) + "(" + str(my_tuple[1]) + ")" + " - "
        if nt_string is not None:
            res_prop += str(my_tuple[0]) + "(" + str(round(float(my_tuple[1] * 100) / len(nt_string), 1)) + ")" + " - "

    res_str = res_str[:len(res_str)-3]
    if nt_string is not None:
        res_prop = res_prop[:len(res_prop) - 3]
        return res_str, res_prop
    return res_str

def get_group_nt_info(list_aa):
    """

    :param list_aa: (list of string) list of amino acid
    :return: the number of nt in all the codon coded by all the amino acid in the list, their proportion and their
    pondered proportion
    example : list_aa = K, W
    codon list = AAA, AAG, TGG:
    count_str = A(5) - G(3) - T(1) - C(0) - Y(1) - R(8) - S(3) - W(6) - K(4) - M(5) -D....
    count_prop = A(55.6) - G(33.3) - ....
    count_pond = A((5./6 + 0./3) * 100 /2 = 41.6) - G(41.6)
    """
    res = {"A":0, "T":0, "G":0, "C":0}
    ponderate = {"A":0., "T":0., "G":0., "C":0.}
    curstr = ""
    for aa in list_aa:
        cur_codon = amino_acid2codon[aa].replace(",","")
        for nt in ponderate.keys():
           ponderate[nt] += float(cur_codon.count(nt)) / len(cur_codon)
        curstr += cur_codon

    for nt in ponderate.keys():
        ponderate[nt] = round(ponderate[nt] * 100 / len(list_aa), 1)

    count_pond = sorted_string(ponderate)

    for key in res.keys():
        res[key] = curstr.count(key)

    count_str, count_prop = sorted_string(res, curstr)

    res_iupac = create_iupac_dic(res)
    pond_iupac = create_iupac_dic(ponderate)

    count_str_iupac, count_prop_iupac = sorted_string(res_iupac, curstr)
    count_pond_iupac = sorted_string(pond_iupac)
    count_str = count_str + " - " + count_str_iupac
    count_prop = count_prop + " - " + count_prop_iupac
    count_pond = str(count_pond) + " - " + str(count_pond_iupac)

    return count_str, count_prop, count_pond





def get_content_group_enrichment2(control_frequencies, interest_frequencies, interest_frequencies_5p, interest_frequencies_3p, dic_p_val, set_number, name, reg_dic):
    """
    :param control_frequencies: (dictionary of floats) a dictionary containing the amino acid nature frequencies of the
    control sets
    :param interest_frequencies: (dictionary of floats) a dictionary frequency of each amino acid nature in the user set
    of exons
    :param dic_p_val: (dictionary of floats a dictionary containing the p_values
    :param set_number: (int) the number of set to create
    :param name: (string) : the name of the first column
    :param reg_dic: (dictionary of list of string) : dictionary having keys corresponding to the group of interest and
    a list associated to those keys corresponding to the amino acid aggregated in those groups
    :return: (list of list of strings) the content of the nature sheet ! Each sublist correspond to a row in the
    nature sheet of the enrichment_report.xlsx file
    """
    dic_padjust = {}
    content = [[name, "frequencies_of_the_interest_set", "frequencies_interest_set_5p", "frequencies_interest_set_3p",
                "average_frequencies_of_the_"+str(set_number)+"_sets", "IC_95_of_the_"+str(set_number)+"_sets",
                "p_values_like", "FDR", "regulation_(p<=0.05)", "regulation(fdr<=0.05)", "nb_nt_group", "prop_nt_group", "ponderate_nt_group"]]
    ic_95 = calculate_ic_95(control_frequencies)
    p_vals = list()
    for nature in dic_p_val.keys():
        p_vals.append(dic_p_val[nature])
    rstats = importr('stats')
    p_adjust = rstats.p_adjust(FloatVector(p_vals), method="BH")
    i = 0
    for nature in dic_p_val.keys():
        info_count, info_prop, count_pond = get_group_nt_info(reg_dic[nature])
        regulation, regulation_fdr = check_regulation(interest_frequencies[nature], ic_95[nature],
                                                      dic_p_val[nature], p_adjust[i])
        content.append([str(nature), str(interest_frequencies[nature]), str(interest_frequencies_5p[nature]), str(interest_frequencies_3p[nature]),
                        str(np.mean(control_frequencies[nature])), str(ic_95[nature]), str(dic_p_val[nature]),
                        str(p_adjust[i]), str(regulation), str(regulation_fdr), str(info_count), str(info_prop), str(count_pond)])
        dic_padjust[nature] = p_adjust[i]
        i += 1

    return content, dic_padjust


def writing_enrichment_report_file(control_frequencies_codon, codon_frequencies, codon_frequencies_5p, codon_frequencies_3p, dic_p_val_codon,
                                           control_frequencies_aa, aa_frequencies, aa_frequencies_5p, aa_frequencies_3p, dic_p_val_aa,
                                           control_schain_frequencies, schain_frequency, schain_frequency_5p, schain_frequency_3p, dic_p_val_schain,
                                           control_hydro_frequencies, hydro_frequency, hydro_frequency_5p, hydro_frequency_3p, dic_p_val_hydro,
                                           control_charge_frequencies, charge_frequency, charge_frequency_5p, charge_frequency_3p, dic_p_val_charge,
                                           control_polarity_frequencies, polarity_frequency, polarity_frequency_5p, polarity_frequency_3p, dic_p_val_polarity,
                                           control_misc_frequencies, misc_frequency, misc_frequency_5p, misc_frequency_3p, dic_p_val_misc,
                                           control_chimical_frequencies, chimical_frequency, chimical_frequency_5p, chimical_frequency_3p, dic_p_val_chimical,
                                           control_structural_frequencies, structural_frequency, structural_frequency_5p, structural_frequency_3p, dic_p_val_structural,
                                           control_nucleic_acid_frequencies, nucleic_acid_frequency, nucleic_acid_frequency_5p, nucleic_acid_frequency_3p, dic_p_val_nt,
                                           control_ntp_frequencies, ntp_frequency, ntp_frequency_5p, ntp_frequency_3p, dic_p_val_ntp,
                                           control_dnt_frequencies, dnt_frequency, dnt_frequency_5p, dnt_frequency_3p, dic_p_val_dnt,
                                           control_hexa_frequencies, hexa_frequency, hexa_frequency_5p, hexa_frequency_3p,
                                           dic_p_val_hexa,
                                           output_folder, set_number):
    """
    :param control_frequencies_codon:  (dictionary of floats) a dictionary containing the codon frequencies of the
    control sets
    :param interest_frequencies_codon:  (dictionary of floats) a dictionary frequency of each codon in the user set of
    exons
    :param dic_p_val_codon: (dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each codon
    :param control_frequencies_aa: (dictionary of floats) a dictionary containing the amino acid frequencies of the
    control
    sets
    :param interest_frequencies_aa: (dictionary of floats) a dictionary frequency of each amino acid in the user set of
    exons
    :param dic_p_val_aa: (dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each amino acid
    :param control_schain_frequencies: (dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their side chain properties) of the control sets
    :param schain_frequency:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their side chain properties) in the user set of exons
    :param dic_p_val_schain:(dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each sets of amino acid (grouped by their side chain properties)
    :param control_hydro_frequencies:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their hydrophilic/phobic properties) of the control sets
    :param hydro_frequency:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their hydrophilic/phobic properties) in the user set of exons
    :param dic_p_val_hydro:(dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each sets of amino acid (grouped by their hydrophilic/phobic properties)
    :param control_charge_frequencies:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their charge properties) of the control sets
    :param charge_frequency:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their charge properties) in the user set of exons
    :param dic_p_val_charge:(dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each sets of amino acid (grouped by their charge properties)
    :param control_polarity_frequencies:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their polarity properties) of the control sets
    :param polarity_frequency: (dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their polarity properties) in the user set of exons
    :param dic_p_val_polarity:(dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each sets of amino acid (grouped by their polarity properties)
    :param control_misc_frequencies:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by miscellaneous properties) of the control sets
    :param misc_frequency:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by miscellaneous properties) in the user set of exons
    :param dic_p_val_misc:(dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each sets of amino acid (grouped by miscellaneous properties)
    :param control_chimical_frequencies:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their chemical properties) of the control sets
    :param chimical_frequency:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their chemical properties) in the user set of exons
    :param dic_p_val_chimical: (dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each sets of amino acid (grouped by their chemical properties)
    :param control_structural_frequencies: (dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their structural properties) of the control sets
    :param structural_frequency:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their chemical properties) in the user set of exons
    :param dic_p_val_structural:(dictionary of floats) a dictionary containing the frequencies of different
    sets of amino acid (grouped by their structural properties) in the user set of exons
    :param control_nucleic_acid_frequencies:(dictionary of floats) a dictionary containing the  nucleotides frequencies
     of exons of the control sets
    :param nucleic_acid_frequency:(dictionary of floats) a dictionary containing the  nucleotides frequencies
     of exons of the user set of exons
    :param dic_p_val_nt:  (dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each nucleotides of the user set of exons
    :param output_folder: (string) the folder where the report enrichment_report.xlsx will be created
    :param set_number: (int) the number of control set that have been created
    :return:
    """
    workbook = xlsxwriter.Workbook(output_folder + "enrichment_report.xlsx")
    # setting the formats
    header_format = workbook.add_format({'bg_color': '#00DCFF', 'align': 'center', 'valign': 'vcenter', 'border': True})
    normal_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': True})
    # getting the sheet content
    codon_content, dic_padjust_codon = get_content_codon_enrichment(control_frequencies_codon,
                                                                    codon_frequencies, codon_frequencies_5p, codon_frequencies_3p, dic_p_val_codon,
                                                                    set_number)
    aa_content, dic_padjust_aa = get_content_amino_acid_enrichment(control_frequencies_aa, aa_frequencies, aa_frequencies_5p, aa_frequencies_3p,
                                                                   dic_p_val_aa, set_number)
    schain_content, dic_padjust_schain = get_content_group_enrichment2(control_schain_frequencies, schain_frequency, schain_frequency_5p, schain_frequency_3p, dic_p_val_schain,
                                                                      set_number, "side_chain_nature", schain2aa)
    hydro_content, dic_padjust_hydro = get_content_group_enrichment2(control_hydro_frequencies, hydro_frequency, hydro_frequency_5p, hydro_frequency_3p, dic_p_val_hydro,
                                                                      set_number, "hydric_relation", hydro_info2aa)
    charge_content, dic_padjust_charge = get_content_group_enrichment2(control_charge_frequencies, charge_frequency, charge_frequency_5p, charge_frequency_3p, dic_p_val_charge,
                                                                      set_number, "charge_group", charge_info2aa)
    polarity_content, dic_padjust_polarity = get_content_group_enrichment2(control_polarity_frequencies, polarity_frequency, polarity_frequency_5p, polarity_frequency_3p, dic_p_val_polarity,
                                                                      set_number, "polarity_group", polarity_info2aa)
    misc_content, dic_padjust_misc = get_content_group_enrichment2(control_misc_frequencies, misc_frequency, misc_frequency_5p, misc_frequency_3p, dic_p_val_misc,
                                                                      set_number, "misc_group", misc2aa)

    chimical_content, dic_padjust_chimical = get_content_group_enrichment(control_chimical_frequencies, chimical_frequency, chimical_frequency_5p, chimical_frequency_3p, dic_p_val_chimical,
                                                                      set_number, "chimical_group")
    structural_content, dic_padjust_structural = get_content_group_enrichment(control_structural_frequencies, structural_frequency, structural_frequency_5p, structural_frequency_3p, dic_p_val_structural,
                                                                      set_number, "structural_group")

    nt_list = ["A", "G", "C", "T", " ", "Y", "R", " ", "S", "W", " ", "K", "M", " ", "D", "C", " ", "V", "T", " ", "H", "G", " ",
               "B", "A"]
    nt_content, dic_padjust_nt = get_content_group_enrichment(control_nucleic_acid_frequencies,
                                                                      nucleic_acid_frequency,
                                                                      nucleic_acid_frequency_5p,
                                                                      nucleic_acid_frequency_3p,
                                                                      dic_p_val_nt,
                                                                      set_number, "nt_info", nt_list)

    ntp_list = ["A1", "A1n2", "A2", "A3", "C1", "C1n2", "C2", "C3", "G1", "G1n2", "G2", "G3", "T1", "T1n2", "T2", "T3",
                "Y1", "Y2", "Y3", "R1", "R2", "R3", "S1", "S2", "S3", "W1", "W2", "W3", "K1", "K2", "K3", "M1", "M2", "M3",
                "D1", "D2", "D3", "V1", "V2", "V3", "H1", "H2", "H3", "B1", "B2", "B3"]
    ntp_content, dic_padjust_ntp = get_content_group_enrichment(control_ntp_frequencies,
                                                              ntp_frequency,
                                                              ntp_frequency_5p,
                                                              ntp_frequency_3p,
                                                              dic_p_val_ntp,
                                                              set_number, "ntp_info", ntp_list)

    iu = {"Y": ["C", "T"], "R": ["A", "G"], "W": ["A", "T"], "S": ["G", "C"], "K": ["T", "G"],
          "M": ["C", "A"]}
    dnt_list = []
    for k1 in iu.keys():
        for k2 in iu.keys():
            dnt_list.append(k1 + k2)
            dnt_list.append(iu[k1][0] + iu[k2][0])
            dnt_list.append(iu[k1][0] + iu[k2][1])
            dnt_list.append(iu[k1][1] + iu[k2][0])
            dnt_list.append(iu[k1][1] + iu[k2][1])
            dnt_list.append(" ")



    dnt_content, dic_padjust_dnt = get_content_group_enrichment(control_dnt_frequencies,
                                                                      dnt_frequency,
                                                                dnt_frequency_5p,
                                                                dnt_frequency_3p,
                                                                      dic_p_val_dnt,
                                                                      set_number, "dnt_info", dnt_list)

    hexa_content, dic_padjust_hexa = get_content_group_enrichment(control_hexa_frequencies,
                                                                              hexa_frequency,
                                                                              hexa_frequency_5p,
                                                                              hexa_frequency_3p,
                                                                              dic_p_val_hexa,
                                                                              set_number, "structural_group")
    # creating the sheets...
    codon_sheet = workbook.add_worksheet("codon")
    aa_sheet = workbook.add_worksheet("amino_acid")
    schain_sheet = workbook.add_worksheet("side_chain")
    hydro_sheet = workbook.add_worksheet("hydric_relation")
    charge_sheet = workbook.add_worksheet("charge_info")
    polarity_sheet = workbook.add_worksheet("polarity_info")
    misc_sheet = workbook.add_worksheet("misc_info")
    chimical_sheet = workbook.add_worksheet("chimical_info")
    structural_sheet = workbook.add_worksheet("structural_info")
    nt_sheet = workbook.add_worksheet("nt_info")
    ntp_sheet = workbook.add_worksheet("nt_pos_info")
    dnt_sheet = workbook.add_worksheet("dnt_info")
    hexa_sheet = workbook.add_worksheet("hexa_info")

    # filling the sheets...
    sheet_filler(codon_content, codon_sheet, header_format, normal_format)
    sheet_filler(aa_content, aa_sheet, header_format, normal_format)
    sheet_filler(schain_content, schain_sheet, header_format, normal_format)
    sheet_filler(hydro_content, hydro_sheet, header_format, normal_format)
    sheet_filler(charge_content, charge_sheet, header_format, normal_format)
    sheet_filler(polarity_content, polarity_sheet, header_format, normal_format)
    sheet_filler(misc_content, misc_sheet, header_format, normal_format)
    sheet_filler(chimical_content, chimical_sheet, header_format, normal_format)
    sheet_filler(structural_content, structural_sheet, header_format, normal_format)
    sheet_filler(nt_content, nt_sheet, header_format, normal_format)
    sheet_filler(ntp_content, ntp_sheet, header_format, normal_format)
    sheet_filler(dnt_content, dnt_sheet, header_format, normal_format)
    sheet_filler(hexa_content, hexa_sheet, header_format, normal_format)

    workbook.close()
    return dic_padjust_codon, dic_padjust_aa

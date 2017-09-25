##############################################################################
#                                  Description                               #
##############################################################################
# This script allows  the creation of the query_summary.xlsx file and the
# creation of the enrichment report.xlsx file

##############################################################################
#                                 File                                       #
##############################################################################

import xlsxwriter
from dictionnary import codon2rare
from dictionnary import codon2anticodon
from dictionnary import codon2aminoAcid
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


def get_content_codon_enrichment(control_frequencies, interest_frequencies, dic_p_val, set_number):
    """
    :param control_frequencies: (dictionary of float) a dictionary containing the codon frequencies of the control sets
    :param interest_frequencies: (dictionary of float) a dictionary frequency of each codon in the user set of exons
    :param dic_p_val: (dic of floats) a dictionary containing the p_values
    :param set_number: (int) the number of set to create
    :return: (list of list of strings) the content of the codon sheet
    """
    dic_padjust = {}
    content = [["codon", "tRNA", "amino_acid", "frequencies_of_the_interest_set",
                "average_frequencies_of_the_"+str(set_number)+"_sets", "IC_95_of_the_"+str(set_number)+"_sets",
                "p_values_like", "FDR", "regulation_(p<=0.05)", "regulation(fdr<=0.05)", "codon_info"]]
    ic_95 = calculate_ic_95(control_frequencies)
    p_vals = list()
    for codon in dic_p_val.keys():
        p_vals.append(dic_p_val[codon])
    rstats = importr('stats')
    p_adjust = rstats.p_adjust(FloatVector(p_vals), method="BH")
    i = 0
    for codon in dic_p_val.keys():
        if codon2rare[codon] == "+":
            info_codon = "Frequent"
        elif codon2rare[codon] == "-":
            info_codon = "Rare"
        else:
            info_codon = ""
        regulation, regulation_fdr = check_regulation(interest_frequencies[codon], ic_95[codon], dic_p_val[codon],
                                                      p_adjust[i])
        content.append([str(codon), str(codon2anticodon[codon]), str(codon2aminoAcid[codon]),
                        str(interest_frequencies[codon]), str(np.mean(control_frequencies[codon])),
                        str(ic_95[codon]), str(dic_p_val[codon]), str(p_adjust[i]), str(regulation),
                        str(regulation_fdr), str(info_codon)])
        dic_padjust[codon] = p_adjust[i]
        i += 1
    return content, dic_padjust


def get_content_amino_acid_enrichment(control_frequencies, interest_frequencies, dic_p_val, set_number):
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
    content = [["amino_acid", "frequencies_of_the_interest_set", "average_frequencies_of_the_"+str(set_number)+"_sets",
                "IC_95_of_the_"+str(set_number)+"_sets", "p_values_like", "FDR", "regulation_(p<=0.05)",
                "regulation(fdr<=0.05)"]]
    ic_95 = calculate_ic_95(control_frequencies)
    p_vals = list()
    for amino_acid in dic_p_val.keys():
        p_vals.append(dic_p_val[amino_acid])
    rstats = importr('stats')
    p_adjust = rstats.p_adjust(FloatVector(p_vals), method="BH")
    i = 0
    for amino_acid in dic_p_val.keys():
        regulation, regulation_fdr = check_regulation(interest_frequencies[amino_acid], ic_95[amino_acid],
                                                      dic_p_val[amino_acid], p_adjust[i])
        content.append([str(amino_acid), str(interest_frequencies[amino_acid]),
                        str(np.mean(control_frequencies[amino_acid])), str(ic_95[amino_acid]),
                        str(dic_p_val[amino_acid]), str(p_adjust[i]), str(regulation), str(regulation_fdr)])
        dic_padjust[amino_acid] = p_adjust[i]
        i += 1
    return content, dic_padjust


def get_content_nature_enrichment(control_frequencies, interest_frequencies, dic_p_val, set_number):
    """
    :param control_frequencies: (dictionary of floats) a dictionary containing the amino acid nature frequencies of the
    control sets
    :param interest_frequencies: (dictionary of floats) a dictionary frequency of each amino acid nature in the user set
    of exons
    :param dic_p_val: (dictionary of floats a dictionary containing the p_values
    :param set_number: (int) the number of set to create
    :return: (list of list of strings) the content of the nature sheet ! Each sublist correspond to a row in the
    nature sheet of the enrichment_report.xlsx file
    """
    content = [["amino_acid_nature", "frequencies_of_the_interest_set",
                "average_frequencies_of_the_"+str(set_number)+"_sets", "IC_95_of_the_"+str(set_number)+"_sets",
                "p_values_like", "FDR", "regulation_(p<=0.05)", "regulation(fdr<=0.05)"]]
    ic_95 = calculate_ic_95(control_frequencies)
    p_vals = list()
    for nature in dic_p_val.keys():
        p_vals.append(dic_p_val[nature])
    rstats = importr('stats')
    p_adjust = rstats.p_adjust(FloatVector(p_vals), method="BH")
    i = 0
    for nature in dic_p_val.keys():
        regulation, regulation_fdr = check_regulation(interest_frequencies[nature], ic_95[nature],
                                                      dic_p_val[nature], p_adjust[i])
        content.append([str(group_nature2complete_name[nature]), str(interest_frequencies[nature]),
                        str(np.mean(control_frequencies[nature])), str(ic_95[nature]), str(dic_p_val[nature]),
                        str(p_adjust[i]), str(regulation), str(regulation_fdr)])
        i += 1
    return content


def get_content_metabolism_enrichment(control_frequencies, interest_frequencies, dic_p_val, set_number):
    """
    :param control_frequencies: (dictionary of floats) a dictionary containing the frequencies of each amino acid that
    come from different
    metabolism (TCA_cycle, Glycolyse or pentoses)
    :param interest_frequencies: (dictionary of floats) a dictionary frequency of each "metabolism" where the amino acid
     come from, in the user set of exons
    :param dic_p_val: (dictionary of floats) a dictionary containing the p_values
    :param set_number: (int) the number of set to create
    :return: (list of list of strings) the content of the metabolism sheet ! Each sublist correspond to a row in the
    metabolism sheet of the enrichment_report.xlsx file
    """
    content = []
    ic_95 = calculate_ic_95(control_frequencies)
    content.append(["amino_acid_metabolism", "frequencies_of_the_interest_set",
                    "average_frequencies_of_the_"+str(set_number)+"_sets", "IC_95_of_the_"+str(set_number)+"_sets",
                    "p_values_like", "FDR", "regulation_(p<=0.05)", "regulation(fdr<=0.05)"])
    p_vals = list()
    for metabolism in dic_p_val.keys():
        p_vals.append(dic_p_val[metabolism])
    rstats = importr('stats')
    p_adjust = rstats.p_adjust(FloatVector(p_vals), method="BH")
    i = 0
    for metabolism in dic_p_val.keys():
        regulation, regulation_fdr = check_regulation(interest_frequencies[metabolism], ic_95[metabolism],
                                                      dic_p_val[metabolism], p_adjust[i])
        content.append([str(metabolism), str(interest_frequencies[metabolism]),
                        str(np.mean(control_frequencies[metabolism])), str(ic_95[metabolism]),
                        str(dic_p_val[metabolism]), str(p_adjust[i]), str(regulation), str(regulation_fdr)])
        i += 1
    return content


def writing_enrichment_report_file(control_frequencies_codon, interest_frequencies_codon, dic_p_val_codon,
                                   control_frequencies_aa, interest_frequencies_aa, dic_p_val_aa,
                                   control_frequencies_nature, interest_frequencies_nature, dic_p_val_nature,
                                   control_frequencies_metabolism,
                                   interest_frequencies_metabolism, dic_p_val_metabolism, outpath, set_number):
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
    :param control_frequencies_nature: (dictionary of floats) a dictionary containing the amino_acid nature frequencies
    of the control sets
    :param interest_frequencies_nature: (dictionary of floats) a dictionary frequency of each amino acid nature in the
    user set of exons
    :param dic_p_val_nature: (dictionary of floats) a dictionary containing the p_values indicating the significance of
    the enrichment for each amino acid nature
    :param control_frequencies_metabolism: (dictionary of floats) a dictionary frequency of each amino acid metabolism
    in the user set of exons : 'Metabolism' here corresponds to the origin metabolism of
    amino acids in the user set of exons
    :param interest_frequencies_metabolism: (dictionary of floats) a dictionary frequency of each amino acid metabolism
    in the user set of exons : 'Metabolism' here corresponds to the origin metabolism of
    amino acids in the user set of exons
    :param dic_p_val_metabolism: (dictionary of floats) a dictionary containing the p_values indicating the significance
     of the enrichment for each amino acid metabolism.'Metabolism' here corresponds to the origin metabolism of
    amino acids in the user set of exons
    :param outpath: (string) the folder where the report enrichment_report.xlsx will be created
    :param set_number: (int) the number of control set that have been created
    :return:
    """
    workbook = xlsxwriter.Workbook(outpath + "enrichment_report.xlsx")
    # setting the formats
    header_format = workbook.add_format({'bg_color': '#00DCFF', 'align': 'center', 'valign': 'vcenter', 'border': True})
    normal_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': True})
    # getting the sheet content
    codon_content, dic_padjust_codon = get_content_codon_enrichment(control_frequencies_codon,
                                                                    interest_frequencies_codon, dic_p_val_codon,
                                                                    set_number)
    aa_content, dic_padjust_aa = get_content_amino_acid_enrichment(control_frequencies_aa, interest_frequencies_aa,
                                                                   dic_p_val_aa, set_number)
    nature_content = get_content_nature_enrichment(control_frequencies_nature, interest_frequencies_nature,
                                                   dic_p_val_nature, set_number)
    metabolism_content = get_content_metabolism_enrichment(control_frequencies_metabolism,
                                                           interest_frequencies_metabolism, dic_p_val_metabolism,
                                                           set_number)
    # creating the sheets...
    codon_sheet = workbook.add_worksheet("codon")
    aa_sheet = workbook.add_worksheet("amino_acid")
    nature_sheet = workbook.add_worksheet("nature")
    metabolism_sheet = workbook.add_worksheet("metabolism")

    # filling the sheets...
    sheet_filler(codon_content, codon_sheet, header_format, normal_format)
    sheet_filler(aa_content, aa_sheet, header_format, normal_format)
    sheet_filler(nature_content, nature_sheet, header_format, normal_format)
    sheet_filler(metabolism_content, metabolism_sheet, header_format, normal_format)
    workbook.close()
    return dic_padjust_codon, dic_padjust_aa

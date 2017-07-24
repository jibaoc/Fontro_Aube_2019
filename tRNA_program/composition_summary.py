#######################################
#               Script
#######################################
# This script is a part of the main script : main_project.py
# it contains functions that allow the construction of of the file composition_summary.xlsx

# the meta-exon method :
#  It consists in concatenating all the exon in a set of exon. the lists in the dictionaries obtained with this method
#  have at the first position the number of time an element is seen in this meta-exon and in the second position, the
#  frequency of this element in this meta-exon : example :
# lets say we have 5 exons :  Exon 1 : length 17 codon, 1 codon codes for A, Exon 2 : length 17 codon, 3 codons code 'A'
# Exon 3 : length 17 codon, 5 codons code for A, Exon 4 : length 5 codon, 1 codon code for A,
# Exon 5 : length 20 codon, 2 codon code for A
# count : 1 + 3 +5 + 1 + 2
# frequency : (1 + 3 +5 + 1 + 2) / (17 +17 +17 + 5 + 20)


# The weighted method :
#  The lists in the dictionaries obtained with this method
#  have at the first position the number of time an element is seen in a set of exon of interest. It is the same value
#  obtained with the meta-exon method unless an exon is tiny (below 17 codon), in that case the number of occurence seen
#  in this exon will be penalized... example for the amino acid A
# lets say we have 5 exons :  Exon 1 : length 17 codon, 1 codon codes for A, Exon 2 : length 17 codon, 3 codons code 'A'
# Exon 3 : length 17 codon, 5 codons code for A, Exon 4 : length 5 codon, 1 codon code for A,
# Exon 5 : length 20 codon, 2 codon code for A
# the number of occurence is : 1 + 3 + 5 + 1 * 5./17 (tiny exon) + 2
# the second value of the lists of dictionaries obtained with this method are the frequency of an element (codon,
# anticodon of the set). Again for 'A' with the previous set of exon, its frequency is :
# (1./17 + 3./17  + 5./17 + 1./5 * 5./17 + 2./20) * 100 / (1 + 1 + 1 + 5./17 + 1)

#######################################
#               Imports
#######################################

from dictionnary import amino_acid2codon_list, codon2rareness
import string
import xlsxwriter


#######################################
#              Functions
#######################################


def get_content_first_sheet(dic_up_codon, dic_up_aa, dic_down_codon, dic_down_aa, dic_up_down_codon,
                            dic_up_down_aa, dic_ace_codon, dic_ace_aa, dic_cce_codon, dic_cce_aa):
    """
    This function allow to create the content of the sheets names "meta_exon_recap" and weighted_exon_recap" in the
    composition_summary.xlsx file from many dictionaries. Those previously mentioned dictionaries are described below :
    :param dic_up_codon: (dictionary of list of 2 floats) : gives for each codon of the exon set up-regulated in a given
    condition, its number of occurrence and its frequency calculated with the meta-exon method or with the weighted
    method
    :param dic_up_aa:(dictionary of list of 2 floats) : gives for each amino acid of the exon set up-regulated in a
    given condition, its number of occurrence and its frequency calculated with the meta-exon method or with the
    weighted method
    :param dic_down_codon:(dictionary of list of 2 floats) : gives for each codon of the exon set down-regulated in a
    given condition, its number of occurrence and its frequency calculated with the meta-exon method or with the
    weighted method
    :param dic_down_aa: (dictionary of list of 2 floats) : gives for each amino acid of the exon set down-regulated in a
    given condition, its number of occurrence and its frequency calculated with the meta-exon method or with the
    weighted method
    :param dic_up_down_codon:(dictionary of list of 2 floats) : gives for each codon of the exon set down-regulated or
    up-regulated in a given condition, its number of occurrence and its frequency calculated with the meta-exon method
    or with the weighted method
    :param dic_up_down_aa:(dictionary of list of 2 floats) : gives for each amino acid of the exon set down-regulated or
    up-regulated in a given condition, its number of occurrence and its frequency calculated with the meta-exon method
    or with the weighted method
    :param dic_ace_codon:(dictionary of list of 2 floats) : gives for each codon of the ACE control exon set its
    number of occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :param dic_ace_aa:(dictionary of list of 2 floats) : gives for each amino acid of the ACE control exon set its
    number of occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :param dic_cce_codon:(dictionary of list of 2 floats) : gives for each codon of the CCE control exon set its
    number of occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :param dic_cce_aa:(dictionary of list of 2 floats) : gives for each amino acid of the CCE control exon set its
    number of occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :return:  1 - main_content :  a list of list corresponding to the content of the first or third sheet of the
    composition_summary file. Each sublist within the main list corresponds to a line of the 1st or the 3rd sheet of
    the composition_summary file.
              2 - merging_value : list of string and int. It corresponds to the cells we want to merge in the
              composition_summary file
              3 - merging_content : list of strings that correspond to the content of the merged cells
    """
    # definition of some variable
    header = None
    merging_value = list()
    merging_content = None
    main_content = list()
    # definition of the header, merging values and merging_content
    if dic_down_codon is not None and dic_up_down_codon is not None:
        header = ["amino_acid", "up_exon", '  ', "down_exon", ' ', "up_down_exon", '', "ACE(up)", '', "CCE(up)",
                  '', "up_vs_down", "up_vs_ACE", "down_vs_ACE", "up_down_vs_ACE", "up_vs_CCE", "down_vs_CCE", "up_down_vs_CCE",
                  "codon", "up_exon_c", '', "down_exon_c", '', "up_down_exon_c", '', "ACEc(up)", '', "CCEc(up)", '',
                  "up_vs_down", "up_vs_ACE", "down_vs_ACE", "up_down_vs_ACE", "up_vs_CCE", "down_vs_CCE", "up_down_vs_CCE"]
        merging_value = ["B1:C1", "D1:E1", "F1:G1", "H1:I1", "J1:K1", "T1:U1", "V1:W1", "X1:Y1", "Z1:AA1",
                         "AB1:AC1"]
        merging_content = [header[1], header[3], header[5], header[7], header[9], header[19],
                           header[21], header[23], header[25], header[27]]
    elif dic_down_codon is not None and dic_up_down_codon is None:
        header = ["amino_acid", "up_exon", '  ', "down_exon", ' ', "ACE(up)", '', "CCE(up)", '',
                  "up_vs_down", "up_vs_ACE", "down_vs_ACE", "up_vs_CCE", "down_vs_CCE", "codon",
                  "up_exon", '', "down_exon", '', "ACE(up)", '', "CCE(up)", '', "up_vs_down", "up_vs_ACE",
                  "down_vs_ACE", "up_vs_CCE", "down_vs_CCE"]
        merging_content = [header[1], header[3], header[5], header[7], header[15], header[17], header[19], header[21]]
        merging_value = ["B1:C1", "D1:E1", "F1:G1", "H1:I1", "P1:Q1", "R1:S1", "T1:U1", "V1:W1"]
    elif dic_down_codon is None and dic_up_down_codon is None:
        header = ["amino_acid", "input_exon", '  ', "ACE(up)", '', "CCE(up)", '',
                  "input_vs_ACE", "input_vs_CCE", "codon",
                  "input_exon", '', "ACE(up)", '', "CCE(up)", '', "input_vs_ACE", "up_vs_CCE"]
        merging_content = [header[1], header[3], header[5], header[10], header[12], header[14]]
        merging_value = ["B1:C1", "D1:E1", "F1:G1", "K1:L1", "M1:N1", "O1:P1"]
    main_content.append(header)
    counter_tot = 1
    # creating the content of the sheet
    for aa in amino_acid2codon_list.keys():
        for i in range(len(amino_acid2codon_list[aa])):
            counter_tot += 1
            codon = amino_acid2codon_list[aa][i]
            my_row = None
            if i == 0:
                merging_value.append(counter_tot)
            if dic_down_codon is not None and dic_up_down_codon is not None:
                my_row = [aa, dic_up_aa[aa][0], dic_up_aa[aa][1], dic_down_aa[aa][0], dic_down_aa[aa][1],
                          dic_up_down_aa[aa][0], dic_up_down_aa[aa][1], dic_ace_aa[aa][0], dic_ace_aa[aa][1],
                          dic_cce_aa[aa][0], dic_cce_aa[aa][1], dic_up_aa[aa][1] - dic_down_aa[aa][1],
                          dic_up_aa[aa][1] - dic_ace_aa[aa][1],
                          dic_down_aa[aa][1] - dic_ace_aa[aa][1], dic_up_down_aa[aa][1] - dic_ace_aa[aa][1],
                          dic_up_aa[aa][1] - dic_cce_aa[aa][1], dic_down_aa[aa][1] - dic_cce_aa[aa][1],
                          dic_up_down_aa[aa][1] - dic_cce_aa[aa][1], codon +  codon2rareness[codon], dic_up_codon[codon][0],
                          dic_up_codon[codon][1], dic_down_codon[codon][0], dic_down_codon[codon][1],
                          dic_up_down_codon[codon][0], dic_up_down_codon[codon][1], dic_ace_codon[codon][0],
                          dic_ace_codon[codon][1], dic_cce_codon[codon][0], dic_cce_codon[codon][1],
                          dic_up_codon[codon][1] - dic_down_codon[codon][1],
                          dic_up_codon[codon][1] - dic_ace_codon[codon][1],
                          dic_down_codon[codon][1] - dic_ace_codon[codon][1],
                          dic_up_down_codon[codon][1] - dic_ace_codon[codon][1],
                          dic_up_codon[codon][1] - dic_cce_codon[codon][1],
                          dic_down_codon[codon][1] - dic_cce_codon[codon][1],
                          dic_up_down_codon[codon][1] - dic_cce_codon[codon][1]]
            elif dic_down_codon is not None and dic_up_down_codon is None:
                my_row = [aa, dic_up_aa[aa][0], dic_up_aa[aa][1], dic_down_aa[aa][0], dic_down_aa[aa][1],
                          dic_ace_aa[aa][0], dic_ace_aa[aa][1], dic_cce_aa[aa][0],
                          dic_cce_aa[aa][1], dic_up_aa[aa][1] - dic_ace_aa[aa][1],
                          dic_up_aa[aa][1] - dic_ace_aa[aa][1],
                          dic_down_aa[aa][1] - dic_ace_aa[aa][1],
                          dic_up_aa[aa][1] - dic_cce_aa[aa][1], dic_down_aa[aa][1] - dic_cce_aa[aa][1],
                          codon + codon2rareness[codon], dic_up_codon[codon][0],
                          dic_up_codon[codon][1], dic_down_codon[codon][0], dic_down_codon[codon][1],
                          dic_ace_codon[codon][0], dic_ace_codon[codon][1], dic_cce_codon[codon][0],
                          dic_cce_codon[codon][1], dic_up_codon[codon][1] - dic_down_codon[codon][1],
                          dic_up_codon[codon][1] - dic_ace_codon[codon][1],
                          dic_down_codon[codon][1] - dic_ace_codon[codon][1],
                          dic_up_codon[codon][1] - dic_cce_codon[codon][1],
                          dic_down_codon[codon][1] - dic_cce_codon[codon][1]]
            elif dic_down_codon is None and dic_up_down_codon is None:
                my_row = [aa, dic_up_aa[aa][0], dic_up_aa[aa][1],
                          dic_ace_aa[aa][0], dic_ace_aa[aa][1], dic_cce_aa[aa][0],
                          dic_cce_aa[aa][1], dic_up_aa[aa][1] - dic_ace_aa[aa][1],
                          dic_up_aa[aa][1] - dic_cce_aa[aa][1],
                          codon + codon2rareness[codon], dic_up_codon[codon][0],
                          dic_up_codon[codon][1],
                          dic_ace_codon[codon][0], dic_ace_codon[codon][1], dic_cce_codon[codon][0],
                          dic_cce_codon[codon][1], dic_up_codon[codon][1] - dic_ace_codon[codon][1],
                          dic_up_codon[codon][1] - dic_cce_codon[codon][1]]
            main_content.append(my_row)
    merging_value.append(merging_value[-1] + 2)
    return main_content, merging_value, merging_content


def get_content_second_sheet(dic_up_composition, dic_up_last_nt, dic_down_composition, dic_down_last_nt,
                             dic_up_down_composition, dic_up_down_last_nt, dic_ace_composition, dic_ace_last_nt,
                             dic_cce_composition, dic_cce_last_nt):
    """
    This function allow to create the content of the sheets names "meta_exon_composition" and weighted_exon_composition"
    in the composition_summary.xlsx file from many dictionaries. Those previously mentioned dictionaries are described
    below :
    :param dic_up_composition: (dictionary of list of 2 floats) for each of the following nucleotides or group of
    nucleotides (A, T, G, C, AT, GC, purine, pyrimidine) in the up-regulated user set of exon, give its number of
    occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :param dic_up_last_nt:(dictionary of list of 2 floats) for each of the following nucleotides (A, T, G, C)  give its
    number of occurrence and its frequency at the last position of codons of the up-regulated user set of exon
    (calculated with the meta-exon method or with the weighted method)
    :param dic_down_composition:(dictionary of list of 2 floats) for each of the following nucleotides or group of
    nucleotides (A, T, G, C, AT, GC, purine, pyrimidine) in the down-regulated user set of exon, give its number of
    occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :param dic_down_last_nt:(dictionary of list of 2 floats) for each of the following nucleotides (A, T, G, C) give its
    number of occurrence and its frequency at the last position of codons of the down-regulated user set of exon
    (calculated with the meta-exon method or with the weighted method)
    :param dic_up_down_composition:(dictionary of list of 2 floats) for each of the following nucleotides or group of
    nucleotides (A, T, G, C, AT, GC, purine, pyrimidine) in the down-regulated and up_regulated user set of exon,
    give its number of occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :param dic_up_down_last_nt:(dictionary of list of 2 floats) for each of the following nucleotides (A, T, G, C) give
    its number of occurrence and its frequency at the last position of codons of the down-regulated and up_regulated
    user set of exon (calculated with the meta-exon method or with the weighted method)
    :param dic_ace_composition:(dictionary of list of 2 floats) for each of the following nucleotides or group of
    nucleotides (A, T, G, C, AT, GC, purine, pyrimidine) in the ACE control set of exon,
    give its number of occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :param dic_ace_last_nt:(dictionary of list of 2 floats) for each of the following nucleotides (A, T, G, C) give
    its number of occurrence and its frequency at the last position of codons of the ACE control set of exon
    (calculated with the meta-exon method or with the weighted method)
    :param dic_cce_composition:(dictionary of list of 2 floats) for each of the following nucleotides or group of
    nucleotides (A, T, G, C, AT, GC, purine, pyrimidine) in the CCE control set of exon,
    give its number of occurrence and its frequency calculated with the meta-exon method or with the weighted method
    :param dic_cce_last_nt:(dictionary of list of 2 floats) for each of the following nucleotides (A, T, G, C) give
    its number of occurrence and its frequency at the last position of codons of the CCE control set of exon
    (calculated with the meta-exon method or with the weighted method)
    :return:  1 - main_content :  a list of list corresponding to the content of the second or fourth sheet of the
    composition_summary file. Each sublist within the main list corresponds to a line of the 1st or the 3rd sheet of
    the composition_summary file.
              2 - merging_value : list of string. It corresponds to the cells we want to merge in the
              composition_summary file
              3 - merging_content : list of strings that correspond to the content of the merged cells
    """
    # initialisation of some variables
    merging_content = None
    list_comp = ['A>=2', 'C>=2', 'G>=2', 'T>=2', 'AC>=2', 'AG>=2', 'AT>=2', 'CG>=2', 'CT>=2', 'GT>=2', 'A=0', 'C=0', 'G=0', 'T=0', 'AC=0', 'AG=0', 'AT=0', 'CG=0', 'CT=0', 'GT=0']
    content = list()
    header1 = None
    merging_value = None
    # filling the headers and the merging_value and merging_content variables
    if dic_down_composition is not None and dic_up_down_composition is not None:
        header1 = ["codon_composition", "up_exon", '  ', "down_exon", ' ', "up_down_exon", '', "ACE(up)", '',
                   "CCE(up)", '', "up_vs_down", "up_vs_ACE", "down_vs_ACE", "up_down_vs_ACE", "up_vs_CCE", "down_vs_CCE",
                   "up_down_vs_CCE"]
        merging_value = ["B1:C1", "D1:E1", "F1:G1", "H1:I1", "J1:K1", "B22:C22", "D22:E22", "F22:G22", "H22:I22",
                         "J22:K22"]
        merging_content = [header1[1], header1[3], header1[5], header1[7], header1[9]] * 2

    elif dic_down_composition is not None and dic_up_down_composition is None:
        header1 = ["codon_composition", "up_exon", '  ', "down_exon", ' ', "ACE(up)", '', "CCE(up)",
                   '', "up_vs_down", "up_vs_ACE",
                   "down_vs_ACE", "up_vs_CCE", "down_vs_CCE"]
        merging_value = ["B1:C1", "D1:E1", "F1:G1", "H1:I1", "B22:C22", "D22:E22", "F22:G22", "H22:I22"]
        merging_content = [header1[1], header1[3], header1[5], header1[7]] * 2
    elif dic_down_composition is None and dic_up_down_composition is None:
        header1 = ["codon_composition", "input_exon", '  ', "ACE(up)", '', "CCE(up)",
                   '', "input_vs_ACE", "input_vs_CCE"]
        merging_value = ["B1:C1", "D1:E1", "F1:G1", "B22:C22", "D22:E22", "F22:G22"]
        merging_content = [header1[1], header1[3], header1[5]] * 2
    header2 = ["last_nt"] + header1[1:]
    content.append(header1)
    # getting the content of the sheet
    # 1 - First part : codon rich in A, T, C, G, AT ,GC, purine, pyrimidine
    for comp in list_comp:
        my_row = None
        if dic_down_composition is not None and dic_up_down_composition is not None:
            my_row = [comp, dic_up_composition[comp][0], dic_up_composition[comp][1],
                      dic_down_composition[comp][0], dic_down_composition[comp][1], dic_up_down_composition[comp][0],
                      dic_up_down_composition[comp][1], dic_ace_composition[comp][0], dic_ace_composition[comp][1],
                      dic_cce_composition[comp][0], dic_cce_composition[comp][1],
                      dic_up_composition[comp][1] - dic_down_composition[comp][1],
                      dic_up_composition[comp][1] - dic_ace_composition[comp][1],
                      dic_down_composition[comp][1] - dic_ace_composition[comp][1],
                      dic_up_down_composition[comp][1] - dic_ace_composition[comp][1],
                      dic_up_composition[comp][1] - dic_cce_composition[comp][1],
                      dic_down_composition[comp][1] - dic_cce_composition[comp][1],
                      dic_up_down_composition[comp][1] - dic_cce_composition[comp][1]]
        elif dic_down_composition is not None and dic_up_down_composition is None:
            my_row = [comp, dic_up_composition[comp][0], dic_up_composition[comp][1],
                      dic_down_composition[comp][0], dic_down_composition[comp][1], dic_ace_composition[comp][0],
                      dic_ace_composition[comp][1], dic_cce_composition[comp][0], dic_cce_composition[comp][1],
                      dic_up_composition[comp][1] - dic_down_composition[comp][1],
                      dic_up_composition[comp][1] - dic_ace_composition[comp][1],
                      dic_down_composition[comp][1] - dic_ace_composition[comp][1],
                      dic_up_composition[comp][1] - dic_cce_composition[comp][1],
                      dic_down_composition[comp][1] - dic_cce_composition[comp][1]]
        elif dic_down_composition is None and dic_up_down_composition is None:
            my_row = [comp, dic_up_composition[comp][0], dic_up_composition[comp][1],
                      dic_ace_composition[comp][0], dic_ace_composition[comp][1], dic_cce_composition[comp][0],
                      dic_cce_composition[comp][1], dic_up_composition[comp][1] - dic_ace_composition[comp][1],
                      dic_up_composition[comp][1] - dic_cce_composition[comp][1]]
        content.append(my_row)
    content.append(header2)
    # 1 - second part: codon rich in A, T, C, G, AT ,GC, purine, pyrimidine
    for nt in dic_up_last_nt.keys():
        my_row = None
        if dic_down_composition is not None and dic_up_down_composition is not None:
            my_row = [nt, dic_up_last_nt[nt][0], dic_up_last_nt[nt][1], dic_down_last_nt[nt][0],
                      dic_down_last_nt[nt][1], dic_up_down_last_nt[nt][0], dic_up_down_last_nt[nt][1],
                      dic_ace_last_nt[nt][0], dic_ace_last_nt[nt][1], dic_cce_last_nt[nt][0],
                      dic_cce_last_nt[nt][1], dic_up_last_nt[nt][1] - dic_down_last_nt[nt][1],
                      dic_up_last_nt[nt][1] - dic_ace_last_nt[nt][1],
                      dic_down_last_nt[nt][1] - dic_ace_last_nt[nt][1],
                      dic_up_down_last_nt[nt][1] - dic_ace_last_nt[nt][1],
                      dic_up_last_nt[nt][1] - dic_cce_last_nt[nt][1],
                      dic_down_last_nt[nt][1] - dic_cce_last_nt[nt][1],
                      dic_up_down_last_nt[nt][1] - dic_cce_last_nt[nt][1]]
        elif dic_down_composition is not None and dic_up_down_composition is None:
            my_row = [nt, dic_up_last_nt[nt][0], dic_up_last_nt[nt][1], dic_down_last_nt[nt][0],
                      dic_down_last_nt[nt][1],
                      dic_ace_last_nt[nt][0],
                      dic_ace_last_nt[nt][1],
                      dic_cce_last_nt[nt][0], dic_cce_last_nt[nt][1], dic_up_last_nt[nt][1] - dic_down_last_nt[nt][1],
                      dic_up_last_nt[nt][1] - dic_ace_last_nt[nt][1],
                      dic_down_last_nt[nt][1] - dic_ace_last_nt[nt][1],
                      dic_up_last_nt[nt][1] - dic_cce_last_nt[nt][1], dic_down_last_nt[nt][1] - dic_cce_last_nt[nt][1]]
        elif dic_down_composition is None and dic_up_down_composition is None:
            my_row = [nt, dic_up_last_nt[nt][0], dic_up_last_nt[nt][1],
                      dic_ace_last_nt[nt][0],
                      dic_ace_last_nt[nt][1],
                      dic_cce_last_nt[nt][0], dic_cce_last_nt[nt][1], dic_up_last_nt[nt][1] - dic_ace_last_nt[nt][1],
                      dic_up_last_nt[nt][1] - dic_cce_last_nt[nt][1]]
        content.append(my_row)
    return content, merging_value, merging_content


def get_dimention(sheet_list):
    """
    :param sheet_list: (list of list of integer and string) each list of the main list corresponds to a line that
    must be written in the excel result file
    :return: (string) the size of the table "sheet list" described using the cells identifiers in a excel file
    """
    letters = string.ascii_uppercase
    col_len = len(sheet_list[0])
    row_len = len(sheet_list)
    return "A1:" + letters[col_len - 1] + str(row_len + 1)


def formating_merging_value(merging_value, content, merging_content_input, leng_merg):
    """
    :param merging_value: list of string and int, indicate the cells we want to merge together. The string part of
    the list is already ok but the int part must be processed by this function for merging
    :param content: (list of list) the content of the sheet (will be used to get the information we need to fill
    the merged cells)
    :param merging_content_input: list of string : the first part of the content oof the cells to merged.
    :param leng_merg: int, allow to merge the correct cells : indeed the cells to merge will be different if
    the user don't want the information on the down, the up or the up and down sets of exon
    :return: 1 - merging_value_formated: list of string : the correct format of the cells to merge together
             2 - merging_content: (list of string) the content of the cells to merge
    """
    merging_value_formated = list()
    merging_content = list()
    k = 0
    for k in range(len(merging_content_input)):
        my_str = merging_value[k]
        merging_content.append(merging_content_input[k])
        merging_value_formated.append(my_str)
    for i in range(k + 1, len(merging_value) - 1):
        for j in range(len(string.letters[0:leng_merg].upper())):
            my_str = str(string.letters[0:leng_merg].upper()[j]) + str(merging_value[i]) + ":" + \
                     str(string.letters[0:leng_merg].upper()[j]) + str(merging_value[i + 1] - 1)
            merging_content.append(content[merging_value[i] - 1][j])
            merging_value_formated.append(my_str)

    return merging_value_formated, merging_content


def sheet_filer_aa_codon(workbook, your_sheet, content, merging_value, merging_content_input):
    """
    :param workbook: (Workbook object from xlsxwriter) => use to create the file composition_summary
    :param your_sheet: (Sheet object from xlsxwriter) => use to fill the sheet
    :param content:  (list of list) the content of the sheet
    :param merging_value: (list of string) the cells we want to merge together
    :param merging_content_input: (list of string) the content of the string merged together
    """
    merge_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': True})
    header_color = workbook.add_format({'bg_color': '#00DCFF', 'align': 'center', 'valign': 'vcenter', 'border': True})
    red = workbook.add_format()
    red.set_bg_color('#FF0000')
    green = workbook.add_format()
    green.set_bg_color('#00CC00')

    i = 1
    row_name = "A" + str(i)
    for row in content:
        if i == 1:
            your_sheet.write_row(row_name, row, header_color)
        else:
            your_sheet.write_row(row_name, row, merge_format)
        i += 1
        row_name = "A" + str(i)
    if len(content[0]) == 36:
        leng_merg = 18
    elif len(content[0]) == 28:
        leng_merg = 14
    else:
        leng_merg = 9
    # merging procedure...
    merging_value_formated, merging_content = formating_merging_value(merging_value, content, merging_content_input,
                                                                      leng_merg)
    j = 0
    for j in range(len(merging_content_input)):
        your_sheet.merge_range(merging_value_formated[j], merging_content_input[j], header_color)
    for i in range(j + 1, len(merging_value_formated)):
        your_sheet.merge_range(merging_value_formated[i], merging_content[i], merge_format)

    # resize cells
    your_sheet.set_column(0, 0, 11)
    if len(content[0]) == 36:
        your_sheet.set_column(11, 17, 15)
        your_sheet.set_column(29, 35, 15)
    if len(content[0]) == 28:
        your_sheet.set_column(9, 13, 15)
        your_sheet.set_column(23, 26, 15)
    if len(content[0]) == 18:
        your_sheet.set_column(7, 8, 15)
        your_sheet.set_column(16, 17, 15)
    # freezing panes
    your_sheet.freeze_panes(1, 1)

    # applying conditional coloration
    if len(content[0]) == 36:
        cell_list = ['L2:L65', 'M2:M65', 'N2:N65', 'O2:O65', 'P2:P65', 'Q2:Q65', 'R2:R65', 'AD2:AD65', 'AE2:AE65', 'AF2:AF65',
                     'AG2:AG65', 'AH2:AH65', 'AI2:AI65', 'AJ2:AJ65']
    elif len(content[0]) == 28:
        cell_list = ['J2:J65', 'K2:K65', 'L2:L65', 'M2:M65', 'N2:N65', 'X2:X65', 'Y2:Y65', 'Z2:Z65',
                     'AA2:AA65', 'AB2:AB65']
    else:
        cell_list = ['H2:H65', 'I2:I65', 'Q2:Q65', 'R2:R65']
    for cell in cell_list:
        your_sheet.conditional_format(cell, {'type': 'cell',
                                             'criteria': 'between',
                                             'minimum': 1,
                                             'maximum': 100,
                                             'format': green})
        your_sheet.conditional_format(cell, {'type': 'cell',
                                             'criteria': 'between',
                                             'minimum': -100,
                                             'maximum': -1,
                                             'format': red})


def sheet_filer_exon_composition(workbook, your_sheet, content, merging_values, merging_content_input):
    """
    :param workbook: (Workbook object from xlsxwriter) => use to create the file composition_summary
    :param your_sheet: (Sheet object from xlsxwriter) => use to fill the sheet
    :param content:  (list of list) the content of the sheet
    :param merging_values: (list of string) the cells we want to merge together
    :param merging_content_input: (list of string) the content of the string merged together
    """
    merge_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': True})
    header_color = workbook.add_format({'bg_color': '#00DCFF', 'align': 'center', 'valign': 'vcenter', 'border': True})
    red = workbook.add_format()
    red.set_bg_color('#FF0000')
    green = workbook.add_format()
    green.set_bg_color('#00CC00')

    i = 1
    row_name = "A" + str(i)
    for row in content:
        if i == 1 or i == 22:
            your_sheet.write_row(row_name, row, header_color)
        else:
            your_sheet.write_row(row_name, row, merge_format)
        i += 1
        row_name = "A" + str(i)

    # merging procedure
    merging_content = merging_content_input
    for i in range(len(merging_values)):
        your_sheet.merge_range(merging_values[i], merging_content[i], header_color)

    # resize the column content
    your_sheet.set_column(0, 0, 20)

    if len(content[0]) == 18:
        your_sheet.set_column(11, 17, 20)
    elif len(content[0]) == 14:
        your_sheet.set_column(9, 13, 20)
    else:
        your_sheet.set_column(7, 8, 20)
    # freezing panes
    your_sheet.freeze_panes(0, 1)

    # conditional coloration
    if len(content[0]) == 18:
        cell_list = ['L2:L26', 'M2:M26', 'N2:N26', 'O2:O26', 'P2:P26', 'Q2:Q26', 'R2:R26']
    elif len(content[0]) == 14:
        cell_list = ['J2:J26', 'K2:K26', 'L2:L26', 'M2:M26','N2:N26']
    else:
        cell_list = ['H2:H14', 'I2:I14']
    for cell in cell_list:
        your_sheet.conditional_format(cell, {'type': 'cell',
                                             'criteria': 'between',
                                             'minimum': 1,
                                             'maximum': 100,
                                             'format': green})
        your_sheet.conditional_format(cell, {'type': 'cell',
                                             'criteria': 'between',
                                             'minimum': -100,
                                             'maximum': -1,
                                             'format': red})


def write_xls(output_folder, content_first_sheet, merging_value, merging_content_first_sheet, content_second_sheet,
              merging_value_second_sheet, merging_content_second_sheet, content_third_sheet, merging_value_third_sheet,
              merging_content_third_sheet, content_fourth_sheet, merging_value_fourth_sheet,
              merging_content_fourth_sheet):
    """
    :param output_folder: (string), the output folder where yuo want t create the composition_summary file
    :param content_first_sheet: (list of list), the content of the first sheet
    :param merging_value: (list of string), the cells we want to merge together for the first sheet
    :param merging_content_first_sheet:  (list of string) the content of the string we want to merge together
    :param content_second_sheet:(list of list), the content of the second sheet
    :param merging_value_second_sheet:(list of string), the cells we want to merge together for the second sheet
    :param merging_content_second_sheet:(list of string) the content of the string we want to merge together
    :param content_third_sheet: (list of list), the content of the third sheet
    :param merging_value_third_sheet:(list of string) the content of the string we want to merge together
    :param merging_content_third_sheet:(list of string) the content of the string we want to merge together
    :param content_fourth_sheet: (list of list), the content of the fourth sheet
    :param merging_value_fourth_sheet:(list of string) the content of the string we want to merge together
    :param merging_content_fourth_sheet:(list of string) the content of the string we want to merge together
    """
    # Creation of the document
    workbook = xlsxwriter.Workbook(output_folder + "composition_summary.xlsx")
    meta_sheet = workbook.add_worksheet("meta_exon_recap")
    sheet_filer_aa_codon(workbook, meta_sheet, content_first_sheet, merging_value, merging_content_first_sheet)
    meta_sheet2 = workbook.add_worksheet("meta_exon_composition")
    sheet_filer_exon_composition(workbook, meta_sheet2, content_second_sheet, merging_value_second_sheet,
                                 merging_content_second_sheet)
    weighted_sheet3 = workbook.add_worksheet("weighted_exon_recap")
    sheet_filer_aa_codon(workbook, weighted_sheet3, content_third_sheet, merging_value_third_sheet,
                         merging_content_third_sheet)
    weighted_sheet4 = workbook.add_worksheet("weighted_exon_composition")
    sheet_filer_exon_composition(workbook, weighted_sheet4, content_fourth_sheet, merging_value_fourth_sheet,
                                 merging_content_fourth_sheet)
    workbook.close()


def getting_summary_dictionaries(exon_list, penalty_size):
    """
    :param exon_list: (list of exon object defined in exon_class.py), an exon list
    :param penalty_size: (int) the size below which an exon is considerate as small
    :return:
        1 - dic_codon : dictionary of list of floats, give for each codon, its proportion and its number of occurence in
        the exon_list (calculated with the meta-exon method)
        2 - dic_aa : dictionary of list of floats, give for each aa, its proportion and its number of occurence in
        the exon_list (calculated with the meta-exon method)
        3 - dic_composition : dictionary of list of floats, give the number of exon rich in  A, T, C, G, AT, GC,
        purine, pyrimidine and their frequencies in the exon_list  (calculated with the meta-exon method)
        4 - dic_last_nt : dictionary of list of floats : give for each nucleotide the number and the frequency of codons
        that have their 3rd base corresponding to that nucleotide  (calculated with the meta-exon method)
        5 - weighted_aa: dictionary of list of floats, give for each codon, its proportion and its number of occurence
        in the exon_list (calculated with the weighted method)
        6 - weighted_codon:dictionary of list of floats, give for each aa, its proportion and its number of occurence in
        the exon_list (calculated with the weighted method)
        7 - weighted_composition:dictionary of list of floats, give the number of exon rich in  A, T, C, G, AT, GC,
        purine, pyrimidine and their frequencies in the exon_list  (calculated with the weighted method)
        8 - weighted_last_nt: dictionary of list of floats : give for each nucleotide the number and the frequency of
        codons that have their 3rd base corresponding to that nucleotide  (calculated with the weighted method)
    """
    dic_aa = exon_list.meta_exon_aa_calculator()
    dic_codon = exon_list.meta_exon_codon_calculator()
    dic_composition = exon_list.meta_exon_composition_analyser()
    dic_last_nt = exon_list.meta_exon_last_nt_proportion_calculator()
    weighted_aa = exon_list.amino_acid_frequency_calculator_recap(penalty_size / 3)
    weighted_codon = exon_list.weighted_frequency_calculator_recap(penalty_size / 3)
    weighted_composition = exon_list.weighted_exon_composition_recap(penalty_size / 3)
    weighted_last_nt = exon_list.weighted_last_nt_frequency_calculator_recap(penalty_size / 3)
    return dic_codon, dic_aa, dic_composition, dic_last_nt, weighted_codon, weighted_aa, weighted_composition, \
           weighted_last_nt

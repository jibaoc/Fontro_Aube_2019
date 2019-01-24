#!/usr/bin/python3


# The goal of this script if to calculate the statistics of the Fig4B of the article interplay between \
# coding and exonic  splicing regulatory sequence

# imports
import pandas as pd
import subprocess
import rpy2.robjects
import rpy2.robjects.vectors as v

# global variables:

folder = os.path.realpath(os.path.dirname(__file__).replace("src", "data/")
output = os.path.realpath(os.path.dirname(__file__).replace("src", "result/")

feature_dic = {
    "Very-small": ["A", "C", "G", "S"],
    "Small#2": ["A", "C", "D", "G", "N", "P", "S", "T"],
    "Large" : ["F", "I", "K", "L", "M", "R", "W", "Y"],
    "Disorder-promoting#1": ["A", "E", "G", "K", "P", "Q", "R", "S"],
    "Order-promoting#1": ["C", "F", "I", "L", "N", "W", "V", "Y"],
    "Disorder-promoting#2": ["A", "E", "G", "K", "P", "Q", "S"],
    "Order-promoting#2": ["C", "F", "H", "I", "L", "M", "N", "W", "V", "Y"],
    "Polar-uncharged#1": ["C", "N", "Q", "S", "T", "Y"],
    "Polar-uncharged#2": ["N", "Q", "S", "T", "Y"],
    "Charged#1": ["R", "H", "K", "D", "E"],
    "Charged#2": ["R", "K", "D", "E"],
    "Hydrophilic#1": ["D", "E", "K", "N", "Q", "R"],
    "Hydrophobic#1": ["A", "C", "F", "I", "L", "M", "V"],
    "Neutral": ["G", "H", "P", "S", "T", "Y"],
    "Hydroxylic": ["S", "T", "Y"],
    "Negatively-charged": ["D", "E"],
    "Positively-charged#1": ["R", "H", "K"],
    "Positively-charged#2": ["R", "K"]
}

ft_i = ['Very-small', 'Neutral', 'Positively-charged#2', 'Polar-uncharged#2',
              'Negatively-charged', 'Large', 'Hydrophilic#1', 'Charged#2', 'Hydrophobic#1']


# functions


def read_sequence(excel_file):
    """
    :param excel_file: (string) an excel_file containing exons sequence
    :return: (list of strings) the list of amino acid sequence of interest exons
    """
    # opening the excel file
    xl = pd.ExcelFile(excel_file)
    df = "NA"
    # opening the sheet of interest
    for sheet in xl.sheet_names:
        if "sequence" == sheet:
            df = xl.parse(sheet)
    # if the sheet "exon_skipping*" doesn't exist, the faRLine file cannot be used so we stop the program
    if str(df) == "NA":
        print("the sheet names sequence wasn't found")
        print("exiting...")
        exit(1)

    pep = []
    for row in df.itertuples():
        if not isinstance(row.CDS_peptide_sequence, float):
            if len(row.CDS_peptide_sequence.replace("*","")) >= 1:
                pep.append(row.CDS_peptide_sequence.replace("*",""))
    return pep


def frequency_calculator_ft(list_seq, feature):
    """
    :param list_seq: (list of string) list of peptide sequence
    :param feature: (string) a feature
    :return: (list of float) the frequency of every sequence in list_seq for teh feature
    """
    freq_ft = []
    for seq in list_seq:
        count = 0
        for aa in seq:
            if aa in feature_dic[feature]:
                count += 1
        freq_ft.append(float(count) / len(seq))
    return freq_ft


def couple_up_and_down_file_finder(folder):
    """
    From a folder containing folders that contains an up and down regulated list of exons, return a \
    dictionary that links each folder (corresponding to a project) to a list containing the files \
    with the up-regulated and down regulated exons.

    :param folder: (string) A folder containing subfolders containing to files: one file corresponds to \
    the up-regulated exons files and is named 'query_results_up' the other contains the down regulated exons \
    and is named 'query_results_down'
    :return: (dictionary of list of string) links each project (subfolders of `folder`) to the file it contains.
    """
    list_file = subprocess.check_output(["find", folder, "-name", "*.xlsx", "-type", "f"],
                                        stderr=subprocess.STDOUT).decode("ascii").split("\n")[:-1]
    dic = {}
    for myfile in list_file:
        sf = myfile.split("/")[-2].split("_")[2]
        if sf not in dic.keys():
            dic[sf] = [myfile]
        else:
            dic[sf].append(myfile)
    return dic


def dic_frequency(dic_file, feature):
    """
    :param dic_file: (dictionary of list of string) links each project (subfolders) to the file it contains.
    :param feature: (string) the feature of interest
    :return: (dictionary of list of float) link each project the the list of frequency of `feature` \ for the \
    up and down-regulated set of exons for this project
    """
    new_dic = {}
    for key in dic_file.keys():
        list_seq = read_sequence(dic_file[key][0])
        new_dic[key] = [frequency_calculator_ft(list_seq, feature)]
        list_seq = read_sequence(dic_file[key][1])
        new_dic[key].append(frequency_calculator_ft(list_seq, feature))
    return new_dic


def shapiro_test(list_val):
    """
    :param list_val: (list of float) list of frequency for a feature corresponding to a set of exon
    :return: (float) p-value of the shapiro-wilk test
    """

    shapiro_test = rpy2.robjects.r(
        """
        function(list_val){
            return(shapiro.test(list_val)$p.value)
        }
        """
    )
    return shapiro_test(v.FloatVector(list_val))[0]


def comparison_test(list_val1, list_val2, test):
    """

    :param list_val1: (list of float) list of frequency for a feature corresponding to a set of exon
    :param list_val2: (list of float) list of frequency for a feature corresponding to another set of exon
    :param test: (string) the type of test to use
    :return: (float) pvalue of the comparison test
    """
    ttest = rpy2.robjects.r(
        """
        function(list_val1, list_val2){
            return(t.test(list_val1, list_val2)$p.value)
        }
        """
    )
    wilcox = rpy2.robjects.r(
        """
        function(list_val1, list_val2){
            return(wilcox.test(list_val1, list_val2)$p.value)
        }
        """
    )
    if test == "wilcoxon":
        return wilcox(v.FloatVector(list_val1), v.FloatVector(list_val2))[0]
    else:
        return ttest(v.FloatVector(list_val1), v.FloatVector(list_val2))[0]


def main():
    dic_file = couple_up_and_down_file_finder(folder)
    test = "t-test"
    # Test what test we want to do
    for feature in ft_i:
        dic_freq = dic_frequency(dic_file, feature)
        for key in dic_freq.keys():
            if shapiro_test(dic_freq[key][0]) < 0.05 or shapiro_test(dic_freq[key][1]) < 0.05:
                test = "wilcoxon"

    with open(output + "stat_Fig6B.csv", "w") as outfile:
        outfile.write("Up vs Down comparison : test used : " + str(test) + "\n")
        outfile.write("Feature\tSF_Union\tPvalue\n")
        for feature in ft_i:
            dic_freq = dic_frequency(dic_file, feature)

            for sf in dic_freq.keys():
                pval = comparison_test(dic_freq[sf][0], dic_freq[sf][1], test)
                outfile.write(str(feature) + "\t" + str(sf) + "\t" + str(pval) + "\n")

if __name__ == "__main__":
    main()

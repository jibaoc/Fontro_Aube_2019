#!/usr/bin/python3

# Calculates the pvalue of somme amino acids and dinucleotide (man withney test)

# imports

import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import pandas as pd
import subprocess

# Global variables

# the statistical analysis for those di nucleotide will be perform
dnt_i = ["CA", "AC", "CT", "TC"]
# the statistical analysis for those amino acid will be perform
aa_i = ["H", "Q", "T", "S", "L"]

# folder interest
folder="/media/nicolas/DD_2/Projects/Fig2H_stat_calculator/data"

# SF of interest
sf_i = ["PTBP1", "hnRNPL"]

# output
output="/media/nicolas/DD_2/Projects/Fig2H_stat_calculator/result/"

# function

def file_finder(folder, sf_i):
    """
    :param fasta: the folder containing enrichment file of interest
    :param sf_i: the list containing the 2 sf of interest
    :return: 2 lists, one containing the file related to the first splicing factor of interest and the second \
    the files related to the second splicing factor of interest
    """
    a = subprocess.check_output(["find", folder, "-name", "enrichment_report.xlsx", "-type", "f"])
    a = a.decode("ascii").split("\n")[:-1]

    listsf1 = []
    listsf2= []
    for mfile in a:
        if sf_i[0] in mfile:
            listsf1.append(mfile)
        else:
            listsf2.append(mfile)
    return listsf1, listsf2



def get_interest_frequency(excel_file, dic):
    """
    :param excel_file: (string) path to a query file given by the tRNA program
    :param dic: (dictionary) contains the value for the interets dinucleotide and amino acids
    :return: (2 list of strings) the list of cds sequence and the list of amino acid sequence
    """
    # opening the excel file
    xl = pd.ExcelFile(excel_file)
    df = "NA"
    # opening the sheet of interest
    for sheet in xl.sheet_names:
        if "amino_acid" == sheet:
            df = xl.parse(sheet)
    # if the sheet "amino_acid" doesn't exist, the faRLine file cannot be used so we stop the program
    if str(df) == "NA":
        print("the sheet names amino_acid wasn't found")
        print("exiting...")
        exit(1)

    for row in df.itertuples():
        if row.amino_acid in aa_i:
            val = float(row.frequencies_of_the_interest_set) * 100
            if row.amino_acid not in dic.keys():
                dic[row.amino_acid] = [val]
            else:
                dic[row.amino_acid].append(val)

    # opening the sheet of interest
    for sheet in xl.sheet_names:
        if "dnt_info" == sheet:
            df = xl.parse(sheet)
    # if the sheet "amino_acid" doesn't exist, the faRLine file cannot be used so we stop the program
    if str(df) == "NA":
        print("the sheet names dnt_info wasn't found")
        print("exiting...")
        exit(1)

    dic_test = {i:0 for i in dnt_i + aa_i}
    for row in df.itertuples():
        if row.dnt_info in dnt_i:
            val = float(row.frequencies_of_the_interest_set) * 100
            if dic_test[row.dnt_info] == 0:
                if row.dnt_info not in dic.keys():
                    dic[row.dnt_info] = [val]
                else:
                    dic[row.dnt_info].append(val)
            dic_test[row.dnt_info] += 1

    return dic


def pvalue_getter(list1, list2):
    """
    :param list1: list of float
    :param list2: list of float
    :return: a p value
    """
    print(list1)
    print(list2)
    wilcox = robj.r(
        """
        function(list1, list2){
            return(wilcox.test(list1,list2)$p.value)
        }
        """)
    return(float(wilcox(v.FloatVector(list1), v.FloatVector(list2))[0]))


def main():
    listsf1, listsf2 = file_finder(folder, sf_i)
    dic_1 = {}
    dic_2 = {}
    for myfile in listsf1:

        dic_1 = get_interest_frequency(myfile, dic_1)
    for myfile in listsf2:
        dic_2 = get_interest_frequency(myfile, dic_2)

    with open(output + "Fig5E.stat.txt", "w") as outfile:
        outfile.write("Comparaison des SF " + sf_i[0] + " et " + sf_i[1] + "\n")
        outfile.write("order of file for " + sf_i[0] + str(listsf1) + "\n")
        outfile.write("order of file for " + sf_i[1] + str(listsf2) + "\n")
        for key in dic_1.keys():
            pval = pvalue_getter(dic_1[key], dic_2[key])
            outfile.write(str(key) + " - Valeurs " + str(sf_i[0])  + " : "  + str(dic_1[key]) + "Valeurs " + str(sf_i[1])  + " : " + str(dic_2[key]) + "  - pval = " + str(pval) + "\n")


if __name__=="__main__":
    main()
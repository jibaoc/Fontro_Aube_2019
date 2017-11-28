import subprocess
import argparse
import os
import pandas as pd
from Bio import SeqIO


def fasta_maker(dnt_list, output, name_file):
    """
    Create a fasta file with the given dnt
    :param dnt_list: (dic of float value) : the proportion (value) of each dnt (key)
    :param output: (string) the folder where the fasta file will be created
    :param name_file: (string) the name of the fasta file
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    subprocess.check_call(["python3", dir_path + "/fasta_generator.py", "--AA", dnt_list[0], "--AC", dnt_list[1],
                           "--AG", dnt_list["AG"], "--AT", dnt_list["AT"], "--CA", dnt_list["CA"], "--CC",
                           dnt_list["CC"], "--CG", "--CT", dnt_list["CT"], "--GA", dnt_list["GA"], "--GC",
                           dnt_list["GC"], "--GG", dnt_list["GG"], "--GT", dnt_list["GT"], "--TA", dnt_list["TA"],
                           "--TC", dnt_list["TC"], "--TG", dnt_list["TG"],
                           "--output", output, "--filename", name_file, "--size_max", "270"],
                          stderr=subprocess.STDOUT)


def proportion_getter(excel_file):
    """
    :param excel_file: (string) an excel file corresponding to an enrichment file produce by the tRNA program
    :return: (dic of float value) : the proportion (value) of each dnt (key)
    """
    dnt_list = {"AA": 0., "AC": 0., "AG": 0., "AT": 0., "CA": 0., "CC": 0., "CG": 0., "CT": 0., "GA": 0., "GC": 0.,
                "GG": 0., "GT": 0., "TA": 0., "TC": 0., "TG": 0.}
    xl = pd.ExcelFile(excel_file)
    df = "NA"
    for sheet in xl.sheet_names:
        if "dnt_info" == sheet:
            df = xl.parse(sheet)
    # if the sheet nt info doesn't exist, we end the program
    if str(df) == "NA":
        print("The sheet dnt_info was not fing in " + str(excel_file))
        print("Terminating...")
        exit(1)
    for row in df.itertuples():
        if row.dnt_info in dnt_list.keys():
            dnt_list[row.dnt_info] = row.frequencies_of_the_interest_set
    return dnt_list


def fasta_creator(input_xls, bg_xls, output):
    """
    :param input_xls: (string) the excel input file
    :param bg_xls: (string) the excel background file
    """
    dnt_list = proportion_getter(input_xls)
    fasta_maker(dnt_list, output, "input")
    dnt_list = proportion_getter(bg_xls)
    fasta_maker(dnt_list, output, "background")


def calcul_dic(dic, seq, length, count):
    """
    :param dic: (a dictionary of float) freq of the word af interest
    :param seq: (string) the nt sequence
    :param length: (int) the length of the word that we want to find
    :param count: (int) the number of words in the dic
    :return: 'dic' completed
    """
    for i in range(len(seq)-(length-1)):
        if seq[i:i+length] not in dic.keys():
            dic[seq[i:i+length]] = 1
            count += 1
        else:
            dic[seq[i:i+length]] += 1
            count += 1
    return dic, count

def get_frequencies(dic, count):
    """
    :param dic: (dictionary of int) each key correspond to an oligo-mer and it's value : the number of time it was
    seen in the dictionary
    :param count: (int) the total number of word in that dictionary
    :return: (dic of float) the frequencies of each word
    """
    for key in dic.keys():
        dic[key] = float(dic[key]) / count
    return dic

def counting_letters_frequencies(fasta_file):
    """
    :param fasta_file: (string)  a fasta file
    :return: (5 dic of float):
         - dic_6l : link each word of 6 letters found in the fasta files with its frequencies in the fasta sequence
         - dic_5l : link each word of 5 letters found in the fasta files with its frequencies in the fasta sequence
         - dic_4l : link each word of 4 letters found in the fasta files with its frequencies in the fasta sequence
         - dic_3l : link each word of 3 letters found in the fasta files with its frequencies in the fasta sequence
    """
    dic_6l = {}
    count_6 = 0
    dic_5l = {}
    count_5 = 0
    dic_4l = {}
    count_4 = 0
    dic_3l = {}
    count_3 = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        dic_6l, count_6 = calcul_dic(dic_6l, record.seq, 6, count_6)
        dic_5l, count_5 = calcul_dic(dic_5l, record.seq, 6, count_5)
        dic_4l, count_4 = calcul_dic(dic_4l, record.seq, 6, count_4)
        dic_3l, count_3 = calcul_dic(dic_3l, record.seq, 6, count_3)

    dic_6l = get_frequencies(dic_6l, count_6)
    dic_5l = get_frequencies(dic_5l, count_5)
    dic_4l = get_frequencies(dic_4l, count_4)
    dic_3l = get_frequencies(dic_3l, count_3)

    return dic_6l, dic_5l, dic_4l, dic_3l






def main():
    "This is the main function that will parse the given argument to the program"
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""this script takes an 'input' and a 'background' trna enrichment
                                     report file. For each of this files, it will create in a given output directory
                                     a fasta file with random sequences having the same dinucleotides proportions than
                                     the ones given in the input file
    """,
                                     usage='%(prog)s --input input_file.xlsx  --background bg.xlsx '
                                           '[--output an output folder] ')

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--input', dest='input', help="the input file for the tRNA program",
                        required=True)
    requiredNamed.add_argument('--background', dest='background', help="the background file of the tRNA program",
                        required=True)
    parser.add_argument('--output', dest="output", help="an output folder, default : your working directory",
                        default=".")
    args = parser.parse_args()


# launching the program
if __name__ == "__main__":
    main()

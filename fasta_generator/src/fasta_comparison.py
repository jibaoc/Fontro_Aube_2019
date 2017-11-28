import subprocess
import argparse
import os
import pandas as pd


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
                           "--output", output, "--filename", name_file],
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

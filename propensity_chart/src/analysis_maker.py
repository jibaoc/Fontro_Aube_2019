"""Script to make all the propensity graphics given a data folder."""

# Import
import argparse
import os
import subprocess


def find_all_folder(data_folder):
    """
    Find all sub-folder in data_folder.

    :param data_folder: (string) the path where the data folder are located
    :return: the list of sub-folders in data_folder
    """
    list_of_file = subprocess.check_output(["ls", data_folder]).decode("utf-8")
    list_of_file = list_of_file.strip("\n").split("\n")
    return list_of_file


def return_up_and_down_file(data_folder, sub_folder):
    """
    Return the file containing up and down-regulated peptide in data_folder.

    :param sub_folder: (string) a sub-folder within data_folder
    :param data_folder: (string) the path where the data folder are located
    :return: (list of 2 string) the up or down peptide file
    """
    res_file = [data_folder + sub_folder + "/query_result_up.xlsx",
                data_folder + sub_folder + "/query_result_down.xlsx"]
    return res_file


def chart_launcher(up_file, down_file, output):
    """
    Launch the propensity_chart_maker.py file.

    :param up_file: (string) the file containing up-regulated exons encoded peptides
    :param down_file: (string) the file containing down-regulated exons encoded peptides
    :param output: (string) folder where the graph will be created
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    subprocess.check_call(["python2", file_dir + "/propensity_chart_maker.py", "--up",
                           up_file, "--down", down_file, "--output", output],
                          stderr=subprocess.STDOUT)


def result_creator(data_folder, result_folder):
    """
    Create the propensity chart for every project (subfolder) in data_folder.

    :param result_folder: (path where the result will be created
    :param data_folder: (string) the path where the data folder are located
    """
    sub_folders = find_all_folder(data_folder)
    for folder in sub_folders:
        print("Creating figures of project " + str(folder))
        output_folder = result_folder + folder + "/"
        os.mkdir(output_folder)
        files = return_up_and_down_file(data_folder, folder)
        chart_launcher(files[0], files[1], output_folder)


def launcher():
    """Function that contains a parser to launch the program."""
    # description on how to use the program
    desc = """From a folder with the following structure :
    folder/
     |
     folder1/
     |  |
     |  query_result_up.xlsx
     |  query_result_down.xlsx
     folder2/
     |  |
     |  query_result_up.xlsx
     |  query_result_down.xlsx
     ...
     foldern/

    create in a given output_folder all the sub-folders in your given folder (i.e.
    folder1, folder2..., foldern here) and create in those sub_folders the propensity charts
    for the 2 sets of peptides (contained in query_result_up.xlsx, query_result_down.xlsx)
    given in the sequence sheet of query_result_up/down.xlsx
    """
    format_arg = argparse.RawDescriptionHelpFormatter
    usage = '%(prog)s --data_folder a_data_folder --result_folder your_result_folder'
    parser = argparse.ArgumentParser(formatter_class=format_arg,
                                     description=desc,
                                     usage=usage)
    # Arguments for the parser
    required_named = parser.add_argument_group('required arguments')
    required_named.add_argument('--data_folder', dest='data_folder', required=True,
                                help="the folder containing the data result")
    required_named.add_argument('--result_folder', dest='result_folder', required=True,
                                help="the file where the output will be created")

    args = parser.parse_args()  # parsing arguments
    # checking if the output and the result folder exists
    if not os.path.isdir(args.data_folder):
        print("ERROR : your data folder doesn't exist !")
        print("Exiting...")
        exit(1)
    if not os.path.isdir(args.result_folder):
        print("ERROR : your result_folder doesn't exist !")
        print("Exiting...")
        exit(1)

    result_creator(args.data_folder, args.result_folder)


if __name__ == "__main__":
    launcher()

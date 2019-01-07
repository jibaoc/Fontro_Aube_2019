#!/usr/bin/python3

"""Description:
    The goal of this script is to launch ``enrichment_of_exon_fixating_SF.py`` \
    for every file in ``tRNA_input_FarLine`` folder
"""

import os
import subprocess
import argparse


def main(trna_launcher, folder_input, folder_bed, output):
    """

    :param trna_launcher: (string) file corresponding to the tRNA launcher
    :param folder_input: (string) folder containing the tRNA input file
    :param folder_bed: (string) folder containing the bed files
    :param output:  (string) path where the output will be created
    """
    cur_folder = os.path.realpath(os.path.dirname(__file__))
    input_files = sorted(os.listdir(folder_input))
    input_files = [folder_input + my_file for my_file in input_files]
    bed_files = sorted(os.listdir(folder_bed))
    bed_files = [folder_bed + my_file for my_file in bed_files]
    for my_input in input_files:
        print("Working on %s" % my_input)
        name_file = os.path.basename(my_input).split(".")[0]
        my_bed = None
        for bed_file in bed_files:
            if name_file.split("_")[0] == os.path.basename(bed_file).split(".")[0]:
                my_bed = bed_file
                break
        output_folder = output + name_file
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)

        cmd = "python3 %s/enrichment_of_exon_fixating_SF.py  --output %s --name %s --clip_bed %s  --trna_input %s \
        --trna_launcher %s" % (cur_folder, output_folder, name_file, my_bed, my_input, trna_launcher)
        print(cmd)
        subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Launch the script ``enrichment_of_exon_fixating_sf`` multiple time

    """)
    # Arguments for the parser
    required_args = parser.add_argument_group("Required argument")


    parser.add_argument('--output', dest='output',
                        help="""path where the result will be created - default : current directory""",
                        default=".")
    required_args.add_argument('--clip_folder', dest='clip_folder',
                               help="The bed folder containing bed file",
                               required=True)

    required_args.add_argument('--trna_input', dest='trna_input',
                               help="""the name of the tRNA input folder selected""",
                               required=True)
    required_args.add_argument('--trna_launcher', dest='trna_launcher',
                               help="""file corresponding the tRNA launcher""",
                               required=True)
    args = parser.parse_args()  # parsing arguments

    # Defining global parameters
    if args.output[-1] != "/":
        args.output += "/"
    main(args.trna_launcher, args.trna_input, args.clip_folder, args.output)


if __name__ == "__main__":
    launcher()


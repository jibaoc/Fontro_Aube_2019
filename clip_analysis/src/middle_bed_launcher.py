#!/usr/bin/env python3

"""
Description:
    This script will launch ``middle_bed_enrichment`` for every bed in a given folder
"""


import os
import subprocess
import argparse


def main(trna_launcher, folder_bed, fasterdb_bed, output):
    """

    :param trna_launcher: (string) file corresponding to the tRNA launcher
    :param folder_bed: (string) folder containing the bed files
    :param output:  (string) path where the output will be created
    """
    cur_folder = os.path.realpath(os.path.dirname(__file__))
    bed_files = sorted(os.listdir(folder_bed))
    bed_files = [folder_bed + my_file for my_file in bed_files]
    for my_bed in bed_files:
        print("Working on %s" % my_bed)
        name_file = os.path.basename(my_bed).split(".")[0]
        output_folder = output + name_file
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)

        cmd = "python3 %s/middle_bed_enrichment.py  --output %s --name %s --clip_bed %s  --fasterdb_bed %s \
        --trna_launcher %s" % (cur_folder, output_folder, name_file, my_bed, fasterdb_bed, trna_launcher)
        if "SRSF3" in my_bed:
            cmd += " --overlap 2"
        if "SRSF1" in my_bed:
            cmd += " --overlap 5"
        print(cmd)
        subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Launch the script ``middle_bed_enrichment`` multiple time

    """)
    # Arguments for the parser
    required_args = parser.add_argument_group("Required argument")


    parser.add_argument('--output', dest='output',
                        help="""path where the result will be created - default : current directory""",
                        default=".")
    required_args.add_argument('--clip_folder', dest='clip_folder',
                               help="The bed folder containing bed file",
                               required=True)

    required_args.add_argument('--fasterdb_bed', dest='fasterdb_bed',
                               help="""the bed containing all fasterDB exons""",
                               required=True)
    required_args.add_argument('--trna_launcher', dest='trna_launcher',
                               help="""file corresponding the tRNA launcher""",
                               required=True)
    args = parser.parse_args()  # parsing arguments

    # Defining global parameters
    if args.output[-1] != "/":
        args.output += "/"
    main(args.trna_launcher, args.clip_folder, args.fasterdb_bed, args.output)


if __name__ == "__main__":
    launcher()


import pandas as pd
import argparse


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

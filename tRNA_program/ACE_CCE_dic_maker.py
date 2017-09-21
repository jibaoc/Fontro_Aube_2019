from main_project import *


import argparse

def main(input_file, input_file2, out_path, complete_exon, correspondence):
    full_exon_list, input_content = execute_request(input_file, out_path, complete_exon, correspondence, True)
    res = full_exon_list.codon_rich_exon_analyser()
    fichier = open(out_path + "ACE_rich.py","w")
    fichier.write("dic_ACE = " + str(res))
    fichier.close()
    poor_ace = full_exon_list.codon_poor_exon_analyser()
    fichier = open(out_path + "ACE_poor.py","w")
    fichier.write("dic_ACE = " + str(poor_ace))
    fichier.close()
    last_nt_ace = full_exon_list.codon_last_nt_proportion()
    fichier = open(out_path + "ACE_last_nt.py","w")
    fichier.write("dic_ACE = " + str(last_nt_ace))
    fichier.close()
    only_ace = full_exon_list.codon_plus_exon_analyser()
    fichier = open(out_path + "ACE_only.py","w")
    fichier.write("dic_ACE = " + str(only_ace))
    fichier.close()

    full_exon_list, input_content = execute_request(input_file2, out_path, complete_exon, correspondence, True)

    fichier = open(out_path + "CCE_rich.py", "w")
    cce = full_exon_list.codon_rich_exon_analyser()
    fichier.write("dic_CCE = " + str(cce))
    fichier.close()

    fichier = open(out_path + "CCE_poor.py", "w")
    poor_cce = full_exon_list.codon_poor_exon_analyser()
    fichier.write("dic_CCE = " + str(poor_cce))
    fichier.close()

    fichier = open(out_path + "CCE_last_nt.py", "w")
    last_nt_cce = full_exon_list.codon_last_nt_proportion()
    fichier.write("dic_CCE = " + str(last_nt_cce))
    fichier.close()

    fichier = open(out_path + "CCE_only.py", "w")
    only_cce = full_exon_list.codon_plus_exon_analyser()
    fichier.write("dic_CCE = " + str(only_cce))
    fichier.close()


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=""" display the dictionary of enriched codon (in each nucleotide)
                                      within the all set of exons given in input
    """,
                                     usage='%(prog)s --input input_file.txt --input2 input_file.txt '
                                           '[--output an output folder] ')
    # Arguments for the parser

    parser.add_argument('--input', dest='input',
                        help="""your request containing an exon name, its chromosome number and its
                        chromosomal coordinates""", default=None)
    parser.add_argument('--input2', dest='input2',
                        help="""your request containing an exon name, its chromosome number and its
                            chromosomal coordinates""", default=None)

    parser.add_argument('--output', dest='output', help="An output folder",
                        default=".")


    args = parser.parse_args()  # parsing arguments
    if args.output[-1] != "/":
        args.output += "/"
    main(args.input, args.input2, args.output, True, True)

launcher()
#!/usr/bin/env python3

"""
Description:
    The goal of this script is to select the exon for which have have an SRSF1 peak falling exactly within the exons. \
    We will then intersect this set of exons with the exons found down regulated by this splicing factor with Farline. \
    The we will se if the intersect list of exon present some particular enrichment.
"""

import os
import subprocess
import argparse


def create_bed_from_tran_input(input_file, output):
    """
    Transform a tRNA input in a bed file
    :param input_file: (string) a tRNA program input_file
    :param output: (string) path where the bed will bed created
    :return: (string) the name of the bed file created
    """
    new_name = os.path.basename(input_file).split(".")[0]
    new_file = output + new_name
    with open(input_file, "r") as infile:
        with open(new_file, "w") as outfile:
            for line in infile:
                line = line.replace("\n", "").split("\t")
                new_line = [line[1], str(int(line[2]) - 1), line[3],  line[0]]
                outfile.write("\t".join(new_line) + "\n")
    return new_file


def intersectbed(farline_bed, clip_bed, output):
    """
    Make the intersection using 2 bed.

    :param farline_bed: (string) the bed file build from a tRNA program input file.
    :param clip_bed: (string)  the bed file corresponding to the bed analysis
    :param output: (string) folder where the output file will be created
    :return: (string) the name of the intersect bed file.
    """
    outfile = output + "intersect.bed"
    cmd = "intersectBed -a %s -b %s -wa -u -F 1 > %s" % (farline_bed, clip_bed, outfile)
    print(cmd)
    subprocess.check_call(cmd, shell=True,
                          stderr=subprocess.STDOUT)
    return outfile


def intersect2input(intersect_bed, output):
    """
    Transform a bed file into an
    :param intersect_bed: (string) a bed file that just have been intersected
    :param output: (string) path where the tRNA input will bed created
    :return: (string) the name of the tRNA input
    """
    name_input = output + "input_file.txt"
    with open(intersect_bed, "r") as bedfile:
        with open(name_input, "w") as outfile:
            for line in bedfile:
                line = line.replace("\n", "").split("\t")
                new_line = [line[3], line[0], str(int(line[1]) + 1), line[2]]
                outfile.write("\t".join(new_line) + "\n")
    return name_input


def clip_bed_adapter(clip_bed, output):
    """
    Adapt the format of the bed obtain from clip
    :param clip_bed: (string) a bed file obtain after a clip experiment
    :param output: (string) path where the adapted file will be created
    :return: (string) the name of the clip adapted
    """
    new_file = output + os.path.basename(clip_bed).replace(".bed", "_1based.bed")
    with open(clip_bed, "r") as bedfile:
        with open(new_file, "w") as outfile:
            for line in bedfile:
                if line[0] != "#":
                    line = line.replace("\n", "").replace("#", "-")
                    line = line.split("\t")
                    line[0] = line[0].replace("chr", "")
                    line[1] = str(int(line[1]) + 1)
                    outfile.write("\t".join(line[0:6]) + "\n")
    return new_file


def enrichment_analysis(input_file, output, exon_type, launcher_file, name_analysis):
    """
    Launch the enrichment analysis.

    :param input_file: (string) input file to launch the tRNA program
    :param output:  (string) path where the result will be created
    :param exon_type: (string) the type of control exon used
    :param launcher_file: (string) location of the file that will launch the tRNA program
    :param name_analysis: (string) project name
    :return:
    """

    cmd1 = "python2 " + launcher_file + " \
            --input " +  input_file + " \
            --exon_type " + exon_type + " \
            --thread 2 \
            --output " + output + " \
            --enrichment True \
            --set_number 10000 \
            --alt_figure no \
            --summary no \
            --set_name \"" + str(name_analysis) + "\""
    cmd1 = cmd1.replace("            ", "")
    print(cmd1)
    cmd1 = cmd1.split(" ")
    subprocess.check_call(cmd1, stderr=subprocess.STDOUT)


def main(output, clip_bed, trna_input, trna_launcher, exon_type, name_analysis):
    """
    Create a tRNA input_file by intersecting a bed composed of the exon found down-regulated by a splicing factor \
    by the clip found within an exon for this splicing factor. Then an enrichment analysis is performed on \
     the exons in this intersection

    :param output: (string) path where the output will be created
    :param clip_bed: (string) a file corresponding to a bed file produced from clip analysis
    :param trna_input: (string) file containing a input of the tRNA program
    :param trna_launcher: (string) path to the tRNA launcher
    :param exon_type: (string) the control exon type used for the enrichment
    :param name_analysis: (string) the name of the analysis
    """
    bed_farline = create_bed_from_tran_input(trna_input, output)
    new_clip = clip_bed_adapter(clip_bed, output)
    bed_intersect = intersectbed(bed_farline, new_clip, output)
    final_trna_input = intersect2input(bed_intersect, output)
    new_output = output + "Exon_analysis_%s" % exon_type + "/"
    if not os.path.isdir(new_output):
        os.mkdir(new_output)
    enrichment_analysis(final_trna_input, new_output, exon_type, trna_launcher, name_analysis)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    From 2 files : one corresponding to a tRNA input file containing some exon
    found differentially splicing in a farline analysis and the other is a bed 
    file optained from clip data intersect them to obtain an new input for the
    tRNA program and launch the anlysis from this new input.

    """)
    # Arguments for the parser
    required_args = parser.add_argument_group("Required argument")

    parser.add_argument('--exon_type', dest='exon_type',
                        help="""the exon type of the control exon used for the enrichment analysis""",
                        default="CCE")
    parser.add_argument('--output', dest='output',
                        help="""path where the result will be created - default : current directory""",
                        default=".")
    parser.add_argument('--name', dest='name',
                        help="""the name of the analysis""",
                        default=".")
    required_args.add_argument('--clip_bed', dest='clip_bed',
                               help="The bed file obtained from clip data",
                               required=True)

    required_args.add_argument('--trna_input', dest='trna_input',
                               help="""the name of the tRNA input selected""",
                               required=True)
    required_args.add_argument('--trna_launcher', dest='trna_launcher',
                               help="""file corresponding the tRNA launcher""",
                               required=True)
    args = parser.parse_args()  # parsing arguments

    # Defining global parameters
    if args.output[-1] != "/":
        args.output += "/"
    main(args.output, args.clip_bed, args.trna_input, args.trna_launcher, args.exon_type, args.name)


if __name__ == "__main__":
    launcher()


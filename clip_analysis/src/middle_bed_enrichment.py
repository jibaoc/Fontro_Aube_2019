#!/usr/bin/env python3

"""
Description:
    The goal of this script is from a bed containing the peaks for a specific factors concatenated from different \
    clip experiment, to see if those peaks (falling in an exon) have a particular enrichment of features \
    using the tRNA program
"""

import enrichment_of_exon_fixating_SF
import os
import math
import subprocess
import argparse


def get_middle_bed(intersect_bed, output, wanted_size):
    """
    Modify the ``intersect_bed`` file to get only the middle part of the bed file with the size \
    ``wanted_size``

    :param intersect_bed:  (string) path to a bed file containing the peaks of a particular factor \
    falling exactly within an exon.
    :param output: (string) path where the bed will be created.
    :param wanted_size:  (int) the size of the new wanted interval
    :return: (string) the middle bed name
    """
    out_file = output + os.path.basename(intersect_bed).replace(".bed", "_middle_seq.bed")
    with open(intersect_bed, "r") as inbed:
        with open(out_file, "w") as outbed:
            for line in inbed:
                line = line.replace("\n", "")
                line = line.split("\t")
                line[1] = int(line[1])
                line[2] = int(line[2])
                if line[2] - line[1] + 1 > wanted_size :
                    middle = math.ceil((line[2] + line[1]) / 2)
                    start = middle - math.floor(wanted_size / 2)
                    stop = middle + math.ceil(wanted_size / 2) - 1
                    if start < 1:
                        print("Warning on line : %s, start below 1" % str(line))
                    if stop > line[2]:
                        print("Warning on line : %s, stop greater than line[2] " % str(line))
                    line[1] = start
                    line[2] = stop
                    outbed.write("\t".join(list(map(str, line))) + "\n")
                elif line[2] - line[1] + 1 == wanted_size:
                    outbed.write("\t".join(list(map(str, line))) + "\n")
    return out_file


def remove_duplicate(middle_bed, output):
    """
    Remove the duplicate sequences in ``middle_bed``.

    :param middle_bed: (string) a bed containing middle sequence
    :param output: (string) path where the output will be created
    :return: (string) the bed files without duplicate
    """
    # As the tRNA program only need coordinates and not the strand we don't take it into account to removeduplicate
    new_file = output + os.path.basename(middle_bed).replace(".bed", "_noDup.bed")
    subprocess.check_call("sort -k1,1 -k2,2n -k3,3n -u %s > %s" % (middle_bed, new_file), shell=True,
                          stderr=subprocess.STDOUT)
    return new_file


def enrichment_analysis(input_file, output, exon_type, launcher_file, name_analysis, size_wanted):
    """
    Launch the enrichment analysis.

    :param input_file: (string) input file to launch the tRNA program
    :param output:  (string) path where the result will be created
    :param exon_type: (string) the type of control exon used
    :param launcher_file: (string) location of the file that will launch the tRNA program
    :param name_analysis: (string) project name
    :param size_wanted: (int) the size of the new wanted interval (used to select control sequence with the same size)
    :return:
    """
    cmd1 = "python2 " + launcher_file + " \
            --input " + input_file + " \
            --exon_type " + exon_type + " \
            --complete_exon False \
            --duplication True \
            --penalty_size 3 \
            --size_control False \
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


def intersectbed(bed_clip, fasterdb_bed, output):
    """
    Make the intersection using 2 bed.

    :param bed_clip: (string) the bed file corresponding to the clip bed analysis
    :param fasterdb_bed: (string)  the bed file corresponding to all fasterdb exons
    :param output: (string) folder where the output file will be created
    :return: (string) the name of the intersect bed file.
    """
    outfile = output + "intersect.bed"
    cmd = "intersectBed -a %s -b %s -wa -u -f 1 > %s" % (bed_clip, fasterdb_bed, outfile)
    print(cmd)
    subprocess.check_call(cmd, shell=True,
                          stderr=subprocess.STDOUT)
    return outfile


def bed_sort(bed_file, output):
    """
    sort the bed ``bed_file``.

    :param bed_file: (string) a bed file
    :param output: (string) path where the bed will be created
    :return: (string) the name of the sorted bed
    """
    new_file = output + os.path.basename(bed_file).replace(".bed", ".sort.bed")
    cmd = "sort -k1,1 -k2,2n %s > %s" % (bed_file, new_file)
    subprocess.check_call(cmd, shell=True,
                          stderr=subprocess.STDOUT)
    return new_file


def merge_bed(bed_file, output):
    """
    Merge overlapping feature and only keep those that merge 3 features.

    :param bed_file:  (string) a bed file
    :param output: (string) path where the bed will be created
    :return: (string) the name of the sorted bed
    """
    new_file = output + os.path.basename(bed_file).replace(".bed", ".merged.bed")
    cmd = "mergeBed -i %s -s -c 1 -o count | awk '{if ($5 > 2) print $0}' > %s" % (bed_file, new_file)
    subprocess.check_call(cmd, shell=True,
                          stderr=subprocess.STDOUT)
    return new_file



def main(output, clip_bed, fasterdb_bed, trna_launcher, exon_type, name_analysis, wanted_size):
    """
    Create a tRNA input_file by intersecting a bed composed of the exon found down-regulated by a splicing factor \
    by the clip found within an exon for this splicing factor. Then an enrichment analysis is performed on \
     the exons in this intersection

    :param output: (string) path where the output will be created
    :param clip_bed: (string) a file corresponding to a bed file produced from clip analysis
    :param fasterdb_bed: (string) file containing a input of the tRNA program
    :param trna_launcher: (string) path to the tRNA launcher
    :param exon_type: (string) the control exon type used for the enrichment
    :param name_analysis: (string) the name of the analysis
    :param wanted_size: (string) the wanted size
    """
    sorted_bed = bed_sort(clip_bed, output)
    adapted_clip = enrichment_of_exon_fixating_SF.clip_bed_adapter(sorted_bed, output)
    bed_intersect = intersectbed(adapted_clip, fasterdb_bed, output)
    middle_bed = get_middle_bed(bed_intersect, output, wanted_size)
    sorted_bed = bed_sort(middle_bed, output)
    merged_bed = merge_bed(sorted_bed, output)
    nodup_bed = remove_duplicate(merged_bed, output)  # remove potential feature with the same coordinates on different strand
    final_trna_input = enrichment_of_exon_fixating_SF.intersect2input(nodup_bed, output)
    new_output = output + "Exon_analysis_%s" % exon_type + "/"
    if not os.path.isdir(new_output):
        os.mkdir(new_output)
    enrichment_analysis(final_trna_input, new_output, exon_type, trna_launcher, name_analysis, wanted_size)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Program that will test then enrichment of a set of peaks (in a bed file and located in exons), 
    obtained from clip experiment, in nucleotide, codon amino acid...

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
                        default="")
    parser.add_argument('--size', dest='size',
                        help="""size of the peaks selected""",
                        default=20)
    required_args.add_argument('--clip_bed', dest='clip_bed',
                               help="The bed file obtained from clip data",
                               required=True)

    required_args.add_argument('--fasterdb_bed', dest='fasterdb_bed',
                               help="""bed files contaiing fasterdb (i.e hg19) exons""",
                               required=True)
    required_args.add_argument('--trna_launcher', dest='trna_launcher',
                               help="""file corresponding the tRNA launcher""",
                               required=True)

    args = parser.parse_args()  # parsing arguments

    # Defining global parameters
    if args.output[-1] != "/":
        args.output += "/"
    try:
        args.size = int(args.size)
    except ValueError:
        parser.error("Error : wrong value for argument size")
    main(args.output, args.clip_bed, args.fasterdb_bed, args.trna_launcher, args.exon_type, args.name, args.size)


if __name__ == "__main__":
    launcher()


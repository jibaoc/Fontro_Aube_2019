#!/usr/bin/python3


"""
Description:

The goal of this script is to:
    1. Create a lot of random fasta enriched in a particular unit (*nucleotide*, *di-nucleotide* or *feature*) and calculates the mean frequency of this unit in the fasta file
    2. Create a lot of random fasta impoverished in **the same unit** (*nucleotide*, *di-nucleotide* or *feature*) and calculates the mean frequency of this unit in the fasta file

And then compare their frequency for this unit for each couple of random fasta enriched and impoverished \
(for this unit). To compare their frequency, the relative frequency is computed

The relative frequency is calculated as follow:

.. math::

  F_{relative} = \frac{F_{interest} - F_{control}}{F_{control}}

Where:
  * :math:`F_{relative}` is the relative frequency of a unit :math:`F`
  * :math:`F_{interest}` is the frequency of :math:`F` in the interest set of exons
  * :math:`F_{control}` is the frequency of :math:`F` in the control sets of exons

We then calculate the mean and the standard deviation of:
    1. The list of frequencies in **unit** obtained by creating random fasta files enriched in this **unit**
    2. The list of frequencies in **unit** obtained by creating random fasta files impoverished in this **unit**
    3. The list of relative frequencies between random fasta files enriched and impoverished in this **unit**

"""

import subprocess
from Bio import SeqIO
import numpy as np
import argparse
import os
import sys

feature_dic = {
    "Small#1": ["A", "C", "D", "G", "N", "P", "S", "T", "V"],
    "Small#2": ["A", "C", "D", "G", "N", "P", "S", "T"],
    "Large" : ["F", "I", "K", "L", "M", "R", "W", "Y"],
    "Disorder_promoting#1": ["A", "E", "G", "K", "P", "Q", "R", "S"],
    "Order_promoting#1": ["C", "F", "I", "L", "N", "W", "V", "Y"],
    "Disorder_promoting#2": ["A", "E", "G", "K", "P", "Q", "S"],
    "Order_promoting#2": ["C", "F", "H", "I", "L", "M", "N", "W", "V", "Y"],
    "Polar_uncharged#1": ["C", "N", "Q", "S", "T", "Y"],
    "Polar_uncharged#2": ["N", "Q", "S", "T", "Y"],
    "Charged": ["R", "H", "K", "D", "E"],
    "Hydrophilic#1": ["D", "E", "K", "N", "Q", "R"],
    "Hydrophobic#1": ["A", "C", "F", "I", "L", "M", "V"],
    "Hydrophilic#2": ["D", "E", "H", "K", "N", "Q", "R", "S", "T"],
    "Hydrophobic#2": ["A", "C", "F", "I", "L", "M", "P", "V", "W", "Y"],
    "Hydroxylic": ["S", "T", "Y"],
    "Negatively_charged": ["D", "E"],
    "Positively_charged": ["R", "H", "K"],
}



def fasta_generation_nt_dnt(nt_dnt, freq, filename, output, iscub):
    """
    Launch the genesis of random fasta with fixed nucleotide or dinucleotide frequency. \
    The random sequences can by computed from codon usage of CCE exons or from real CCE exons. \
    The generated sequences are stored in a fasta file located in the ``output`` directory

    :param nt_dnt: (string) the nucleotide or dinucleotide for which we want to have a fixed frequency in random \
    sequences
    :param freq: (float) the frequency of the ``nt_dnt``
    :param filename: (string) the name of the file to create
    :param output: (string) the path where the random fasta will be created
    :param iscub: (boolean)
        * **True** if the random sequences are computed from the CCE codon usage bias or \
        * **False** if the random sequences are computed from real CCE exons
    """
    cmd = ["python"]
    if not iscub:
        cmd += ["/home/nicolas/PycharmProjects/fasta_generator/src/fasta_generator_dinucleotide_from_real_exon.py"]
    else:
        cmd += ["/home/nicolas/PycharmProjects/fasta_generator/src/fasta_generator.py"]
    cmd += ["--output", output]
    cmd += ["--nt_dnt", nt_dnt]
    cmd += ["--freq", str(freq)]
    cmd += ["--ctrl",  "CCE"]
    cmd += ["--filename", filename]
    subprocess.check_output(cmd, stderr=subprocess.STDOUT)


def fasta_generation_feature(feature, freq, output, iscub):
    """
    Launch the genesis of random fasta with fixed feature frequency. \
    The random sequences can by computed from codon usage of CCE exons or from real CCE exons. \
    The generated sequences are stored in a fasta file located in the ``output`` directory

    :param feature: (string) the feature for which we want to have a fixed frequency in random \
    sequences
    :param freq: (float) the frequency of the ``feature``
    :param output: (string) the path where the random fasta will be created
    :param iscub: (boolean)
        * **True** if the random sequences are computed from the CCE codon usage bias or \
        * **False** if the random sequences are computed from real CCE exons
    """
    cmd = ["python"]
    if not iscub:
        cmd += ["/home/nicolas/PycharmProjects/fasta_generator/src/fasta_generator_from_real_exons.py"]
    else:
        cmd += ["/home/nicolas/PycharmProjects/fasta_generator/src/fasta_reverse_generator.py"]
    cmd += ["--output", output]
    cmd += ["--ctrl", "CCE"]
    cmd += ["--feature", feature]
    cmd += ["--prop", str(freq)]
    subprocess.check_output(cmd, stderr=subprocess.STDOUT)


def translator(sequence):
    """
    Translate a nucleotide sequence. The ORF of the sequence starts at the first nucleotide and codons can't \
    be truncated.

    :param sequence: (string a nucleotide sequence
    :return: (string an amino acid sequence
    """
    codon2aminoacid = dict(TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L", ATT="I", ATC="I",
                           ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V", TCT="S", TCC="S", TCA="S", TCG="S",
                           CCT="P", CCC="P", CCA="P", CCG="P", ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A",
                           GCA="A", GCG="A", TAT="Y", TAC="Y", TAA="*", TAG="*", CAT="H", CAC="H", CAA="Q", CAG="Q",
                           AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E", TGT="C", TGC="C",
                           TGA="*", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R", AGT="S", AGC="S", AGA="R", AGG="R",
                           GGT="G", GGC="G", GGA="G", GGG="G")
    new_seq = []
    for i in range(0, len(sequence), 3):
        if i+3 <= len(sequence):
            new_seq.append(sequence[i:i+3])
    aa_seq = ""
    for codon in new_seq:
        aa_seq += codon2aminoacid[codon]
    return aa_seq


def get_feature_frequency(sequence, feature):
    """
    Get the feature frequency of a nucleotide sequence. It's ORF starts at the first nucleotide and codons can't \
    be truncated.

    :param sequence: (string) a nucleotide sequence
    :param feature: (string) the feature for which we want to compute the frequency in the sequence
    :return: (float) the frequency of the feature ``feature``
    """
    aa_seq = translator(sequence)
    count = 0.
    for aa in feature_dic[feature]:
            count += aa_seq.count(aa)
    return (float(count) / len(aa_seq)) * 100


def get_dinucleotide_frequency(sequence, dnt):
    """
    Get the di-nucleotide frequency of a sequence

    :param sequence: (string) a nucleotide sequence
    :param dnt: (string) a di-nucleotide for which we want to calculate the frequency
    :return: (float) the frequency of the di-nucleotide ``dnt`` in ``sequence``
    """
    count = 0.
    for i in range(len(sequence) - 1):
        if sequence[i:i+2] == dnt:
            count += 1
    return (count / (len(sequence) - 1)) * 100


def get_nucleotide_frequency(sequence, nt):
    """
    Get the nucleotide frequency of a nucleotide sequence

    :param sequence: (string) a nucleotide sequence
    :param nt: (string) a nucleotide for which we want to compute teh frequency
    :return: (float) the frequency of the nucleotide ``nt`` in ``sequence``
    """
    iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
             'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}
    count = 0.
    if nt in ["A", "C", "G", "T"]:
        for i in range(len(sequence)):
            if sequence[i] == nt:
                count += 1
    else:
        for i in range(len(sequence)):
            if sequence[i] in iupac[nt]:
                count += 1
    return (count / len(sequence)) * 100


def get_mean_frequency_in_fasta_file(fasta_file, type_unit, unit):
    """
    Give the mean frequency of the ``type_unit`` ``unit`` in the all set of sequence \
    contains in the fasta files.

    :param fasta_file: (string) path to a fasta files
    :param type_unit: (string) the unit type we want to study : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :return: the mean frequency of ``unit`` of all the sequence in ``fasta_file``
    """
    mean_list = []
    while len(mean_list) == 0:
        mean_list = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq)
            if type_unit == "feature":
                mean_list.append(get_feature_frequency(seq, unit))
            if type_unit == "dnt":
                mean_list.append(get_dinucleotide_frequency(seq, unit))
            if type_unit == "nt":
                mean_list.append(get_nucleotide_frequency(seq, unit))
    return np.mean(mean_list)


def get_mean_frequency_of_multiple_fasta(type_unit, unit, freq, iteration, output, iscub):
    """
    Create ``iteration`` fasta files having a frequency of ``freq`` of ``unit`` and gets the frequency of \
    ``unit`` in all the fasta files generated

    :param type_unit: (string) the unit type we want to study : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :param freq: (float) the frequency of ``unit`` in the fasta file that will be produced
    :param iteration: (int) the number of fasta files we want to create
    :param output: (string) path where the fasta_file will be created
    :param iscub: (bollean) True if we want to create random fasta respecting CCE codon bias usage \
    False if we want to mutate fasta sequence
    :return: (list of ``itaration`` float value) list of the mean of ``unit`` in the ``iteration`` fasta \
    files generated.
    """
    list_mean_fasta = []
    filename_fasta = "CCE_" + unit + "_" + str(freq)
    for i in range(iteration):
        sys.stdout.write("\rprogress: " + str(i + 1) + " / " + str(iteration))
        if i == iteration - 1:
            sys.stdout.write("\n")
        sys.stdout.flush()
        if type_unit == "feature":
            fasta_generation_feature(unit, freq, output, iscub)
        else:
            fasta_generation_nt_dnt(unit, freq, filename_fasta, output, iscub)
        list_mean_fasta.append(get_mean_frequency_in_fasta_file(output + filename_fasta + ".fasta", type_unit, unit))
    subprocess.check_call(["rm", output + filename_fasta + ".fasta"], stderr=subprocess.STDOUT)
    return list_mean_fasta


def get_relative_freq_values(list_high, list_low):
    """
    Get the relative frequencies of fasta_files having up and low frequency for a particular unit.

    :param list_high: (list of float) list of mean frequencies in a particular unit in a list of fasta file enriched \
    for this unit
    :param list_low: (list of float) list of mean frequencies in a particular unit in a list of fasta file impoverished \
    for this unit
    :return: (list of float) relative frequencies for each value of list_high and list low. i.e \

    .. code-block:: bash

        relative\_freq[i] = \frac{list_high[i] - list_low[i]}{list_low[i]}


    """
    relative_freq = []
    for i in range(len(list_high)):
        relative_freq.append((list_high[i] - list_low[i]) / list_low[i])
    return relative_freq


def write_tsv_file(unit_type, unit, freq_high, freq_low, list_high_freq, list_low_freq, iteration, output, iscub):
    """
    Create a tsv file.

    :param unit_type: (string) the unit type we want to study : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :param freq_high: (float) the frequency of ``unit`` for the fasta file having an high content of the unit of \
    interest
        :param freq_low: (float) the frequency of ``unit`` for the fasta file having a low content of the unit of \
    interest
    :param list_high_freq: (list of float) list of mean frequencies in a particular unit in a list of fasta file enriched \
    for this unit
    :param list_low_freq: (list of float) list of mean frequencies in a particular unit in a list of fasta file impoverished \
    for this unit
    :param iteration: (int) the number of fasta files we want to create
    :param output: (string) path where the fasta_file will be created
    :param iscub: (bollean) True if we want to create random fasta respecting CCE codon bias usage \
    False if we want to mutate fasta sequence
    """
    if iscub:
        fname = "CUB"
    else:
        fname = "mutated"
    relative_freq = get_relative_freq_values(list_high_freq, list_low_freq)
    mean_rel_freq = np.mean(relative_freq)
    sd_rel_freq = np.std(relative_freq)

    mean_high = np.mean(list_high_freq)
    std_high = np.std(list_high_freq)

    mean_low = np.mean(list_low_freq)
    std_low = np.std(list_low_freq)
    unit = unit.replace("_", "-")
    with open(output + unit_type + "_" + unit + "_frequency_comparison_between_" + str(iteration) + "_" + fname + "_fasta_file-high:"+str(freq_high) + "_low:" + str(freq_low) + ".tsv", "w") as outfile:
        outfile.write("frequency_fasta" + fname + "_" + str(unit) + ":" + str(freq_high) + "\t")
        outfile.write("frequency_fasta" + fname + "_" + str(unit) + ":" + str(freq_low) + "\t")
        outfile.write("relative_frequency\n")
        for i in range(len(list_high_freq)):
            outfile.write(str(list_high_freq[i]) + "\t" + str(list_low_freq[i]) + "\t" + str(relative_freq[i]) + "\n")
        outfile.write(str(mean_high) + "\t" + str(mean_low) + "\t" + str(mean_rel_freq) + "\t" + "mean\n")
        outfile.write(str(std_high) + "\t" + str(std_low) + "\t" + str(sd_rel_freq) + "\t" + "std")


def main(type_unit, unit, freq_high, freq_low, iteration, output, iscub):
    """
    Create file file that the frequencies of ``unit`` of many (``ieration``) random fasta enriched for this unit to \
    the frequencies for that ``unit`` of many (``ieration``) random fasta imoverished for this unit

    :param type_unit: (string) the unit type we want to study : it can be *dnt* or *feature or *nt*
    :param unit: (string) the unit for which we want to calculate the median frequency
    :param freq_high: (float) the frequency of ``unit`` for the fasta file having an high content of the unit of \
    interest
        :param freq_low: (float) the frequency of ``unit`` for the fasta file having a low content of the unit of \
    interest
    :param iteration: (int) the number of fasta files we want to create
    :param output: (string) path where the fasta_file will be created
    :param iscub: (bollean) True if we want to create random fasta respecting CCE codon bias usage \
    False if we want to mutate fasta sequence
    """
    list_high_freq = get_mean_frequency_of_multiple_fasta(type_unit, unit, freq_high, iteration, output, iscub)
    list_low_freq = get_mean_frequency_of_multiple_fasta(type_unit, unit, freq_low, iteration, output, iscub)
    write_tsv_file(type_unit, unit, freq_high, freq_low, list_high_freq, list_low_freq, iteration, output, iscub)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
The goal of this script is to:
     1. Create a lot of random fasta enriched in a particular unit (*nucleotide*, *di-nucleotide* or *feature*)
     2. Create a lot of random fasta impoverished in **the same unit** (*nucleotide*, *di-nucleotide* or *feature*)

And then compare their frequency for this unit for each couple of random fasta enriched and impoverished \
(for this unit). To compare their frequency, the relative frequency is computed

The relative frequency is calculated as follow:

  F_{relative} = frac{F_{interest} - F_{control}}{F_{control}}

Where:
  * :math:`F_{relative}` is the relative frequency of a unit :math:`F`
  * :math:`F_{interest}` is the frequency of :math:`F` in the interest set of exons
  * :math:`F_{control}` is the frequency of :math:`F` in the control sets of exons

We then calculate the mean and the standard deviation of:
    1. The list of frequencies in **unit** obtained by creating random fasta files enriched in this **unit**
    2. The list of frequencies in **unit** obtained by creating random fasta files impoverished in this **unit**
    3. The list of relative frequencies between random fasta files enriched and impoverished in this **unit**
    """)
    # Arguments for the parser

    parser.add_argument('--iteration', dest='iteration', help="the number of random sequence to cerate",
                        default=100)

    parser.add_argument('--iscub', dest='iscub', help="True to create cub sequence False to create mutated sequence",
                        default=False)
    req_arg = parser.add_argument_group("required arguments")

    req_arg.add_argument('--type_unit', dest='type_unit', help="the type of unit of interest (nt, dnt or feature)",
                         required=True)
    req_arg.add_argument('--unit', dest='unit', help="the unit for which we want to create random fasta having a low "
                                                     "content of if and a high content of it",
                         required=True)
    req_arg.add_argument('--freq_high', dest='freq_high', help="the frequency of unit in your enriched fasta file "
                                                              "enriched for this unit",
                        required=True)
    req_arg.add_argument('--freq_low', dest='freq_low', help="the frequency of unit in your impoverished fasta file "
                                                              "enriched for this unit",
                        required=True)
    req_arg.add_argument('--output', dest='output', help="path where the result will be created",
                        required=True)



    args = parser.parse_args()

    if args.output[-1] != "/":
        args.output += "/"
    if not os.path.isdir(args.output):
        print("ERROR : the given output folder does not exist")
        exit(1)

    if args.iscub == "False":
        args.iscub = False

    if args.iscub == "True":
        args.iscub = True

    try:
        args.freq_high = float(args.freq_high)
    except ValueError:
        print("Wrong value for freq_high argument...")
        exit(1)

    try:
        args.freq_low = float(args.freq_low)
    except ValueError:
        print("Wrong value for freq_high argument...")
        exit(1)

    if args.type_unit not in ["nt", "dnt", "feature"]:
        print("Unrecognized type unit")
        exit(1)

    try:
        args.iteration = int(args.iteration)
    except ValueError:
        print("Wrong value for iteration argument...")
        exit(1)

    main(args.type_unit, args.unit, args.freq_high, args.freq_low, args.iteration, args.output, args.iscub)

if __name__ == "__main__":
    launcher()

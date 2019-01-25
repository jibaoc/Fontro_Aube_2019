#!/user/bin/python3

"""
The goal of this script is to make a barplot indicating the proportion of stretch different kind of stretches for
a nucleotide in the different sequence given in input.
It will allow to decide what type of stretches to use to better discriminate the input files
"""


# Set the environment :
from Bio import SeqIO
import pandas as pd
from math import sqrt
from ncephes import cprob
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import scipy.stats

codon2aminoAcid = dict(TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L", ATT="I", ATC="I",
                       ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V", TCT="S", TCC="S", TCA="S", TCG="S",
                       CCT="P", CCC="P", CCA="P", CCG="P", ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A",
                       GCA="A", GCG="A", TAT="Y", TAC="Y", TAA="", TAG="", CAT="H", CAC="H", CAA="Q", CAG="Q",
                       AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E", TGT="C", TGC="C",
                       TGA="", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R", AGT="S", AGC="S", AGA="R", AGG="R",
                       GGT="G", GGC="G", GGA="G", GGG="G")

feature_dic = {
    "Very-small": ["A", "C", "G", "S"],
    "Small#2": ["A", "C", "D", "G", "N", "P", "S", "T"],
    "Large" : ["F", "I", "K", "L", "M", "R", "W", "Y"],
    "Disorder-promoting#1": ["A", "E", "G", "K", "P", "Q", "R", "S"],
    "Order-promoting#1": ["C", "F", "I", "L", "N", "W", "V", "Y"],
    "Disorder-promoting#2": ["A", "E", "G", "K", "P", "Q", "S"],
    "Order-promoting#2": ["C", "F", "H", "I", "L", "M", "N", "W", "V", "Y"],
    "Polar-uncharged#1": ["C", "N", "Q", "S", "T", "Y"],
    "Polar-uncharged#2": ["N", "Q", "S", "T", "Y"],
    "Charged#1": ["R", "H", "K", "D", "E"],
    "Charged#2": ["R", "K", "D", "E"],
    "Hydrophilic#1": ["D", "E", "K", "N", "Q", "R"],
    "Hydrophobic#1": ["A", "C", "F", "I", "L", "M", "V"],
    "Neutral": ["G", "H", "P", "S", "T", "Y"],
    "Hydroxylic": ["S", "T", "Y"],
    "Negatively-charged": ["D", "E"],
    "Positively-charged#1": ["R", "H", "K"],
    "Positively-charged#2": ["R", "K"]
}

###
# Functions
###


def read_fasta_files(fasta_file, stretch_len):
    """
    :param fasta_file: (string) the name of a fasta file of interest
    :param stretch_len: (int) the len of the stretch of interest
    :return: the list of sequence in the fasta file, nt and amino
    """
    list_sequence = []
    nt_seq = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        seq_aa = translator(seq)
        if len(seq_aa) > stretch_len-1:
            list_sequence.append(seq_aa)
        if len(seq) > stretch_len-1:
            nt_seq.append(seq)
    return nt_seq, list_sequence


def read_sequence(excel_file, stretch_len):
    """
    :param stretch_len: (int) the len of the stretch of interest
    :param excel_file: (string) path to a query file given by the tRNA program
    :return: (2 list of strings) the list of cds sequence and the list of amino acid sequence
    """
    # opening the excel file
    xl = pd.ExcelFile(excel_file)
    df = "NA"
    # opening the sheet of interest
    for sheet in xl.sheet_names:
        if "sequence" == sheet:
            df = xl.parse(sheet)
    # if the sheet "exon_skipping*" doesn't exist, the faRLine file cannot be used so we stop the program
    if str(df) == "NA":
        print("the sheet names sequence wasn't found")
        print("exiting...")
        exit(1)

    cds = []
    pep = []
    for row in df.itertuples():
        if not isinstance(row.CDS_genomic_sequence, float):
            if len(row.CDS_genomic_sequence) > stretch_len-1:
                cds.append(row.CDS_genomic_sequence)
        if not isinstance(row.CDS_peptide_sequence, float):
            if len(row.CDS_peptide_sequence) > stretch_len-1:
                pep.append(row.CDS_peptide_sequence)
    return cds, pep


def translator(seq):
    """
    :param seq: (string) a nucleotide sequence
    :return: the translated sequence
    """
    res = ""
    for i in range(0, len(seq)-2, 3):
        res += codon2aminoAcid[seq[i:i+3]]
    return res


def stretch_finder_feature(sequence, feature, stretch_len, stretch_content):
    """
    :param sequence: (string) an amino acid sequences
    :param feature: (string) the name of the feature to return
    :param stretch_len: (int) the length of the stretch of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there ise a stretch in the sub-sequence
    :return: the number of stretch of the feature "feature" here
    """
    nb_stretch = 0
    for i in range(len(sequence) - stretch_len - 1):
        count = 0
        for letter in sequence[i:i+stretch_len]:
            if letter in feature_dic[feature]:
                count += 1
        if count >= stretch_content:
            nb_stretch += 1
    return nb_stretch


def stretch_finder_aa(sequence, aa, stretch_len, stretch_content):
    """
    :param sequence: (string) an amino acid sequences
    :param aa: (string) the name of the aa to return
    :param stretch_len: (int) the length of the stretch of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there ise a stretch in the sub-sequence
    :return: the number of stretch of the aa "aa" here
    """
    nb_stretch = 0
    for i in range(len(sequence) - stretch_len - 1):
        count = 0
        for letter in sequence[i:i+stretch_len]:
            if letter == aa:
                count += 1
        if count >= stretch_content:
            nb_stretch += 1
    return nb_stretch


def stretch_finder_nt(sequence, nt, stretch_len, stretch_content):
    """
    :param sequence: (string) an amino acid sequences
    :param nt: (string) the name of the nt to return
    :param stretch_len: (int) the length of the stretch of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there ise a stretch in the sub-sequence
    :return: the number of stretch of the nucleotide "nt" here
    """
    nb_stretch = 0
    if nt in ["A", "T", "G", "C"]:
        for i in range(len(sequence) - stretch_len - 1):
            count = 0
            for letter in sequence[i:i+stretch_len]:
                if letter == nt:
                    count += 1
            if count >= stretch_content:
                nb_stretch += 1
    else:
        iupac = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
                 'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T']}
        for i in range(len(sequence) - stretch_len - 1):
            count = 0
            for letter in sequence[i:i+stretch_len]:
                if letter in iupac[nt]:
                    count += 1
            if count >= stretch_content:
                nb_stretch += 1
    return nb_stretch


def get_stretch(list_of_sequence, unit, unit_type, stretch_len, stretch_content):
    """

    :param list_of_sequence: (list of string) list of amino acid sequences
    :param unit: (string) he name of the unit of interest
    :param stretch_len: (int) the length of the stretch of interest
    :param unit_type: (string) the name of the unit of interest (aa, feature, nt)
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there ise a stretch in the sub-sequence
    :return: (dictionary of int) the number of sequences having 0 to 10+ stretch in the list of sequences
    """
    st_dic = {i: 0 for i in range(0, 11)}
    for sequence in list_of_sequence:
        if unit_type == "feature":
            nb_stretch = stretch_finder_feature(sequence, unit, stretch_len, stretch_content)
        elif unit_type == "aa":
            nb_stretch = stretch_finder_aa(sequence, unit, stretch_len, stretch_content)
        else:
            nb_stretch = stretch_finder_nt(sequence, unit, stretch_len, stretch_content)
        if nb_stretch > 10:
            nb_stretch = 10
        st_dic[nb_stretch] += 1
    counter = 0
    for key in st_dic.keys():
        counter += st_dic[key]
    st_dic["all"] = counter
    return st_dic


def frequency_test(obs1, tot1, obs2, tot2):
    """

    :param obs1: (int) the count number of an amino acid X in the set of protein 1.
    :param tot1: (int) the total number of amino acids in the set of protein 1.
    :param obs2: (int) the count number of an amino acid X in the set of protein 2.
    :param tot2: (int) the total number of amino acids in the set of protein 2.
    :return: proportion test p-value
    """
    if obs1 == 0 and obs2 == 0:
        return float('nan')
    mean1 = float(obs1) / tot1
    mean2 = float(obs2) / tot2

    var1 = float(obs1) * (1 - mean1) * (1 - mean1) + (tot1 - obs1) * mean1 * mean1
    var2 = float(obs2) * (1 - mean2) * (1 - mean2) + (tot2 - obs2) * mean2 * mean2

    df = tot1 + tot2 - 2
    svar = (var1 + var2) / df
    t = (mean1-mean2) / sqrt(svar*(1.0/tot1 + 1.0/tot2))

    return cprob.incbet(0.5*df, 0.5, df/(df+t*t))


def r_ttest(x, y):
    """

    :param x: (list of float value)
    :param y: (list of float value)
    :return: R t test p value
    """
    import rpy2.robjects as robj
    import rpy2.robjects.vectors as v
    ttestmaker = robj.r("""
    function(x,y){
        test = t.test(x,y)
        return(test$p.value)
    }""")
    try:
        pval = ttestmaker(v.FloatVector(x), v.FloatVector(y))
        pval =  float(pval[0])
    except:
        pval = float("nan")
    return pval


def get_more(dic, num):
    """
    :param dic: dic of int) the number of exons having 0, 1, ..., 10 or more stretch in
    an exon set.
    :param num: (int) the number of exons having "num" stretch or more
    :return: the prop of exon having "num" stretch or more
    """
    count = 0
    for key in dic.keys():
        if key != "all":
            if key >= num:
                count += dic[key]
    return [float(count) / dic["all"], float(count)]


def barplot_maker(strecth_count_list, name_strecth_list, name_file, iunit_name, iunit, outfile, size_fasta, iteration):

    new_strecth_count_list = []
    for i in range(0, len(strecth_count_list), 3):
        new_strecth_count_list.append(strecth_count_list[i])
        new_strecth_count_list.append(strecth_count_list[i + 1])
    fig, ax = plt.subplots(figsize=(48. / 2.54 * 1.11, 27. * 0.91 / 2.54))
    # set the abscissa values for bars and tick_abscissa value for legend
    abscissa = []
    tick_abscissa = []
    i = 1
    for j in range(1, len(new_strecth_count_list) + 1):
        abscissa.append(i)
        if j % 2 == 0:
            i += 2
        elif j % 2 == 1:
            tick_abscissa.append(i + 0.5)
            i += 1
        else:
            i += 1

    colors = ["blue", "red"] * int(len(new_strecth_count_list) / 2)
    val = [np.mean(new_strecth_count_list[i]) for i in range(len(new_strecth_count_list))]
    std = [np.std(new_strecth_count_list[i]) for i in range(len(new_strecth_count_list))]
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray', linestyle='dashed', alpha=0.3)
    ymax = ax.get_ylim()[1]
    ax.set_ylim([0, ymax + ymax * 0.1])
    ax.bar(abscissa, val, width=0.8, color=colors, yerr=std)
    ax.set_xticks(tick_abscissa)
    ax.set_xticklabels(["stretches - " + str(cur_name) for cur_name in name_strecth_list])
    ax.set_title(u"Proportion of exons having more than 1\nstretches of " +
                 str(iunit_name) + " " + str(iunit) + " in CCE, high and low fasta files\n" +
                 "hig and low fasta file generated : " + str(iteration) + "\n"
                 + str(name_file[0]) + "\n" +  str(name_file[1]))
    up_patch = mpatches.Patch(color=colors[0], label=name_file[0])
    down_patch = mpatches.Patch(color=colors[1], label=name_file[1])
    for i in range(len(val)):
        prop = str(round(val[i] * size_fasta,1)) + "/" + str(size_fasta)
        ax.text(abscissa[i], val[i] + val[i] * 0.01 + std[i],
                prop,
                fontsize=10, horizontalalignment='center')

    for i in range(0, len(abscissa), 2):
        max_bar = max(val[i] + val[i] * 0.01 + std[i], val[i+1] + val[i+1] * 0.01 + std[i+1])
        max_bar = max_bar + (ymax * 0.035)
        ax.plot([abscissa[i], abscissa[i + 1]], [max_bar, max_bar], 'k-')
        try:
            pval = r_ttest(new_strecth_count_list[i], new_strecth_count_list[i+1])
            pval = round(pval, 3)
        except ValueError:
            pval = "nan"
        ax.text(abscissa[i], max_bar + max_bar * 0.03, "p=" + str(pval))


    plt.legend(handles=[up_patch, down_patch])
    ax.set_ylabel(u"Proportion of exons")
    ax.set_xlabel(u"Number of stretches")
    fig.add_subplot(ax)
    plt.savefig(outfile, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()



def stretch_calculator(my_file, fasta, unit_type, cur_unit, output, stretches):
    """
    :param file_up: (string) the excel file containing the up exons sequences
    :param file_down: (string) the excel file containing the down exons sequences
    :param fasta: (string) a fasta file containing random sequences
    :param unit_type: (string) the name of the unit of interest (aa, feature, nt)
    :param output: (string) the path where the graphics will be created
    """
    size = []
    st_count = []
    name_stretch = []
    #stretches=[[4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12]]
    for stretch_content, stretch_len in stretches:
        name_stretch.append(str(stretch_content) + "/" + str(stretch_len))

        nt_seq, aa_seq = read_fasta_files(my_file, stretch_len)
        if unit_type == "feature":
            seq = aa_seq

        elif unit_type == "aa":
            seq = aa_seq
        else:
            seq = nt_seq
        dic = get_stretch(seq, cur_unit, unit_type, stretch_len, stretch_content)
        st_count.append(get_more(dic,1)[0])
        size.append(dic["all"])
    if size[0] != sum(size) / len(size):
        print("Warining all size are not the same !")
    return st_count, size



#!/usr/bin/python3

# set environment
import pymysql
import os
import config

codon2aminoAcid = dict(TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L", ATT="I", ATC="I",
                       ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V", TCT="S", TCC="S", TCA="S", TCG="S",
                       CCT="P", CCC="P", CCA="P", CCG="P", ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A",
                       GCA="A", GCG="A", TAT="Y", TAC="Y", TAA="", TAG="", CAT="H", CAC="H", CAA="Q", CAG="Q",
                       AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E", TGT="C", TGC="C",
                       TGA="", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R", AGT="S", AGC="S", AGA="R", AGG="R",
                       GGT="G", GGC="G", GGA="G", GGG="G")

feature_dic = {
    "Small": ["A", "C", "D", "G", "N", "P", "S", "T", "V"], "Tiny": ["A", "C", "G", "S", "T"],
    "Aliphatic": ["A", "G", "I", "L", "V"], "Aliphatic_s": ["I", "L", "V"],
    "side_chain_aliphatic_polar": ["C", "M", "S", "T"], "Aromatic": ["F", "W", "Y", "H"], "Aromatic_s": ["F", "W", "Y"],
    "Aromatic_NP": ["F", "W"], "Sulfuric": ["C", "M"], "Hydroxylic": ["S", "T", "Y"], "Amidic": ["N", "Q"],
    "Acidic_side_chain": ["D", "N", "E", "Q"], "Basic_amino_acid": ["H", "K", "R"],
    "Hydrophobic": ["A", "C", "I", "L", "M", "F", "P", "W", "Y", "V"],
    "Hydrophobic_NP": ["A", "G", "I", "L", "M", "F", "P", "W", "V"],
    "Hydrophobic_side_chain": ["A", "I", "L", "M", "F", "W", "Y", "V"],
    "Hydrophobic_Alkyl": ["A", "G", "I", "L", "M", "P", "V"],
    "Hydrophobic_aromatic": ["F", "W"],
    "Hydrophilic": ["E", "D", "H", "K", "N", "Q", "R", "S", "T"],
    "Hydrophilic_polar": ["N", "C", "Q", "S", "T", "Y", "E", "D", "R", "H", "K"],
    "Hydrophylic_side_chain_polar": ["N", "Q", "S", "T", "Y", "E", "D", "R", "H", "K"],
    "Hydrophilic_neutral": ["N", "C", "Q", "S", "T", "Y"],
    "Hydrophilic_side_chain_uncharged": ["N", "Q", "S", "T", "Y"],
    "Hydrophilic_charged": ["E", "D", "R", "H", "K"],
    "Hydrophilic_Acidic_negative_charged": ["D", "E"],
    "Hydrophilic_Basic_positive_charged": ["R", "H", "K"],
    "Hydrophilic_positively_charged": ["R", "K"],
    "Neutral": ["A", "C", "F", "G", "I", "L", "M", "N", "P", "Q", "S", "T", "V", "W"],
    "Neutral_s": ["A", "C", "N", "Q", "S", "T", "Y"],
    "Charged": ["R", "H", "K", "D", "E"],
    "Positively_charged": ["R", "H", "K"],
    "Positively_charged_s": ["R", "K"],
    "Negatively_charged": ["D", "E"],
    "Non_polar_1": ["G", "A", "V", "L", "I", "M", "P", "F", "W"],
    "Non_polar_2": ["A", "I", "L", "M", "P", "V", "F", "W"],
    "Non_polar_1s": ["G", "A", "V", "L", "I", "M"],
    "Non_polar_alkyl": ["G", "A", "V", "L", "I", "M", "P"],
    "Non_polar_aromatic": ["F", "W"],
    "Polar": ["Y", "S", "T", "C", "Q", "N", "E", "D", "K", "H", "R"],
    "Polar_uncharged1": ["G", "S", "T", "C", "Y", "N", "Q"],
    "Polar_uncharged2": ["S", "T", "Q", "N", "C", "P"],
    "Polar_uncharged3": ["Y", "S", "T", "C", "Q", "N"],
    "Polar_uncharged4": ["S", "T", "N", "Q"],
    "Polar_charged": ["E", "D", "R", "H", "K"],
    "Polar_positively_charged": ["R", "H", "K"],
    "Polar_positively_charged_s": ["R", "K"],
    "Polar_negatively_charged": ["D", "E"],
    "Low_complexity": ["S", "P", "G", "R", "K", "Y"],
    "Disorder_promoting": ["A", "R", "G", "Q", "S", "E", "K", "P"],
    "Disorder_promoting_s": ["S", "P", "G", "R"],
    "Order_promoting": ["W", "Y", "F", "I", "L", "V", "C", "N"],
    "Thiolation": ["K", "Q", "E"],
    "EPRS": ["P", "E"],
    "PEVK": ["P", "E", "V", "K"],
    "Serine": ["S"],
    "Threonine": ["T"]
}


####################
#    functions     #
####################


def connection():
    """
    :return: an object that contains all the information you need to connect to fasterDB
    """
    cnx = pymysql.connect(user=config.user, password=config.password, host=config.host, db=config.db)
    return cnx


def translator(seq, offset):
    """
    :param seq: (string) a nucleotide sequence
    :param offset: (int) the number of nucleotide of
    the first codon (of the studied exon) contained in
    the previous exons
    :return: the translated sequence
    """
    res = ""
    if offset == 1:
        offset = 2
    elif offset == 2:
        offset = 1
    for i in range(offset, len(seq)-2, 3):
        res += codon2aminoAcid[seq[i:i+3]]
    return res


def dic_sequence(cnx, exon_type, stretch_len):
    """
    :param cnx: the information necessary to connect to fasterDB
    :param exon_type: the type of the exon of the control sets
    :param stretch_len: (int) the length of the low complexity sequence
    return a list containing the nt sequences and the amino acid sequence
    """

    cursor = cnx.cursor()
    if exon_type == "ALL":
        query = ("""SELECT cds_sequence, offset_before_exon FROM
        hsapiens_exonpeptides_filtered """)
    else:
        query = ("""SELECT t1.cds_sequence, t1.offset_before_exon FROM
                    hsapiens_exonpeptides_filtered t1, hsapiens_exonsstatus_improved t2 where t1.gene_id = t2.id_gene
                    and t1.exon_position_on_gene = t2.pos_sur_gene and t2.exon_types LIKE "%""" + str(
            exon_type) + """%" """)
    cursor.execute(query)
    query_result = cursor.fetchall()  # saving all the table "hsapiens_exonpeptides_filtered" from fasterDB
    res = []
    nt = []
    for mytuple in query_result:
        res.append(translator(mytuple[0], mytuple[1]))
        nt.append(mytuple[0])
    # we only keep the amino acids sequences greater than 5 amino acids
    final_results = []
    for sequence in res:
        if len(sequence) > stretch_len-1:
            final_results.append(sequence)

    final_nt = []
    for sequence in nt:
        if len(sequence) > stretch_len-1:
            final_nt.append(sequence)
    return final_nt, final_results


def stretch_finder_feature(sequence, feature, stretch_len, stretch_content):
    """
    :param sequence: (string) an amino acid sequences
    :param feature: (string) the name of the feature to return
    :param stretch_len: (int) the length of the low complexity sequence of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there ise a low complexity in the sub-sequence
    :return: the number of low complexity of the feature "feature" here
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
    say that there is a low complexity sequence in the sub-sequence
    :return: the number of low complexity of the aa "aa" here
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
    say that there ise a low complexity in the sub-sequence
    :return: the number of low complexity of the nucleotide "nt" here
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
    :param stretch_len: (int) the length of the low complexity of interest
    :param unit_type: (string) the name of the unit of interest (aa, feature, nt)
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there ise a low complexity in the sub-sequence
    :return: (dictionary of int) the number of sequences having 0 to 10+ low complexity sequences \
     in the list of sequences
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


def dic_creator_features(stretch_len, stretch_content):
    """
    :param stretch_len: (int) the length of the low complexity sequence of interest
    :param stretch_content: (int) the number of amino acids participating to the feature
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    if not os.path.isdir(file_dir):
        os.mkdir(file_dir + "/control_dir/")
    output = file_dir + "/control_dir/"
    cnx = connection()
    for exon_type in ["CCE"]:
        print("Creating stretches dictionaries for " + str(exon_type) + " exons : ")
        sequence_nt, sequences_aa = dic_sequence(cnx, exon_type, stretch_len)
        with open(output + exon_type + "_feature_" + str(stretch_len) + "_" + str(stretch_content) + "_dic.py", "w") \
                as outfile:
            outfile.write("# stretch_len = " + str(stretch_len) + " - stretch_content : " + str(stretch_content) + "\n")
            for feature in feature_dic.keys():
                print("   * For the feature " + str(feature))
                dic_ft = get_stretch(sequences_aa, feature, "feature", stretch_len, stretch_content)
                outfile.write(feature + "_dic = " + str(dic_ft) + "\n")
    print("Finished !")


def dic_creator_amino_acid(stretch_len, stretch_content):
    """
    :param stretch_len: (int) the length of the low complexity sequence of interest
    :param stretch_content: (int) the number of unit of interest find in the subsequence of stretch_len
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    output = file_dir + "/control_dir/"
    cnx = connection()
    for exon_type in ["CCE"]:
        print("Creating stretches dictionaries for " + str(exon_type) + " exons : ")
        sequence_nt, sequences_aa = dic_sequence(cnx, exon_type, stretch_len)
        with open(output + exon_type + "_aa_" + str(stretch_len) + "_" + str(stretch_content) + "_dic.py", "w") \
                as outfile:
            outfile.write("# stretch_len = " + str(stretch_len) + " - stretch_content : " + str(stretch_content) + "\n")
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                print("   * For the amino acid " + str(aa))
                dic_ft = get_stretch(sequences_aa, aa, "aa", stretch_len, stretch_content)
                outfile.write(aa + "_dic = " + str(dic_ft) + "\n")
    print("Finished !")


def dic_creator_nt(stretch_len, stretch_content):
    """
    :param stretch_len: (int) the length of the low complexity sequence of interest
    :param stretch_content: (int) the number of unit of interest find in the subsequence of stretch_len
    """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    output = file_dir + "/control_dir/"
    cnx = connection()
    for exon_type in ["CCE"]:
        print("Creating stretches dictionaries for " + str(exon_type) + " exons : ")
        sequence_nt, sequences_aa = dic_sequence(cnx, exon_type, stretch_len)
        with open(output + exon_type + "_nt_" + str(stretch_len) + "_" + str(stretch_content) + "_dic.py", "w") \
                as outfile:
            outfile.write("# stretch_len = " + str(stretch_len) + " - stretch_content : " + str(stretch_content) + "\n")
            for nt in "ACTGYRSWMK":
                print("   * For the nucleotide " + str(nt))
                dic_ft = get_stretch(sequence_nt, nt, "nt", stretch_len, stretch_content)
                outfile.write(nt + "_dic = " + str(dic_ft) + "\n")
    print("Finished !")


def dic_creator():
    """Create all the control dictionaries we need"""
    for stretch in [[7, 5], [8, 6], [9, 7], [10, 8], [11, 9], [10, 7], [11, 8], [12, 9]]:
        dic_creator_features(stretch[0], stretch[1])
        dic_creator_amino_acid(stretch[0], stretch[1])
        dic_creator_nt(stretch[0], stretch[1])


if __name__ == "__main__":
    dic_creator()

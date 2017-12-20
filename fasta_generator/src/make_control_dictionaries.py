"""Description: create control dictionaries that will contains the frequency of each
hexanucleotides, codons, amino acids in ACE/CCE/ALL exons of fasterDB"""


# import
from dicitonary import *
import mysql.connector
import os

codon2aminoAcid = dict(TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L", ATT="I", ATC="I",
                       ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V", TCT="S", TCC="S", TCA="S", TCG="S",
                       CCT="P", CCC="P", CCA="P", CCG="P", ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A",
                       GCA="A", GCG="A", TAT="Y", TAC="Y", TAA="", TAG="", CAT="H", CAC="H", CAA="Q", CAG="Q",
                       AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E", TGT="C", TGC="C",
                       TGA="", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R", AGT="S", AGC="S", AGA="R", AGG="R",
                       GGT="G", GGC="G", GGA="G", GGG="G")

####################
#    functions     #
####################


def connection():
    """
    :return: an object that contains all the information you need to connect to fasterDB
    """
    cnx = mysql.connector.connect(***REMOVED***, ***REMOVED***, ***REMOVED***, ***REMOVED***,
                                  buffered=True)
    return cnx


def translator(seq, offset):
    """

    :param seq: (string) a nucleotide sequence
    :param offset: (int) the number of nucleotide of
    the first codon (of the studied exon) contained in
    the previous exons
    :return: the list of codon in the sequence and the translated sequence
    """
    res = ""
    if offset == 1:
        offset =2
    elif offset == 2:
        offset =1
    codon = []
    for i in range(offset, len(seq)-2, 3):
        res += codon2aminoAcid[seq[i:i+3]]
        codon.append(seq[i:i+3])
    return [codon, res]


def dic_sequence(cnx, exon_type):
    """
    :param cnx: the information necessary to connect to fasterDB
    :param exon_type: the type of the exon of the control sets
    return a dictionary containing the sequences and the offsets of every exon in fasterDB
    """

    cursor = cnx.cursor()
    if exon_type == "ALL":
        query = ("""SELECT DISTINCT genomic_sequence, cds_sequence, offset_before_exon FROM
        hsapiens_exonpeptides_filtered """)
    else:
        query = ("""SELECT DISTINCT t1.genomic_sequence, t1.cds_sequence, t1.offset_before_exon FROM
                    hsapiens_exonpeptides_filtered t1, hsapiens_exonsstatus_improved t2 where t1.gene_id = t2.id_gene
                    and t1.exon_position_on_gene = t2.pos_sur_gene and t2.exon_types LIKE "%""" + str(
            exon_type) + """%" """)
    cursor.execute(query)
    query_result = cursor.fetchall()  # saving all the table "hsapiens_exonpeptides_filtered" from fasterDB
    res = []
    for tuple in query_result:
        res.append([tuple[0]] + translator(tuple[1], tuple[2]))
    return res


def calcul_dic(dic, seq):
    """
    :param dic: (a dictionary of float) freq of the word af interest
    :param seq: (string) the nt sequence
    :return: 'dic' completed
    """
    for i in range(len(seq)-(6-1)):
        if seq[i:i+6] not in dic.keys():
            dic[seq[i:i+6]] = 1
        else:
            dic[seq[i:i+6]] += 1
    return dic


def create_an_hexanucleotid_dic(tuple_list):
    """

    :param tuple_list: (list of list) each sublist corresponds to an exon.
    Each sublist contains the geneomic sequence, the codon sequence and the
    encoded amino acid sequence of this exons
    :return:
    """
    dic = {}
    for i in range(len(tuple_list)):
        dic = calcul_dic(dic, tuple_list[i][0])
    count = 0
    for key in dic.keys():
        count += dic[key]
    dic["all"] =  count
    return dic


def calcul_dic_codon(dic, seq):
    """
    :param dic: (dictionary of int) the number of codon for a
    given set of sequence
    :param seq: (list) list of codon
    :return: dic with the updated codon content
    """
    for codon in seq:
        if codon not in dic.keys():
            dic[codon] = 1
        else:
            dic[codon] += 1
    return  dic


def create_a_codon_dic(tuple_list):
    """
    :param tuple_list:(list of list) each sublist corresponds to an exon.
    Each sublist contains the geneomic sequence, the codon sequence and the
    encoded amino acid sequence of this exons
    :return: a dictionary that contains the number of codon
    in all exons in tuple list
    """
    dic = {}
    for i in range(len(tuple_list)):
        dic = calcul_dic_codon(dic, tuple_list[i][1])
    count = 0
    for key in dic.keys():
        count += dic[key]
    dic["all"] = count
    return dic


def calcul_dic_aa(dic, seq):
    """
    :param dic: (dictionary of int) the number of amino acid for a
    given set of sequence
    :param seq: (string) amino acid sequence
    :return: dic with the updated codon content
    """
    for aa in seq:
        if aa not in dic.keys():
            dic[aa] = 1
        else:
            dic[aa] += 1
    return  dic


def create_an_aa_dic(tuple_list):
    """
    :param tuple_list:(list of list) each sublist corresponds to an exon.
    Each sublist contains the geneomic sequence, the codon sequence and the
    encoded amino acid sequence of this exons
    :return: a dictionary that contains the number of every amino acids
    in all exons in tuple list
    """
    dic = {}
    for i in range(len(tuple_list)):
        if len(tuple_list[i][2]) > 0:
            dic = calcul_dic_aa(dic, tuple_list[i][2])
    count = 0
    for key in dic.keys():
        count += dic[key]
    dic["all"] = count
    return dic


def create_a_custom_dic(aa_dic, feature_dic):
    """

    :param aa_dic:  a dictionary that contains the number of every amino acids
    in all exons in tuple list (all ACE/CCE/ALL exons of fasterDB)
    :param feature_dic: (dictionary of list of character) link each feature
    to their corresponding amino acids
    :return: a dictionary that link for each amino acid feature their number in fasterDB
    """
    res = {}
    for key in feature_dic.keys():
        res[key] = 0
        for aa in feature_dic[key]:
            if aa in  aa_dic.keys():
                res[key] += aa_dic[aa]
    return res


def get_exons_value(list_seq, dic):
    """
    Calculate the propensity scale given in dic for the sequences in list_seq.

    Give for each sequence in list_seq, its value according to each amino acid
    values in dic.
    :param list_seq: (list of string), list of peptide sequences
    :param dic: (dictionary) each amino acid (key) is associated with a float
    value (value)
    :return: (list of float) the list of propensity value for each sequences
    in list_seq
    """
    list_val = []
    correction = True
    for key in dic.keys():
        if dic[key] < 0:
            correction = False
    for i in range(len(list_seq)):
        val = 0.
        if len(list_seq[i][2]) > 0:
            for j in range(len(list_seq[i][2])):
                val += dic[list_seq[i][2][j]]
            if correction:
                list_val.append(val / len(list_seq[i][2]))
            else:
                list_val.append(val)
    return list_val


def create_dic():
    """
    Function that create all the control dictionaries
    """
    print("start...")
    file_dir = os.path.dirname(os.path.realpath(__file__))
    ctrl = ["ACE", "CCE", "ALL"]
    cnx = connection()
    os.mkdir(file_dir + "/control_dic")
    for exon_type in ctrl:
        print("retrieving_all" + str(exon_type) + " exons")
        tuple_list = dic_sequence(cnx, exon_type)
        print("hexa")
        d6 = create_an_hexanucleotid_dic(tuple_list)
        print("codon")
        dc = create_a_codon_dic(tuple_list)
        print("aa")
        da = create_an_aa_dic(tuple_list)
        print("sh")
        sh = create_a_custom_dic(da, schain2aa)
        print("hy")
        hy = create_a_custom_dic(da, hydro_info2aa)
        print("ch")
        ch = create_a_custom_dic(da, charge_info2aa)
        print("po")
        po = create_a_custom_dic(da, polarity_info2aa)
        print("mi")
        mi = create_a_custom_dic(da, misc2aa)
        list_dic = [aa2kyte_hydrophobicity, aa2eisenberg_hydrophobicity,
                    aa2fauchere_hydrophobicity, aa2zimmerman_polarity,
                    aa2grantham_polarity, aa2deleage_alpha, aa2levitt_alpha,
                    aa2chou_alpha, aa2nagano_beta, aa2deleage_beta, aa2chou_beta,
                    aa2deleage_bturn, aa2levitt_bturn, aa2chou_bturn, aa2nagano_coil,
                    aa2deleage_coil]
        scale = ["hydrophobicity_kyte", "hydrophobicity_eisenberg",
                 "hydrophobicity_fauchere", "polarity_zimmerman",
                 "polarity_grantham", "alpha_helix_prediction_deleage",
                 "alpha_helix_prediction_levitt", "alpha_helix_prediction_chou",
                 "beta_helix_prediction_nagano", "beta_helix_prediction_deleage",
                 "beta_helix_prediction_chou", "beta_turn_prediction_deleage",
                 "beta_turn_prediction_levitt", "beta_turn_prediction_chou",
                 "coil_prediction_nagano", "coil_prediction_deleage"]
        propensities_dic = []
        for i in range(len(list_dic)):
            print(scale[i])
            propensities_dic.append(get_exons_value(tuple_list, list_dic[i]))
        print("writing")
        with open(file_dir + "/control_dic/" + exon_type + "_dic.py", "w") as out_file:
            out_file.write("d6 = " + str(d6) + "\n")
            out_file.write("dc = " + str(dc) + "\n")
            out_file.write("da = " + str(da) + "\n")
            out_file.write("sh = " + str(sh) + "\n")
            out_file.write("hy = " + str(hy) + "\n")
            out_file.write("ch = " + str(ch) + "\n")
            out_file.write("po = " + str(po) + "\n")
            out_file.write("mi = " + str(mi) + "\n")
            for i in range(len(list_dic)):
                out_file.write(str(scale[i]) + " = " + str(list_dic[i]) + "\n")

if __name__ == "__main__":
    create_dic()

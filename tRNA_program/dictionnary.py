# File description :
# this file contains all the dictionaries necessary to retrieve the anticodon, the amino-acid and its nature
# associated with a codon.

# CODON-AMINO_ACID DICTIONARY
# Links each codon to its amino_acid
codon2aminoAcid = dict(TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L", ATT="I", ATC="I",
                       ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V", TCT="S", TCC="S", TCA="S", TCG="S",
                       CCT="P", CCC="P", CCA="P", CCG="P", ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A",
                       GCA="A", GCG="A", TAT="Y", TAC="Y", TAA="*", TAG="*", CAT="H", CAC="H", CAA="Q", CAG="Q",
                       AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E", TGT="C", TGC="C",
                       TGA="*", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R", AGT="S", AGC="S", AGA="R", AGG="R",
                       GGT="G", GGC="G", GGA="G", GGG="G")

# codon-anticodon dictionary
# Links each codon to its respective anticodon
codon2anticodon = dict(TTT="GAA", TTC="GAA", TTA="TAA", TTG="TAA,CAA", TCT="IGA,GGA", TCC="IGA,GGA", TCA="IGA,TGA",
                       TCG="CGA,TGA", TAT="GTA", TAC="GTA", TAA="TTA", TAG="CTA,TTA", TGT="GCA", TGC="GCA", TGA="TCA",
                       TGG="CCA", CTT="IAG,GAG", CTC="IAG,GAG", CTA="IAG,TAG", CTG="CAG,TAG", CCT="IGG,GGG",
                       CCC="IGG,GGG", CCA="IGG,TGG", CCG="CGG,TGG", CAT="GTG", CAC="GTG", CAA="TTG", CAG="CTG,TTG",
                       CGT="ICG,GCG", CGC="ICG,GCG", CGA="ICG,TCG", CGG="CCG,TCG", ATT="IAT,GAT", ATC="IAT,GAT",
                       ATA="IAT,TAT", ATG="CAT", ACT="IGT,GGT", ACC="IGT,GGT", ACA="IGT,TGT", ACG="CGT,TGT", AAT="GTT",
                       AAC="GTT", AAA="TTT", AAG="CTT,TTT", AGT="GCT", AGC="GCT", AGA="TCT", AGG="CCT,TCT",
                       GTT="IAC,GAC", GTC="IAC,GAC", GTA="IAC,TAC", GTG="CAC,TAC", GCT="IGC,GGC", GCC="IGC,GGC",
                       GCA="IGC,TGC", GCG="CGC,TGC", GAT="GTC", GAC="GTC", GAA="TTC", GAG="TTC,CTC", GGT="ICC,GCC",
                       GGC="ICC,GCC", GGA="ICC,TCC", GGG="CCC,TCC")

# codon2rare dictionary
# link each codon to its rareness
# source atgme.org based on the codon usage database for human
codon2rare = dict(TTT="", TTC="", TTA="-", TTG="", CTT="", CTC="", CTA="-", CTG="", ATT="", ATC="",
                  ATA="-", ATG="", GTT="", GTC="", GTA="-", GTG="", TCT="", TCC="", TCA="", TCG="-",
                  CCT="", CCC="", CCA="", CCG="-", ACT="", ACC="", ACA="", ACG="-", GCT="", GCC="",
                  GCA="", GCG="-", TAT="- ", TAC="", TAA="", TAG="", CAT="", CAC="", CAA="", CAG="",
                  AAT="", AAC="", AAA="", AAG="", GAT="", GAC="", GAA="", GAG="", TGT="", TGC="",
                  TGA="", TGG="", CGT="-", CGC="", CGA="-", CGG="", AGT="", AGC="", AGA="", AGG="",
                  GGT="", GGC="", GGA="", GGG="")

# amino_acid2codon
# Links each amino acid to those corresponding codon
amino_acid2codon = {
    "F": "TTT,TTC", "L": "TTA,TTG,CTT,CTC,CTA,CTG", "I": "ATT,ATC,ATA", "M": "ATG", "V": "GTT,GTC,GTA,GTG",
    "S": "TCT,TCC,TCA,TCG,AGT,AGC", "P": "CCT,CCC,CCA,CCG", "T": "ACT,ACC,ACA,ACG", "A": "GCT,GCC,GCA,GCG",
    "Y": "TAT,TAC", "*": "TAA,TAG,TGA", "H": "CAT,CAC", "Q": "CAA,CAG", "N": "AAT,AAC", "K": "AAA,AAG", "D": "GAT,GAC",
    "E": "GAA,GAG", "C": "TGT,TGC", "W": "TGG", "R": "CGT,CGC,CGA,CGG,AGA,AGG", "G": "GGT,GGC,GGA,GGG"}

# Links each amino acid to those corresponding codon
amino_acid2codon_list = {
    "F": ["TTT", "TTC"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "I": ["ATT", "ATC", "ATA"], "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"], "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"], "T": ["ACT", "ACC", "ACA", "ACG"], "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"], "*": ["TAA", "TAG", "TGA"], "H": ["CAT", "CAC"], "Q": ["CAA", "CAG"], "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"], "D": ["GAT", "GAC"], "E": ["GAA", "GAG"], "C": ["TGT", "TGC"], "W": ["TGG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], "G": ["GGT", "GGC", "GGA", "GGG"]}

# nature2amino_acid
# link each defined nature to there corresponding amino acid

nature2amino_acid = {
    "NP": ["G", "A", "V", "L", "I", "M", "P", "F", "W"], "NP-Alkyl": ["G", "A", "V", "L", "I", "M", "P"],
    "NP-aromatic":["F", "W"], "NPR": ["G", "A", "V", "L", "I", "M"],
    "P": ["Y", "S", "T", "C", "Q", "N", "E", "D", "K", "H", "R"], "PN": ["Y", "S", "T", "C", "Q", "N"],
    "PNC": ["S", "T", "N", "Q"], "PNC1": ["S", "T", "Q", "N" "C", "P"], "PNC2": ["S", "T", "Q", "N"],
    "PC": ["E", "D", "R", "H", "K"], "P-NC": ["E", "D"], "P-PC": ["K", "H", "R"],
    "HC": ["A", "V", "I", "L", "M", "F", "Y", "W"], "H": ["A", "V", "L", "I", "M", "F", "W", "P", "G"],
    "Aliphatic": ["G", "A", "V", "L", "I"], "HS": ["S", "C", "T", "M"], "Aromatic": ["F", "Y", "W"]
}

# Metabolism2amino_acid


metabolism2amino_acid = {
    "Non_essential": ["A", "D", "N", "E", "S"], "Essential": ["I", "L", "M", "V", "F", "W", "H", "T", "K"],
    "Conditionally_essential": ["R", "C", "Q", "Y", "G", "P"], "Glycolyse": ["L", "V", "A", "S", "G", "C"],
    "TCA_cycle": ["I", "M", "T", "K", "D", "N", "E", "R", "Q", "P"], "Pentose": ["F", "W", "H", "Y"],
    "Thiolation": ["K", "Q", "E"]
}


# amino_acid-nature dictionary
# Links each amino_acid to its specific nature
amino_acid2nature = {'G': "Nonpolar", 'A': "Nonpolar", 'V': "Nonpolar", 'L': "Nonpolar", 'M': "Nonpolar",
                     'I': "Nonpolar", 'S': "Polar", 'T': "Polar", 'C': "Polar", 'P': "Polar", 'N': "Polar",
                     'Q': "Polar", 'F': "Aromatic", 'Y': "Aromatic", 'W': "Aromatic", 'K': "Positively charged",
                     'R': "Positively charged", 'H': "Positively charged", 'D': 'Negatively charged',
                     'E': "Negatively charged", "*": "None"}

# codon2rareness
# source atg.me for rare codon (--) and fasterDB+the homo sapiens usage database from the codon usage database
# for each codon says if the codon is rare (--) or has a the weakest frequency (-) in comparison of the other codons
# that codes for the same amino acid
codon2rareness = dict(TTT="-", TTC="", TTA="--", TTG="", CTT="", CTC="", CTA="--", CTG="", ATT="", ATC="",
                      ATA="--", ATG="", GTT="", GTC="", GTA="--", GTG="", TCT="", TCC="", TCA="", TCG="--",
                      CCT="", CCC="", CCA="", CCG="--", ACT="", ACC="", ACA="", ACG="--", GCT="", GCC="",
                      GCA="", GCG="--", TAT="", TAC="", TAA="", TAG="", CAT="-", CAC="", CAA="-", CAG="",
                      AAT="-", AAC="", AAA="-", AAG="", GAT="-", GAC="", GAA="-", GAG="", TGT="-", TGC="",
                      TGA="", TGG="", CGT="--", CGC="", CGA="--", CGG="", AGT="", AGC="", AGA="", AGG="",
                      GGT="-", GGC="", GGA="", GGG="")



def dic_sequence(cnx, exon_type):
    """
    :param cnx: the information necessary to connect to fasterDB
    :param exon_type: the type of the exon of the control sets
    return a dictionary containing the sequences and the offsets of every exon in fasterDB
    """

    cursor = cnx.cursor()
    if exon_type == "ALL":
        query = ("""SELECT DISTINCT cds_sequence, offset_before_exon, gene_id, exon_position_on_gene FROM
        hsapiens_exonpeptides_filtered """)
    else:
        query = ("""SELECT DISTINCT t1.cds_sequence, t1.offset_before_exon, t1.gene_id, t1.exon_position_on_gene FROM
                    hsapiens_exonpeptides_filtered t1, hsapiens_exonsstatus_improved t2 where t1.gene_id = t2.id_gene
                    and t1.exon_position_on_gene = t2.pos_sur_gene and t2.exon_types LIKE "%""" + str(
            exon_type) + """%" """)
    cursor.execute(query)
    query_result = cursor.fetchall()  # saving all the table "hsapiens_exonpeptides_filtered" from fasterDB
    d = {}
    for i in range(len(query_result)):
        # filling a dictionary with a key corresponding to "gene_name_exon_position_on_gene" (that identifies
        # an unique exon) linked with its respective sequence
        d[str(query_result[i][2]) + "_" + str(query_result[i][3])] = query_result[i][0], query_result[i][1]
    return d

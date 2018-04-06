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



# feature dictionaries
# 12/03/2018
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


# nature2amino_acid
# link each defined nature to there corresponding amino acid

nature2amino_acid = {
    "NP": ["G", "A", "V", "L", "I", "M", "P", "F", "W"], "NP-Alkyl": ["G", "A", "V", "L", "I", "M", "P"],
    "NP-aromatic": ["F", "W"], "NPR": ["G", "A", "V", "L", "I", "M"],
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

# side chain properties:
schain2aa = {
    "Small": ["A", "C", "D", "G", "N", "P", "S", "T", "V"], "Tiny": ["A", "C", "G", "S", "T"],
    "Aliphatic": ["A", "G", "I", "L", "V"], "Aliphatic_s": ["I", "L", "V"],
    "side_chain_aliphatic_polar": ["C", "M", "S", "T"], "Aromatic": ["F", "W", "Y", "H"], "Aromatic_s": ["F", "W", "Y"],
    "Aromatic_NP": ["F", "W"], "Sulfuric": ["C", "M"], "Hydroxylic": ["S", "T", "Y"], "Amidic": ["N", "Q"],
    "Acidic_side_chain": ["D", "N", "E", "Q"], "Basic_amino_acid": ["H", "K", "R"]
}

# hydro properties
hydro_info2aa = {
    "Hydrophobic": ["A", "C", "I", "L", "M", "F", "P", "W", "Y", "V"],
    "Hydrophobic_NP": ["A", "G", "I", "L", "M", "F", "P", "W", "V"],
    "Hydrophobic_side_chain": ["A", "I", "L", "M", "F", "W", "Y", "V"],
    "Hydrophobic-Alkyl": ["A", "G", "I", "L", "M", "P", "V"],
    "Hydrophobic-aromatic": ["F", "W"],
    "Hydrophilic": ["E", "D", "H", "K", "N", "Q", "R", "S", "T"],
    "Hydrophilic_polar": ["N", "C", "Q", "S", "T", "Y", "E", "D", "R", "H", "K"],
    "Hydrophylic_side_chain_polar": ["N", "Q", "S", "T", "Y", "E", "D", "R", "H", "K"],
    "Hydrophilic_neutral": ["N", "C", "Q", "S", "T", "Y"],
    "Hydrophilic_side_chain_uncharged": ["N", "Q", "S", "T", "Y"],
    "Hydrophilic_charged": ["E", "D", "R", "H", "K"],
    "Hydrophilic_Acidic_negative_charged": ["D", "E"],
    "Hydrophilic_Basic_positive_charged": ["R", "H", "K"],
    "Hydrophilic_positively_charged": ["R", "K"]
}

# charge properties
charge_info2aa = {
    "Neutral": ["A", "C", "F", "G", "I", "L", "M", "N", "P", "Q", "S", "T", "V", "W"],
    "Neutral_s": ["A", "C", "N", "Q", "S", "T", "Y"],
    "Charged": ["R", "H", "K", "D", "E"],
    "Positively_charged": ["R", "H", "K"],
    "Positively_charged_s": ["R", "K"],
    "Negatively_charged": ["D", "E"]
}

# polarity properties of amino acids
polarity_info2aa = {
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
    "Polar_negatively_charged": ["D", "E"]
}

# some miscellaneous properties of amino acids
misc2aa = {
    "Low_complexity": ["S", "P", "G", "R", "K", "Y"],
    "Disorder_promoting": ["A", "R", "G", "Q", "S", "E", "K", "P"],
    "Disorder_promoting_s": ["S", "P", "G", "R"],
    "Order_promoting": ["W", "Y", "F", "I", "L", "V", "C", "N"],
    "Thiolation": ["K", "Q", "E"],
    "EPRS": ["P", "E"],
    "PEVK": ["P", "E", "V", "K"]
}

# scale calculation

# ####################### FIRST GROUP ######################################""
"""
aa2dawson_size = {
    "G": 0.5, "A":2.5, "D":2.5, "C":3.0, "S":3.0, "N":5.0, "E":5.0, "T":5.0, "V":5.0, "I":5.5, "L":5.5, "P":5.5,
    "Q":6.0, "H":6.0, "M":6.0, "F":6.5, "K":7.0, "W":7.0, "Y":7.0, "R": 7.5
}
"""

# Molecular weight of each amino acid - Source : ProtScale - Reference : Most textbook
aa2molecular_weight = {
    "A": 89., "R": 174., "N": 132., "D": 133., "C": 121., "Q": 146., "E": 147., "G": 75., "H": 155., "I": 131.,
    "L": 131., "K": 146., "M": 149., "F": 165., "P": 115., "S": 105., "T": 119., "W": 204., "Y": 181., "V": 117.
}


# The aliphatic index of a protein is defined as the relative volume occupied by aliphatic side chains
# (alanine, valine, isoleucine, and leucine). It may be regarded as a positive factor for the increase of
# thermostability of globular proteins
# Source : ProtParam
# Reference : Ikai, A.J. (1980) Thermostability and aliphatic index of globular proteins. J. Biochem. 88, 1895-1898.
# [PubMed: 7462208]
def aliphatic_index(seq):
    a = 2.9
    b = 3.9
    ala = float(seq.count("A")) / len(seq)
    val = float(seq.count("V")) / len(seq)
    ile = float(seq.count("I")) / len(seq)
    leu = float(seq.count("L")) / len(seq)
    return ala * 100 + a * val * 100 + b * (ile + leu) * 100

# hydrophobicity of amino acids - Source Composition profiler
# Reference :
# Eisenberg D, Schwarz E, Komaromy M, and Wall R. (1984) "Analysis of membrane and surface protein sequences with
# the hydrophobic moment plot." J. Mol. Biol. 179:125-142.
aa2eisenberg_hydrophobicity = {
    "R": -2.53, "K": -1.5, "D": -0.9, "Q": -0.85, "N": -0.78, "E": -0.74, "H": -0.4, "S": -0.18, "T": -0.05,
    "P": 0.12, "Y": 0.26, "C": 0.29, "G": 0.48, "A": 0.62, "M": 0.64, "W": 0.81, "L": 1.06, "V": 1.08, "F": 1.19,
    "I": 1.38
}

# hydrophobicity of amino acids - Source Composition profiler
# Reference : Kyte J, and Doolittle RF. (1982) "A simple method for displaying the hydropathic character of a protein."
#  J. Mol. Biol. 157:105-132
aa2kyte_hydrophobicity = {
    "R": -4.5, "K": -3.9, "D": -3.5, "E": -3.5, "N": -3.5, "Q": -3.5, "H": -3.2, "P": -1.6, "Y": -1.3, "W": -0.9,
    "S": -0.8, "T": -0.7, "G": -0.4, "A": 1.8, "M": 1.9, "C": 2.5, "F": 2.8, "L": 3.8, "V": 4.2, "I": 4.5
}


# hydrophobicity of amino acids - Source Composition profiler
# Reference : Fauchere J-L, and Pliska VE. (1983) "Hydrophobic parameters pi of amino acid side chains from
# partitioning of N-acetyl-amino-acid amides." Eur. J. Med. Chem. 18:369-375.
aa2fauchere_hydrophobicity = {
    "R": -1.01, "K": -0.99, "D": -0.77, "E": -0.64, "N": -0.6, "Q": -0.22, "H": 0.13, "P": 0.72, "Y": 0.96, "W": 2.25,
    "S": -0.04, "T": 0.26, "G": 0., "A": 0.31, "M": 1.23, "C": 1.54, "F": 1.79, "L": 1.7, "V": 1.22, "I": 1.8
}


# polarity of amino acids - Source Composition profiler
# Reference : Zimmerman JM, Eliezer N, and Simha R. (1968) "The characterization of amino acid sequences in
# proteins by statistical methods." J J. Theor. Biol. 21:170-201.
aa2zimmerman_polarity = {
    "A": 0, "G": 0, "I": 0.13, "L": 0.13, "V": 0.13, "F": 0.35, "M": 1.43, "C": 1.48, "P": 1.58, "Y": 1.61, "T": 1.66,
    "S": 1.67, "W": 2.1, "N": 3.38, "Q": 3.53, "K": 49.5, "D": 49.7, "E": 49.9, "H": 51.6, "R": 52
}

# polarity of amino acids - Source Composition profiler
# Reference : Grantham R. Science 185:862-864(1974).
aa2grantham_polarity = {
    "A": 8.100, "R": 10.500, "N": 11.600, "D": 13.000, "C": 5.500, "Q": 10.500, "E": 12.300,
    "G": 9.000, "H": 10.400, "I": 5.200, "L": 4.900, "K": 11.300, "M": 5.700, "F": 5.200,
    "P": 8.000, "S": 9.200, "T": 8.600, "W": 5.400, "Y": 6.200, "V": 5.900
}

# ################################ SECOND GROUP ##########################################################"


# Alpha helix frequency (Nagano, 1973)
# Nagano K. (1973) "Local analysis of the mechanism of protein folding. I. Prediction of helices, loops,
# and beta-structures from primary structure." J. Mol. Biol. 75:401-420
aa2nagano_alpha = {
    "Y": 0.63, "P": 0.70, "G": 0.72, "N": 0.77, "S": 0.78, "R": 0.83, "T": 0.87,
    "C": 0.94, "I": 0.94, "V": 0.97, "D": 1.00, "W": 1.06, "Q": 1.10, "L": 1.23,
    "K": 1.23, "M": 1.23, "F": 1.23, "A": 1.29, "H": 1.29, "E": 1.54
}


# Conformational parameter for alpha helix. - source Protscale
# Reference - Deleage G., Roux B. Protein Engineering 1:289-294(1987).
aa2deleage_alpha = {
    "A": 1.489, "R": 1.224, "N": 0.772, "D": 0.924, "C": 0.966, "Q": 1.164, "E": 1.504,
    "G": 0.510, "H": 1.003, "I": 1.003, "L": 1.236, "K": 1.172, "M": 1.363, "F": 1.195,
    "P": 0.492, "S": 0.739, "T": 0.785, "W": 1.090, "Y": 0.787, "V": 0.990
}

# Normalized frequency for alpha helix.  - source protscale
# Reference :  Levitt M. -  Biochemistry 17:4277-4285(1978).
aa2levitt_alpha = {
    "A": 1.290, "R": 0.960, "N": 0.900, "D": 1.040, "C": 1.110, "Q": 1.270, "E": 1.440,
    "G": 0.560, "H": 1.220, "I": 0.970, "L": 1.300, "K": 1.230, "M": 1.470, "F": 1.070,
    "P": 0.520, "S": 0.820, "T": 0.820, "W": 0.990, "Y": 0.720, "V": 0.910
}

# Conformational parameter for alpha helix (computed from 29 proteins).  - source protscale
# Reference Chou P.Y., Fasman G.D. - Adv. Enzym. 47:45-148(1978).
aa2chou_alpha = {
    "A": 1.420, "R": 0.980, "N": 0.670, "D": 1.010, "C": 0.700, "Q": 1.110, "E": 1.510,
    "G": 0.570, "H": 1.000, "I": 1.080, "L": 1.210, "K": 1.160, "M": 1.450, "F": 1.130,
    "P": 0.570, "S": 0.770, "T": 0.830, "W": 1.080, "Y": 0.690, "V": 1.060
}

# Beta structure frequency (Nagano, 1973) - source composition profiler
# Reference : Nagano K. (1973) "Local analysis of the mechanism of protein folding. I. Prediction of helices, loops,
# and beta-structures from primary structure." J. Mol. Biol. 75:401-420.
aa2nagano_beta = {
    "E": 0.33, "R": 0.67, "N": 0.72, "P": 0.75, "S": 0.77, "K": 0.81, "H": 0.87,
    "D": 0.9, "G": 0.9, "A": 0.96, "Y": 1.07, "C": 1.13, "W": 1.13, "Q": 1.18,
    "T": 1.23, "L": 1.26, "M": 1.29, "F": 1.37, "V": 1.41, "I": 1.54
}


# Conformational parameter for beta-sheet. - source protscale
# Deleage G., Roux B. -  Protein Engineering 1:289-294(1987).
aa2deleage_beta = {
    "A": 0.709, "R": 0.920, "N": 0.604, "D": 0.541, "C": 1.191, "Q": 0.840, "E": 0.567,
    "G": 0.657, "H": 0.863, "I": 1.799, "L": 1.261, "K": 0.721, "M": 1.210, "F": 1.393,
    "P": 0.354, "S": 0.928, "T": 1.221, "W": 1.306, "Y": 1.266, "V": 1.965
}


# Conformational parameter for beta-sheet (computed from 29 proteins).  - source protparam
# Chou P.Y., Fasman G.D. - Adv. Enzym. 47:45-148(1978).
aa2chou_beta = {
    "A": 0.830, "R": 0.930, "N": 0.890, "D": 0.540, "C": 1.190, "Q": 1.100, "E": 0.370,
    "G": 0.750, "H": 0.870, "I": 1.600, "L": 1.300, "K": 0.740, "M": 1.050, "F": 1.380,
    "P": 0.550, "S": 0.750, "T": 1.190, "W": 1.370, "Y": 1.470, "V": 1.700
}


# Conformational parameter for beta-turn. source - protscale
# Deleage G., Roux B. Protein Engineering 1:289-294(1987).
aa2deleage_bturn = {
    "A": 0.788, "R": 0.912, "N": 1.572, "D": 1.197, "C": 0.965, "Q": 0.997, "E": 1.149,
    "G": 1.860, "H": 0.970, "I": 0.240, "L": 0.670, "K": 1.302, "M": 0.436, "F": 0.624,
    "P": 1.415, "S": 1.316, "T": 0.739, "W": 0.546, "Y": 0.795, "V": 0.387
}

# Normalized frequency for beta-turn.  - source protscale
# Levitt M. Biochemistry 17:4277-4285(1978).
aa2levitt_bturn = {
    "A": 0.770, "R": 0.880, "N": 1.280, "D": 1.410, "C": 0.810, "Q": 0.980, "E": 0.990,
    "G": 1.640, "H": 0.680, "I": 0.510, "L": 0.580, "K": 0.960, "M": 0.410, "F": 0.590,
    "P": 1.910, "S": 1.320, "T": 1.040, "W": 0.760, "Y": 1.050, "V": 0.470
}

# Conformational parameter for beta-turn (computed from 29 proteins).
#  Chou P.Y., Fasman G.D.  -  Adv. Enzym. 47:45-148(1978).
aa2chou_bturn = {
    "A": 0.660, "R": 0.950, "N": 1.560, "D": 1.460, "C": 1.190, "Q": 0.980, "E": 0.740,
    "G": 1.560, "H": 0.950, "I": 0.470, "L": 0.590, "K": 1.010, "M": 0.600, "F": 0.600,
    "P": 1.520, "S": 1.430, "T": 0.960, "W": 0.960, "Y": 1.140, "V": 0.500
}


# Coil propensity (Nagano, 1973)
# Nagano K. (1973) "Local analysis of the mechanism of protein folding. I. Prediction of helices, loops,
# and beta-structures from primary structure." J. Mol. Biol. 75:401-420.

aa2nagano_coil = {
    "F": 0.58, "M": 0.62, "L": 0.63, "A": 0.72, "E": 0.75, "H": 0.76, "I": 0.8,
    "Q": 0.81, "V": 0.83, "K": 0.84, "W": 0.87, "C": 1.01, "T": 1.03, "D": 1.04,
    "R": 1.33, "S": 1.34, "G": 1.35, "Y": 1.35, "N": 1.38, "P": 1.43
}


# Conformational parameter for coil.
# Deleage G., Roux B. - Protein Engineering 1:289-294(1987).
aa2deleage_coil = {
    "A": 0.824, "R": 0.893, "N": 1.167, "D": 1.197, "C": 0.953, "Q": 0.947, "E": 0.761,
    "G": 1.251, "H": 1.068, "I": 0.886, "L": 0.810, "K": 0.897, "M": 0.810, "F": 0.797,
    "P": 1.540, "S": 1.130, "T": 1.148, "W": 0.941, "Y": 1.109, "V": 0.772
}


# Transmembrane tendency - source - protscale
# Zhao, G., London E. - Protein Sci. 15:1987-2001(2006).
aa2zhao_mb = {
    "A": 0.380, "R": -2.570, "N": -1.620, "D": -3.270, "C": -0.300, "Q": -1.840, "E": -2.900,
    "G": -0.190, "H": -1.440, "I": 1.970, "L": 1.820, "K": -3.460, "M": 1.400, "F": 1.980,
    "P": -1.440, "S": -0.530, "T": -0.320, "W": 1.530, "Y": 0.490, "V": 1.460
}


# Composition of amino acids in membrane proteins (percent), Cedano et al., J. Mol. Biol. 1997, 266:594-600
aa2cedano_mb = {
    "A": 8.1, "R": 4.6, "N": 3.7, "D": 3.8, "C": 2.0, "Q": 3.1, "E": 4.6,
    "G": 7.0, "H": 2.0, "I": 6.7, "L": 11.0, "K": 4.4, "M": 2.8, "F": 5.6,
    "P": 4.7, "S": 7.3, "T": 5.6, "W": 1.8, "Y": 3.3, "V": 7.7
}

# Transmembrane regions of non-mt-proteins, Nakashima et al., Proteins 1990, 8:173-178
aa2nakashima_mb = {
    "A": 10.17, "L": 16.22, "R": 1.21, "K": 1.04, "N": 1.36, "M": 4.12, "D": 1.18, "F": 9.60, "C": 1.48,
    "P": 2.24, "Q": 1.57, "S": 5.38, "E": 1.15, "T": 5.61, "G": 8.87, "W": 2.67, "H": 1.07, "Y": 2.68,
    "I": 10.91, "V": 11.44
}

# Normalized composition of membrane proteins, Nakashima et al., Proteins 1990, 8:173-178
aa2nakasima_mb2 = {
    "A": 0.34, "L": 0.52, "R": -0.57, "K": -0.75, "N": -0.27, "M": 0.47, "D": -0.56, "F": 1.30, "C": -0.32,
    "P": -0.19, "Q": -0.34, "S": -0.20, "E": -0.43, "T": -0.04, "G": 0.48, "W": 0.77, "H": -0.19, "Y": 0.07,
    "I": 0.39, "V": 0.36
}

########################


def dic_first_group(seq, dic, penalty_size, c_value, group):
    """
    :param seq: (string) an amino acid sequence
    :param dic:  (dictionary) a dictionary with as key, the name of the group we want to compute
    :param penalty_size: (int) size below which, exons will be penalized
    :param c_value: (list of float) the value of penalty attributed for each exon studied
    :param group: (string) "1" for chimiacl group or "2"  for structural group
    :return: the dictionary updated  with the current sequence and the c_value actualised with the current sequence
    """
    if group == "1":
        list_name_dic = ["polarity(Zimmerman, 1968)", "Polarity (Grantham, 1974)"]
        list_dic = [aa2zimmerman_polarity, aa2grantham_polarity]
        mol = 0.
        hfbct_e = 0.
        hfbct_k = 0.
        hfbct_f = 0.
        for letter in seq:
            mol += aa2molecular_weight[letter]
            hfbct_e += aa2eisenberg_hydrophobicity[letter]
            hfbct_k += aa2kyte_hydrophobicity[letter]
            hfbct_f += aa2fauchere_hydrophobicity[letter]

        dic["aliphatic_index"] += round(aliphatic_index(seq), 4)
        dic["molecular_weight"] += mol
        dic["hydrophobicity(Eisenberg, 1984)"] += hfbct_e
        dic["hydrophobicity(Kyte, 1982)"] += hfbct_k
        dic["hydrophobicity(Fauchere, 1983)"] += hfbct_f

    else:
        list_name_dic = ["Alpha_helix_frequency(Nagano)", "Alpha_helix(Deleage&Roux)", "Alpha_helix(Levitt)",
                         "Alpha_helix(Chou&Fasman)", "Beta_structure_frequency(Nagano)",
                         "Beta_sheet(Deleage&Roux)", "Beta_sheet(Chou&Fasman)", "Beta_turn(Deleage&Roux)",
                         "Beta_turn(Levitt)", "Beta_turn(Chou&Fasman)", "Coil_propensity(Nagano)",
                         "Coil(Deleage&Roux)", "AA_composition_of_mb_p",
                         "transmenbrane_region_aa"]
        list_dic = [aa2nagano_alpha, aa2deleage_alpha, aa2levitt_alpha, aa2chou_alpha, aa2nagano_beta, aa2deleage_beta,
                    aa2chou_beta, aa2deleage_bturn, aa2levitt_bturn, aa2chou_bturn, aa2nagano_coil, aa2deleage_coil,
                    aa2cedano_mb, aa2nakashima_mb]
        zhao_mb = 0.
        nakasima_mb2 = 0
        for letter in seq:
            zhao_mb += aa2zhao_mb[letter]
            nakasima_mb2 += aa2nakasima_mb2[letter]

        dic["Transmenbrane_tendancy(Zhao)"] += zhao_mb
        dic["composition_of_mb_p"] += nakasima_mb2

    for i in range(len(list_name_dic)):
        count = 0.
        for letter in seq:
            count += list_dic[i][letter]
        if len(seq) > penalty_size:
            count /= len(seq)
            dic[list_name_dic[i]] += round(count, 4)
        else:
            count = (count / len(seq)) * len(seq) / penalty_size
            dic[list_name_dic[int(i)]] += round(count, 4)
    if len(seq) > penalty_size:
        c_value.append(1)
    else:
        c_value.append(float(len(seq)) / penalty_size)

    return dic, c_value


#################################################################################################################
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

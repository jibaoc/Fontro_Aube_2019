# Imports of dictionary
# Imports of dictionary
from dictionnary import codon2aminoAcid  # Links each codon to its respective anticodon



class ExonSetClass:
    """
    A class corresponding to an exon from the control set. This class allows to extract easily information about an exon
    This class is preferably use for exon of the control sets because a lot of information are not necessary for those
    exons compared to those given by the user
    """

    def __init__(self, init_tuple):
        """
        :param init_tuple: A tuple containing a gene id and the position of the interest exon on this gene.
        """
        self.gene_id = str(init_tuple[0])
        self.exon_number = str(init_tuple[1])
        self.cds_sequence = str()
        self.codon = list()
        self.offset = int()
        self.amino_acid = list()
        self.nature = list()
        self.metabolism = list()
        self.importance = list()

    def retrieve_exon_sequences(self, dic):
        """
        :param dic: a dictionary containing the sequences and the offset (before the exon of interest) of every exon in
        fasterDB. The keys of this dictionary correspond to the concatenation of a gene id and the position of the exon
        in this gene.
        fill the exon object with the cds_sequence of the exon and its offset (before the exon)
        and their offsets
        """
        sequences = dic[str(self.gene_id)+"_"+str(self.exon_number)]
        self.cds_sequence = str(sequences[0])
        self.offset = int(sequences[1])

    def found_codon_anticodon_amino_acid_and_aa_nature(self):
        """
        Fills the list of codons, amino_acids and their nature for each exon
        """
        offset_dict = {1: 2, 2: 1, 0: 0}
        for i in range(offset_dict[self.offset], len(self.cds_sequence), 3):
            if len(self.cds_sequence[i:i + 3]) == 3:
                self.codon.append(self.cds_sequence[i:i + 3])
                self.amino_acid.append(codon2aminoAcid[self.cds_sequence[i:i + 3]])

    def found_codon_amino_acid_and_aa_nature(self):
        """
        Fills the list of codons and amino_acids for each exon
        """
        offset_dict = {1: 2, 2: 1, 0: 0}
        for i in range(offset_dict[self.offset], len(self.cds_sequence), 3):
            if len(self.cds_sequence[i:i + 3]) == 3:
                self.codon.append(self.cds_sequence[i:i + 3])
                self.amino_acid.append(codon2aminoAcid[self.cds_sequence[i:i + 3]])


# Imports of dictionary
from dictionnary import codon2anticodon  # Links each codon to its amino_acid
from dictionnary import codon2aminoAcid  # Links each codon to its respective anticodon
from dictionnary import amino_acid2codon  # Links each amino_acid to its specific codons
# a non essential or a conditionally essential amino acid



class ExonClass:
    """
    A class corresponding to an exon. This class allows to extract easily information about an exon
    """

    def __init__(self, init_tuple, input_name, input_chromosome_number):
        """
        :param init_tuple: a tuple containing many information about the exon in the order specified below :
        (gene_id : the fasterDB id of the exon's gene, exon_number : the exon position within the gene,
        start_on_chromosome : the chromosomal coordinates  where the exon start,
        end_on_chromosome : the chromosomal coordinates where the exon stop,
        cds_start_on_chromosome : the chromosomal coordinates  where the cds of exon start,
        cds_end_on_chromosome : the chromosomal coordinates where the cds of exon stop,
        exon_type : the type of the exon (ACE, FCE, LCE, CCE))
        strand  : strand of the exon
        """

        self.exon_name = str(input_name)
        self.chr = str(input_chromosome_number)
        self.gene_id = str(init_tuple[0])
        self.exon_number = str(init_tuple[1])
        self.exon_start = init_tuple[2]
        self.exon_end = init_tuple[3]
        self.cds_start = init_tuple[4]
        self.cds_end = init_tuple[5]
        self.exon_type = str(init_tuple[6])
        self.strand = str(init_tuple[7])

        # the following attributes will be filled afterward

        self.matching = list()  # a list that will contain the matching status : if the input coordinates given by
        # the user match with one exon the matching status is : "single", otherwise it's multiple
        # the list will also be filled by the input_coverage_on_exon and the exon_coverage_on_input
        self.gene_name = str()
        self.genomic_sequence = str()
        self.cds_sequence = str()
        self.offsets = list()
        self.peptide_sequence = str()
        self.codon = list()
        self.anticodon = list()
        self.amino_acid = list()
        self.nature = list()
        self.importance = list()  # says if a codon code for an essential, conditionally essential or non essential aa
        self.metabolism = list()  # gives the origin of each amino acid of the sequence
        self.possible_codon = list()  # codons that codes for the same amino acid

    def find_matching_status(self, counter):
        """
        :param counter: an integer corresponding to the number exon retrieved in a given input interval
        set the matching status to "single" if only one exon is found on the input interval, or multiple otherwise
        """
        if counter > 1:
            self.matching.append("multiple")
        else:
            self.matching.append("single")

    def retrieve_gene_name(self, cnx):
        """
        :param cnx: the information necessary to connect to fasterDB
        add the attribute "gene_name" of the exon
        """
        cursor = cnx.cursor()
        query = ("""SELECT official_symbol FROM IMPORT_FasterDB_genes WHERE
                    id = """ + str(self.gene_id) + """ ; """)
        cursor.execute(query)
        counter = 0
        for name in cursor:
            counter += 1
            self.gene_name = str(name[0])
        if counter > 1:
            print "Multiple gene_name were received from an unique gene_id"
            print "So died..."
            raise NameError("MultipleGeneNameError")

    def retrieve_exon_sequences(self, cnx):
        """
        :param cnx: the information necessary to connect to fasterDB
        fill the exon object with the genomic sequences, the cds_sequence and the peptide_sequence of the exon
        and their offsets
        """
        cursor = cnx.cursor()
        query = ("""SELECT DISTINCT genomic_sequence, cds_sequence, offset_before_exon, offset_after_exon FROM
        hsapiens_exonpeptides_filtered WHERE gene_id=""" + self.gene_id + """
        AND exon_position_on_gene=""" + self.exon_number + """ ;""")
        cursor.execute(query)
        counter = 0
        for sequences in cursor:
            counter += 1
            self.genomic_sequence = str(sequences[0])
            self.cds_sequence = str(sequences[1])
            self.offsets.append(int(sequences[2]))
            self.offsets.append(int(sequences[3]))
        if counter > 1:
            print "Multiple sequences were received from an unique gene_id"
            print "So died..."
            raise NameError("MultipleSequencesError")

    def calculate_exon_coverage_on_input_and_input_coverage_on_exon(self, input_chromosomal_coordinates_tuple):
        """
        :param input_chromosomal_coordinates_tuple: a tuple containing the chromosomal coordinates of the exon
        fill the attribute matching with the input_coverage_on_exon and with the exon_coverage_on_input
        # resize the genomic sequence of the exon
        """
        if input_chromosomal_coordinates_tuple[0] > self.exon_start:
            start = input_chromosomal_coordinates_tuple[0]
        else:
            start = self.exon_start
        if input_chromosomal_coordinates_tuple[1] < self.exon_end:
            end = input_chromosomal_coordinates_tuple[1]
        else:
            end = self.exon_end
        # input coverage on exon
        self.matching.append(round((float(end - start + 1) / (self.exon_end - self.exon_start + 1)) * 100, 2))
        # exon coverage on input
        self.matching.append(round(
            (float(end - start + 1) / (
                input_chromosomal_coordinates_tuple[1] - input_chromosomal_coordinates_tuple[0] + 1)
             ) * 100, 2))
        if self.strand == "1":
            self.genomic_sequence = self.genomic_sequence[
                                    start - self.exon_start:len(self.genomic_sequence) - (self.exon_end - end)]
        else:
            self.genomic_sequence = self.genomic_sequence[
                                    self.exon_end - end:len(self.genomic_sequence) - (start - self.exon_start)]

    def resize_cds_sequence(self, input_chromosomal_coordinates_tuple):
        """
        :param input_chromosomal_coordinates_tuple: a tuple containing the chromosomal coordinates of the exon
        resize the cds
        """
        if input_chromosomal_coordinates_tuple[0] > self.cds_start:
            start = input_chromosomal_coordinates_tuple[0]
        else:
            start = self.cds_start
        if input_chromosomal_coordinates_tuple[1] < self.cds_end:
            end = input_chromosomal_coordinates_tuple[1]
        else:
            end = self.cds_end

        if self.strand == "1":
            self.cds_sequence = self.cds_sequence[start - self.cds_start:len(self.cds_sequence) - (self.cds_end - end)]
            new_offset_before_exon = (start - self.cds_start + self.offsets[0]) % 3

        else:
            self.cds_sequence = self.cds_sequence[self.cds_end - end:len(self.cds_sequence) - (start - self.cds_start)]
            new_offset_before_exon = (self.cds_end - end + self.offsets[0]) % 3

        # we have to calculate the new offset before exon if the cds_sequence is cut by the user
        self.offsets[0] = new_offset_before_exon

    def found_codon_anticodon_amino_acid_and_aa_nature(self):
        """
        Fills the list of codons, anticodons, amino_acids and their nature for each exon
        """
        offset_dict = {1: 2, 2: 1, 0: 0}
        for i in range(offset_dict[self.offsets[0]], len(self.cds_sequence), 3):
            if len(self.cds_sequence[i:i + 3]) == 3:
                self.codon.append(self.cds_sequence[i:i + 3])
                self.anticodon.append(codon2anticodon[self.cds_sequence[i:i + 3]])
                self.amino_acid.append(codon2aminoAcid[self.cds_sequence[i:i + 3]])
                self.peptide_sequence += self.amino_acid[-1]
                self.possible_codon.append(amino_acid2codon[self.amino_acid[-1]])

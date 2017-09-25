from dictionnary import codon2aminoAcid
from dictionnary import nature2amino_acid
from dictionnary import amino_acid2nature
from dictionnary import metabolism2amino_acid

class ListExon:
    """
    A class that contains a list of exon.
    Designed to wrap a set of exon and to links it with method
    """

    def __init__(self):
        """
        Initialisation of the attribute exon_list that contains a list of exon
        """
        self.exon_list = list()

    ##############################################################################
    #            getting content of the query_results.xlsx file                  #
    ##############################################################################

    def get_content_mapping_sheet(self):
        """
        :return: A list of list of strings : each sublist corresponds to a row that will be written in the mapping sheet
        of the request_results.xlsx file
        """
        content = [["Name", "Matching_status", "Gene", "Exon", "Chromosome", "Exon_start", "Exon_end",
                    "Input_coverage_on_exon", "Exon_coverage_on_input", "CDS_start", "CDS_end", "Exon_type"]]
        for exon in self.exon_list:
            if not isinstance(exon, list):
                content.append([exon.exon_name, exon.matching[0], exon.gene_name, str(exon.exon_number), str(exon.chr),
                                str(exon.exon_start), str(exon.exon_end), str(exon.matching[1]), str(exon.matching[2]),
                                str(exon.cds_start), str(exon.cds_end), exon.exon_type])
            else:
                content.append([exon[0], exon[1]])
        return content

    def get_content_sequence_sheet(self):
        """
        :return: A list of list of strings : each sublist corresponds to a row that will be written in the sequence
        sheet of the request_results.xlsx file
        """
        content = [["Name", " Genomic_sequence", "CDS_genomic_sequence", "CDS_peptide_sequence", "Reading_frame"]]
        for exon in self.exon_list:
            if not isinstance(exon, list):
                orf = str(exon.codon)
                orf = orf.replace(" ", "").replace("'", "").replace("]", "").replace("[", "").replace(",", "")
                content.append([exon.exon_name, exon.genomic_sequence, exon.cds_sequence, exon.peptide_sequence, orf])
            else:
                content.append([exon[0], exon[1]])
        return content

    def get_content_feature_sheet(self):
        """
        :return: A list of list of strings : each sublist corresponds to a row that will be written in the feature sheet
        of the request_results.xlsx file
        """
        content = []
        j = 0
        while j < len(self.exon_list):
            if not isinstance(self.exon_list[j], list):
                content.append([
                    "Name : " + str(self.exon_list[j].exon_name) +
                    " Gene_Name : " + str(self.exon_list[j].gene_name) + " Exon : " + str(
                        self.exon_list[j].exon_number)])
                row = ["Codon"]
                for i in range(len(self.exon_list[j].codon)):
                    row.append(str(self.exon_list[j].codon[i]))
                content.append(row)
                row = ["Anticodon"]
                for i in range(len(self.exon_list[j].anticodon)):
                    row.append(str(self.exon_list[j].anticodon[i]))
                content.append(row)
                row = ["Amino_acid"]
                for i in range(len(self.exon_list[j].amino_acid)):
                    row.append(str(self.exon_list[j].amino_acid[i]))
                content.append(row)
                row = ["Amino_acid_nature"]
                for i in range(len(self.exon_list[j].nature)):
                    row.append(str(self.exon_list[j].nature[i]))
                content.append(row)
                row = ["Codon_with_same_aa"]
                for i in range(len(self.exon_list[j].possible_codon)):
                    row.append(str(self.exon_list[j].possible_codon[i]))
                content.append(row)
            else:
                content.append(["Name : " + self.exon_list[j][0], self.exon_list[j][1]])
                del (self.exon_list[j])  # deletion of the unwanted elements
                j -= 1
            j += 1
        return content

    def write_a_mapping_file(self, path):
        """
        :param path: the path were the output file will be dropped
        Fills a mapping file.
        """
        file_name = path + "mapping.csv"
        map_file = open(file_name, "w")
        map_file.write(
            "Name;Matching_status;Gene;Exon;Chromosome;Exon_start;Exon_end;Input_coverage_on_exon;" +
            "Exon_coverage_on_input;CDS_start;CDS_end;Exon_type\n")
        for exon in self.exon_list:
            if not isinstance(exon, list):
                map_file.write(exon.exon_name + ";" + exon.matching[0] + ";" + exon.gene_name + ";" + str(
                    exon.exon_number) + ";" + str(exon.chr) + ";" + str(exon.exon_start) + ";" +
                               str(exon.exon_end) + ";" + str(exon.matching[1]) + ";" + str(exon.matching[2]) + ";" +
                               str(exon.cds_start) + ";" + str(exon.cds_end) + ";" + exon.exon_type + "\n")
            else:
                map_file.write("{0};{1}\n".format(exon[0], exon[1]))
        map_file.close()

    def write_a_sequence_file(self, path):
        """
        :param path: the path were the output file will be dropped
        Fills a sequence file containing for each exon, its genomic sequence, its cds sequence, its peptide sequences
        and its open reading frame
        """
        file_name = path + "sequence.csv"
        seq_file = open(file_name, "w")
        seq_file.write("Name;Genomic_sequence;CDS_genomic_sequence;CDS_peptide_sequence;Reading_frame\n")
        for exon in self.exon_list:
            if not isinstance(exon, list):
                orf = str(exon.codon)
                orf = orf.replace(" ", "").replace("'", "").replace("]", "").replace("[", "").replace(",", "")
                seq_file.write(
                    exon.exon_name + ";" + exon.genomic_sequence + ";" + exon.cds_sequence + ";" +
                    exon.peptide_sequence + ";" + orf + "\n")
            else:
                seq_file.write("{0};{1}\n".format(exon[0], exon[1]))
        seq_file.close()

    def write_a_feature_file(self, path):
        """
        :param path: the path were the output file will be dropped
        Fills a feature file : a file containing for each exon its list of codon, anticodon, amino_acid and the nature
        of those amino acid
        """
        file_name = path + "feature.csv"
        feature_file = open(file_name, "w")
        j = 0
        while j < len(self.exon_list):
            if not isinstance(self.exon_list[j], list):
                feature_file.write(
                    "Name : " + str(self.exon_list[j].exon_name) +
                    " Gene_Name : " + str(self.exon_list[j].gene_name) + " Exon : " + str(
                        self.exon_list[j].exon_number) + "\n")
                feature_file.write("Codon")
                for i in range(len(self.exon_list[j].codon)):
                    feature_file.write(";" + str(self.exon_list[j].codon[i]))
                feature_file.write("\nAnticodon")
                for i in range(len(self.exon_list[j].anticodon)):
                    feature_file.write(";" + str(self.exon_list[j].anticodon[i]))
                feature_file.write("\nAmino_acid")
                for i in range(len(self.exon_list[j].amino_acid)):
                    feature_file.write(";" + str(self.exon_list[j].amino_acid[i]))
                feature_file.write("\nAmino_acid_nature")
                for i in range(len(self.exon_list[j].nature)):
                    feature_file.write(";" + str(self.exon_list[j].nature[i]))
                feature_file.write("\ncodon_with_same_aa")
                for i in range(len(self.exon_list[j].possible_codon)):
                    feature_file.write(";" + str(self.exon_list[j].possible_codon[i]))
                feature_file.write("\n")
            else:
                feature_file.write(self.exon_list[j][0] + ";" + self.exon_list[j][1] + "\n")
                del (self.exon_list[j])  # deletion of the unwanted elements
                j -= 1
            j += 1
        feature_file.close()

    def return_the_size_of_each_exon_in_the_list(self):
        """
        :return: the size of each exon from a list of exon. Useful when you want to control the size of the control
        set.
        """
        sizes = list()
        for i in range(len(self.exon_list)):
            sizes.append(len(self.exon_list[i].cds_sequence))
        return sizes

    def weighted_frequency_calculator(self, length_penalty):
        """
        weight the frequency of a codon by the length of the exon and by the number of codon
        for instance the set of exons : ATGCCA  ATGATG  GAGCCA CTG
        will have an ATG frequency of ATG = (0.5 + 1 + 0 + 0)/4 = 0.375
        """
        dic = {}
        for key in codon2aminoAcid.keys():
            dic[key] = 0.
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].codon) > 0:
                    if len(self.exon_list[i].codon) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon)
                    else:
                        c += float(len(self.exon_list[i].codon)) / length_penalty
                        dic[key] += float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon) * \
                                    len(self.exon_list[i].codon) / length_penalty
            dic[key] = dic[key] / c
        return dic

    def weighted_nucleotide_calculator(self, length_penalty):
        """
        :param length_penalty: (int) a value that will reduce the weight of the short codons when the frequency of
        A, T, G, C is calculated
        :return: (dictionary of floats) weight the frequency of a nucleotide by the length  of the exon and by the
         the number of codon.
        """
        dic = {}
        for key in ['A', 'T', 'G', 'C']:
            dic[key] = 0.
            c = 0.
            for i in range(len(self.exon_list)):
                seq = str(self.exon_list[i].codon).replace('\'', "").replace('[', "").replace(' ', '')
                seq = seq.replace(']', '').replace(',', '')
                if len(seq) > 0:
                    if len(seq) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(seq.count(key)) / len(seq)
                    else:
                        c += float(len(seq)) / length_penalty
                        dic[key] += float(seq.count(key)) / len(seq) * len(seq) / length_penalty
            dic[key] = round((dic[key] / c) * 100, 2)
        return dic

    def mean_exon_nucleotide_calculator(self):
        """
        :return: (dictionary of list of floats) : a dictionary that links for each nucleotides a list of floats values.
        Each float value corresponds to the frequency of a particular nucleotide in one of the user's exons.

        """
        dic = {}
        for key in ['A', 'T', 'G', 'C']:
            dic[key] = []
            for i in range(len(self.exon_list)):
                seq = str(self.exon_list[i].codon).replace('\'', "").replace('[', "").replace(' ', '')
                seq = seq.replace(']', '').replace(',', '')
                if len(seq) > 0:
                    dic[key].append(float(seq.count(key)) / len(seq) * 100)
        return dic

    def amino_acid_nature_frequency_calculator(self, length_penalty):
        """
        weight the frequency of the nature of amino acids by the length of the exon and by the number of codons
        """
        dic = {}
        for key in nature2amino_acid.keys():
            dic[key] = 0.
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].amino_acid) > 0:
                    count = 0
                    for aa in nature2amino_acid[key]:
                        count += self.exon_list[i].amino_acid.count(aa)
                    if len(self.exon_list[i].amino_acid) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(count) / len(self.exon_list[i].amino_acid)
                    else:
                        c += float(len(self.exon_list[i].amino_acid)) / length_penalty
                        dic[key] += float(count) / len(self.exon_list[i].amino_acid) * len(
                            self.exon_list[i].amino_acid) / length_penalty
            dic[key] = dic[key] / c
        return dic

    def amino_acid_frequency_calculator(self, length_penalty):
        """
        weight the frequency of an amino acid by the length of the exon and by the number of amino acid
        for instance the set of codon : MG  MM  QR S
        will have an M frequency of M = (0.5 + 1 + 0 + 0)/4 = 0.375
        """
        dic = {}
        for key in amino_acid2nature.keys():
            dic[key] = 0.
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].amino_acid) > 0:
                    if len(self.exon_list[i].amino_acid) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(self.exon_list[i].amino_acid.count(key)) / len(self.exon_list[i].amino_acid)
                    else:
                        c += float(len(self.exon_list[i].amino_acid)) / length_penalty
                        dic[key] += float(self.exon_list[i].amino_acid.count(key)) / len(
                            self.exon_list[i].amino_acid) * len(self.exon_list[i].amino_acid) / length_penalty
            dic[key] = dic[key] / c
        return dic


    def amino_acid_metabolism_frequency_calculator(self, length_penalty):
        """
        weight the frequency of the amino acids according to their origin (TCA_cycle, glycolyse, pentoses)
         (essential, non_essential, conditionally essential) by the length of the exon and by the number of amino acid
        """
        dic = {}
        for key in metabolism2amino_acid.keys():
            dic[key] = 0.
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].amino_acid) > 0:
                    count = 0
                    for aa in metabolism2amino_acid[key]:
                        count += self.exon_list[i].amino_acid.count(aa)
                    if len(self.exon_list[i].amino_acid) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(count) / len(self.exon_list[i].amino_acid)
                    else:
                        c += float(len(self.exon_list[i].amino_acid)) / length_penalty
                        dic[key] += float(count) / len(
                            self.exon_list[i].amino_acid) * len(self.exon_list[i].amino_acid) / length_penalty
            dic[key] = dic[key] / c
        return dic

    def weighted_exon_frequency_calculator(self, length_penalty):
        """
        weight the frequency of a codon by the length of the exon and by the number of codon
        for instance the set of codon : ATGCCA  ATGATG  GAGCCA CTG
        will have an ATG frequency of ATG = (0.5 + 1 + 0 + 0)/4 = 0.375
        """
        dic = {}
        dic_freq_exon = {}
        for key in codon2aminoAcid.keys():
            dic[key] = 0.
            dic_freq_exon[key] = list()
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].codon) > 0:
                    dic_freq_exon[key].append(float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon))
                    if len(self.exon_list[i].codon) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon)
                    else:
                        c += float(len(self.exon_list[i].codon)) / length_penalty
                        dic[key] += float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon) * len(
                            self.exon_list[i].codon) / length_penalty
            dic[key] = dic[key] / c
            for i in range(len(self.exon_list)):
                dic_freq_exon[key].append(len(self.exon_list[i].codon))
            dic_freq_exon[key].append(dic[key])
        return dic_freq_exon

    def amino_acid_exon_frequency_calculator(self, length_penalty):
        """
        weight the frequency of an amino acid by the length of the exon and by the number of amino acid
        for instance the set of codon : MG  MM  QR S
        will have an ATG frequency of ATG = (0.5 + 1 + 0 + 0)/4 = 0.375
        """
        dic = {}
        dic_freq_exon = {}
        for key in amino_acid2nature.keys():
            dic[key] = 0.
            c = 0.
            dic_freq_exon[key] = list()
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].amino_acid) > 0:
                    dic_freq_exon[key].append(
                        float(self.exon_list[i].amino_acid.count(key)) / len(self.exon_list[i].amino_acid))
                    if len(self.exon_list[i].amino_acid) > length_penalty - 1:
                        c += 1.
                        dic[key] += float(self.exon_list[i].amino_acid.count(key)) / len(self.exon_list[i].amino_acid)
                    else:
                        c += float(len(self.exon_list[i].amino_acid)) / length_penalty
                        dic[key] += float(self.exon_list[i].amino_acid.count(key)) / len(
                            self.exon_list[i].amino_acid) * len(self.exon_list[i].amino_acid) / length_penalty
            dic[key] = dic[key] / c
            for i in range(len(self.exon_list)):
                dic_freq_exon[key].append(len(self.exon_list[i].amino_acid))
            dic_freq_exon[key].append(dic[key])
        return dic_freq_exon

    ##################################################################################################################
    #                                              Up_and_down graphic with stats
    ##################################################################################################################

    def codon_frequency_per_exon(self):
        """
        frequency for each codon per exons
        """
        dic_freq_exon = {}
        for key in codon2aminoAcid.keys():
            dic_freq_exon[key] = list()
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].codon) > 9:
                    dic_freq_exon[key].append(float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon))
        return dic_freq_exon

    def amino_acid_frequency_per_exon(self):
        """
        frequency for each amino acid per exons
        """
        dic_freq_exon = {}
        for key in amino_acid2nature.keys():
            dic_freq_exon[key] = list()
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].amino_acid) > 9:
                    dic_freq_exon[key].append(
                        float(self.exon_list[i].amino_acid.count(key)) / len(self.exon_list[i].amino_acid))
        return dic_freq_exon

    ##################################################################################################################
    #                                                   RECAP_GRAPH_FUNCTION
    ##################################################################################################################
    # The methods of this section consists of making boxplot and barplot for every codon and amino acid of the user
    # exon sets
    ##################################################################################################################

    # ------------------------------------------------- Boxplot ------------------------------------------------------

    def get_exon_proportion_of_amino_acid(self, amino_acid):
        """
        :param amino_acid: (string) the amino acid of interest
        :return: (list of floats) the proportion for each exon of the list in the amino acids amino_acid
        """
        amino_acid_prop = list()
        for exon in self.exon_list:
            if len(exon.amino_acid) > 0:
                amino_acid_prop.append(float(exon.amino_acid.count(amino_acid)) /
                                       len(exon.amino_acid) * 100)
        return amino_acid_prop

    def get_exon_proportion_of_codon(self, codon):
        """

        :param codon: (string) the codon of interest
        :return: (list of floats) the proportion for each exon of the list in the codons 'codon'
        """
        codon_prop = list()
        for exon in self.exon_list:
            if len(exon.codon) > 0:
                codon_prop.append(float(exon.codon.count(codon)) /
                                  len(exon.codon) * 100)
        return codon_prop

    # ------------------------------------------------- Barplot --------------------------------------------------------

    def get_exon_proportion_count_of_amino_acid(self, amino_acid):
        """
        :param amino_acid: (string) the amino acid of interest
        :return: (list of float) the proportion of exon that have 0, 1, 2, 3, 4, 5 or more amino acids amino_acid in
        the exon list of interest
        """
        amino_acid_prop_count = [0, 0, 0, 0, 0, 0]
        for exon in self.exon_list:
            for i in range(len(amino_acid_prop_count)):
                if exon.amino_acid.count(amino_acid) == i and i < 5:
                    amino_acid_prop_count[i] += 1
                    break
                elif exon.amino_acid.count(amino_acid) >= i > 4:
                    amino_acid_prop_count[-1] += 1
                    break
        for i in range(len(amino_acid_prop_count)):
            amino_acid_prop_count[i] = float(amino_acid_prop_count[i]) / len(self.exon_list) * 100
        return amino_acid_prop_count

    def get_exon_proportion_count_of_codon(self, codon):
        """
        :param codon: (string) the amino acid of interest
        :return: (list of float) the proportion of exon that have 0, 1, 2, 3, 4, 5 or more codons 'codon' in
        the exon list of interest
        """
        codon_prop_count = [0, 0, 0, 0, 0, 0]
        for exon in self.exon_list:
            for i in range(len(codon_prop_count)):
                if exon.codon.count(codon) == i and i < 5:
                    codon_prop_count[i] += 1
                    break
                elif exon.codon.count(codon) >= i > 4:
                    codon_prop_count[-1] += 1
                    break
        for i in range(len(codon_prop_count)):
            codon_prop_count[i] = float(codon_prop_count[i]) / len(self.exon_list) * 100
        return codon_prop_count

    ##################################################################################################################
    #                                                     RECAP
    ##################################################################################################################
    # This part contains functions that allow the construction of the xlsx file names "resultats.xlsx"
    # This file gave the count and the proportion of every amino acid and codon in the user sets (up, down,up and down)
    # of exons and in a ACE and CCE control sets of exons for comparision.
    # The report also gives the number of A, C, T, G, AT, GC, purine and pyrimidine rich codon (i.e. codon with more
    # than one A, C, T, G, AT, GC, purine and pyrimidine) in the user's exon sets and in the ACE and CCE control set of
    # exon. Finally, it also gives the number and the proportion of A, T, C, and G that occur in the last position of
    # the codons in the user's set of exon and in the control sets for comparision
    #################################################################################################################

    def meta_exon_codon_calculator(self):
        """
        :return: (dictionary of lists [ int, float ]) for each codon (key of the dictionary), gives the number of time
        it occurs in the meta-exon sequence ( i.e. sequence formed by the concatenation of every exon sequence from the
        user/control set of exons) and its proportion in this meta-exon sequence
        """
        dic = {}
        codon_sequence_meta_exon = list()
        for i in range(len(self.exon_list)):
            codon_sequence_meta_exon += self.exon_list[i].codon
        for codon in codon2aminoAcid.keys():
            counter = codon_sequence_meta_exon.count(codon)
            dic[codon] = [float(counter), (float(counter) / len(codon_sequence_meta_exon)) * 100]
        return dic

    def meta_exon_aa_calculator(self):
        """
        :return: (dictionary of lists [ int, float ]) for each amino acid (key of the dictionary), gives the number of
        time it occurs in the meta-exon sequence ( i.e. sequence formed by the concatenation of every exon sequence from
        the user/control set of exons) and its proportion in this meta-exon sequence
        """
        dic = {}
        aa_sequence_meta_exon = list()
        for i in range(len(self.exon_list)):
            aa_sequence_meta_exon += self.exon_list[i].amino_acid
        for amino_acid in amino_acid2nature.keys():
            counter = aa_sequence_meta_exon.count(amino_acid)
            dic[amino_acid] = [float(counter), (float(counter) / len(aa_sequence_meta_exon)) * 100]
        return dic

    def meta_exon_composition_analyser(self):
        """
        :return: (dictionary of lists [ int, float ]) gives the number of A, C, T, G, AC, AG, AT, CG, CT, GT  (key
        of the dictionary) rich codons (i.e codon with more than one A, C, T, G, AT, GC, purine or pyrimidine)
        in the meta-exon sequence ( i.e. sequence formed by the concatenation of every
        exon sequences from the user/control set of exons) and their proportions in this meta-exon sequence
        """
        codon_sequence_meta_exon = list()
        for i in range(len(self.exon_list)):
            codon_sequence_meta_exon += self.exon_list[i].codon
        dic = sequence_all_analyser(codon_sequence_meta_exon,
                                    ['A>=2', 'C>=2', 'G>=2', 'T>=2', 'AC>=2', 'AG>=2', 'AT>=2', 'CG>=2', 'CT>=2',
                                     'GT>=2', 'A=0', 'C=0', 'G=0', 'T=0', 'AC=0', 'AG=0', 'AT=0', 'CG=0', 'CT=0', 'GT=0'
                                     ])
        return dic

    def meta_exon_last_nt_proportion_calculator(self):
        """
        :return: (dictionary of lists [ int, float ]) gives the number of A, C, T, G  located at the third position of
        the codon in the meta-exon sequence ( i.e. sequence formed by the
        concatenation of every exon sequence from the user/control set of exons) and its proportion in this meta-exon
        sequence
        """
        codon_sequence_meta_exon = list()
        dic_nt = {'A': [0., 0.], "C": [0., 0.], "G": [0., 0.], 'T': [0., 0.]}
        for i in range(len(self.exon_list)):
            codon_sequence_meta_exon += self.exon_list[i].codon
        for codon in codon_sequence_meta_exon:
            dic_nt[codon[2]][0] += 1
        for nt in dic_nt.keys():
            dic_nt[nt][1] = (float(dic_nt[nt][0]) / len(codon_sequence_meta_exon)) * 100
        return dic_nt

    def weighted_frequency_calculator_recap(self, length_penalty):
        """
        :param length_penalty: (int) the nucleotide size below which an exon will be considered as small
        :return: (dictionary of lists [ int, float ]) for each codon (key of the dictionary), gives the number of time
        it occurs in  the user/control set of exons and its average proportion in this set. If an exon as a size below
        length_penalty,  the count and the proportion of every codon in the exon will be multiplied by
        size_exon/length_penalty (i.e. method called : weighted method)
        """
        dic = {}
        for key in codon2aminoAcid.keys():
            dic[key] = [0., 0.]
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].codon) > 0:
                    if len(self.exon_list[i].codon) > length_penalty - 1:
                        c += 1.
                        dic[key][0] += float(self.exon_list[i].codon.count(key))
                        dic[key][1] += float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon)
                    else:
                        c += float(len(self.exon_list[i].codon)) / length_penalty
                        dic[key][0] += float(self.exon_list[i].codon.count(key)) * len(
                            self.exon_list[i].codon) / length_penalty
                        dic[key][1] += float(self.exon_list[i].codon.count(key)) / len(self.exon_list[i].codon) * len(
                            self.exon_list[i].codon) / length_penalty
            dic[key][1] = (dic[key][1] / c) * 100
        return dic

    def amino_acid_frequency_calculator_recap(self, length_penalty):
        """
        :param length_penalty: (int) the nucleotide size below which an exon will be considered as small
        :return: (dictionary of lists [ int, float ]) for each amino acid (key of the dictionary), gives the number of
        time it occurs in  the user/control set of exons and its average proportion by exon in this set. If an exon as
        size below length_penalty,  the count and the proportion of every amino acid in the exon will be multiplied by
        size_exon/length_penalty (i.e. method called : weighted method)
        """
        dic = {}
        for key in amino_acid2nature.keys():
            dic[key] = [0., 0.]
            c = 0.
            for i in range(len(self.exon_list)):
                if len(self.exon_list[i].amino_acid) > 0:
                    if len(self.exon_list[i].amino_acid) > length_penalty - 1:
                        c += 1.
                        dic[key][0] += float(self.exon_list[i].amino_acid.count(key))
                        dic[key][1] += float(self.exon_list[i].amino_acid.count(key)) / len(
                            self.exon_list[i].amino_acid)
                    else:
                        c += float(len(self.exon_list[i].amino_acid)) / length_penalty
                        dic[key][0] += float(self.exon_list[i].amino_acid.count(key)) * len(
                            self.exon_list[i].amino_acid) / length_penalty
                        dic[key][1] += float(self.exon_list[i].amino_acid.count(key)) / len(
                            self.exon_list[i].amino_acid) * len(self.exon_list[i].amino_acid) / length_penalty
            dic[key][1] = (dic[key][1] / c) * 100
        return dic

    def weighted_last_nt_frequency_calculator_recap(self, length_penalty):
        """
        :param length_penalty: (int) the nucleotide size below which an exon will be considered as small
        :return: (dictionary of lists [ int, float ]) for each nucleotide (key of the dictionary), gives the number of
        time it occurs in the third position of codons in the user/control set of exons and its average proportion
        by exon in this set. If an exon as a size below length_penalty, the count and the proportion of the last
        nucleotide in the exon will be multiplied by size_exon/length_penalty (i.e. method called : weighted method)
        """
        dic_nt = {'A': [0, 0], "C": [0, 0], "G": [0, 0], 'T': [0, 0]}
        c = {'A': 0., "C": 0., "G": 0., 'T': 0.}
        for i in range(len(self.exon_list)):
            if len(self.exon_list[i].codon) > 0:
                dic = {'A': 0, "C": 0, "G": 0, 'T': 0}
                for codon in self.exon_list[i].codon:
                    dic[codon[2]] += 1
                for nt in dic.keys():
                    if len(self.exon_list[i].codon) > length_penalty - 1:
                        c[nt] += 1.
                        dic_nt[nt][0] += float(dic[nt])
                        dic_nt[nt][1] += float(dic[nt]) / len(self.exon_list[i].codon)
                    else:
                        c[nt] += float(len(self.exon_list[i].codon)) / length_penalty
                        dic_nt[nt][0] += float(dic[nt]) * len(self.exon_list[i].codon) / length_penalty
                        dic_nt[nt][1] += float(dic[nt]) / len(self.exon_list[i].codon) * len(
                            self.exon_list[i].codon) / length_penalty
        for nt in dic_nt.keys():
            dic_nt[nt][1] = (dic_nt[nt][1] / c[nt]) * 100
        return dic_nt

    def weighted_exon_composition_recap(self, length_penalty):
        """
        :param length_penalty: (int) the nucleotide size below which an exon will be considered as small
        :return: (dictionary of lists [ int, float ]) gives the number of codon in the user/control set of exons that
        are rich in A, T, G, C, AC, AG, AT, CG, CT, GT (i.e codon with more than one A, C, T, G, AT, GC, purine
        or pyrimidine)and their proportions by exons in this set. If an exon as  a size below length_penalty, the count
        and the proportion of the last nucleotide in the exon will be multiplied by size_exon/length_penalty
        (i.e. method called : weighted method)
        """
        c = {'A>=2': 0., "C>=2": 0., "G>=2": 0., 'T>=2': 0., 'AC>=2': 0., 'AG>=2': 0., 'AT>=2': 0., 'CG>=2': 0.,
             'CT>=2': 0., 'GT>=2': 0.,
             'A=0': 0., "C=0": 0., "G=0": 0., 'T=0': 0., 'AC=0': 0., 'AG=0': 0., 'AT=0': 0., 'CG=0': 0.,
             'CT=0': 0., 'GT=0': 0.}
        dic_res = {'A>=2': [0., 0.], "C>=2": [0., 0.], "G>=2": [0., 0.], 'T>=2': [0., 0.], 'AC>=2': [0., 0.],
                   'AG>=2': [0., 0.],
                   "AT>=2": [0., 0.], "CG>=2": [0, 0], "CT>=2": [0, 0], "GT>=2": [0, 0],
                   'A=0': [0., 0.], "C=0": [0., 0.], "G=0": [0., 0.], 'T=0': [0., 0.], 'AC=0': [0., 0.],
                   'AG=0': [0., 0.],
                   "AT=0": [0., 0.], "CG=0": [0, 0], "CT=0": [0, 0], "GT=0": [0, 0]
                   }
        for i in range(len(self.exon_list)):
            if len(self.exon_list[i].codon) > 0:
                dic = sequence_all_analyser(self.exon_list[i].codon, dic_res.keys())
                for comp in dic.keys():
                    if len(self.exon_list[i].codon) > length_penalty - 1:
                        c[comp] += 1.
                        dic_res[comp][0] += float(dic[comp][0])
                        dic_res[comp][1] += float(dic[comp][1])
                    else:
                        c[comp] += float(len(self.exon_list[i].codon)) / length_penalty
                        dic_res[comp][0] += float(dic[comp][0]) * float(len(self.exon_list[i].codon)) / length_penalty
                        dic_res[comp][1] += float(dic[comp][1]) * float(len(
                            self.exon_list[i].codon)) / length_penalty
        for comp in dic_res.keys():
            dic_res[comp][1] = (dic_res[comp][1] / c[comp])
        return dic_res

    ########################################################################
    #
    #       graphics in the folder nucleotide_information_figures
    #
    ########################################################################

    def codon_rich_exon_analyser(self):
        """
        :return: (dictionary of list of float) give for each exon, their proportion in codons rich  in A,C,G,T,AC, AG,
        AT, CG, CT, GT
        A codon rich in A contains at least 2 A. A codon rich in AC contains at least 2 A or 2 C or one 1 and one C
        """
        dic_res = {'A': [], "C": [], "G": [], 'T': [], 'M': [], 'R': [], 'W': [], 'S': [], 'Y': [], 'K': []}
        for i in range(len(self.exon_list)):
            if len(self.exon_list[i].codon) > 0:
                dic = sequence_rich_analyser(self.exon_list[i].codon, dic_res.keys())
                for comp in dic.keys():
                    dic_res[comp].append(float(dic[comp][1]))
        return dic_res

    def codon_poor_exon_analyser(self):
        """
        :return: (dictionary of list of float) give for each exon, their proportion in codons D, V, H, B (i.e.
        without C, T, G, A respectively)
        A codon poor in A contains no A. A codon poor in AC contains no  A and no C
        """
        dic_res = {'D': [], "V": [], "H": [], 'B': []}
        for i in range(len(self.exon_list)):
            if len(self.exon_list[i].codon) > 0:
                dic = sequence_poor_analyser(self.exon_list[i].codon, dic_res.keys())
                for comp in dic.keys():
                    dic_res[comp].append(float(dic[comp][1]))
        return dic_res

    def codon_plus_exon_analyser(self):
        """
        :return: (dictionary of list of float) give for each exon, their proportion in codons containing only
        Y, R, W, S, K, M  nt
        """
        dic_res = {'Y': [], 'R': [], 'W': [], 'S': [], 'K': [], 'M': []}
        for i in range(len(self.exon_list)):
            if len(self.exon_list[i].codon) > 0:
                dic = sequence_plus_analyser(self.exon_list[i].codon)
                for comp in dic.keys():
                    dic_res[comp].append(float(dic[comp][1]))
        return dic_res

    def codon_last_nt_proportion(self):
        """
        :return: (dictionary of list of float) give for each exon, their proportion in codons that ends with A,C,G,T
        """
        dic_res = {'A': [], "C": [], "G": [], 'T': []}
        for i in range(len(self.exon_list)):
            if len(self.exon_list[i].codon) > 0:
                dic = {'A': 0, "C": 0, "G": 0, 'T': 0}
                for codon in self.exon_list[i].codon:
                    dic[codon[2]] += 1
                for nt in dic.keys():
                    dic_res[nt].append(float(dic[nt]) / len(self.exon_list[i].codon) * 100)
        return dic_res


def sequence_rich_analyser(list_of_codon, list_type_nt):
    """
    :param list_of_codon: (list of string) the list of codon of an exon/ a metaexon ( i.e. sequence formed by the
        concatenation of every exon sequence from the user/control set of exons)
    :param list_type_nt: a list of the possible nucleotide or groups of to nucleotides we want to study i.e A, T, G, C,
    AC, AG, AT, CG, CT, GT
    :return: (dictionary of lists [ int, float ]), the number of codon in the user/control set of exons that
        are rich in A, T, G, C, AT, GC, purine, or pyrimidine (key of the dictionary)
    (i.e codon with more than one A, C, T, G, AT, GC, purine or pyrimidine)
    """
    res_dic = {}
    for type_nt in list_type_nt:
        res_dic[type_nt] = [0, 0]
    for codon in list_of_codon:
        dic_nt = {'A': 0, "C": 0, "G": 0, 'T': 0}
        for letters in codon:
            dic_nt[letters] += 1
        dic_nt['M'] = dic_nt['A'] + dic_nt['C']
        dic_nt['R'] = dic_nt['A'] + dic_nt['G']
        dic_nt['W'] = dic_nt['A'] + dic_nt['T']
        dic_nt['S'] = dic_nt['C'] + dic_nt['G']
        dic_nt['Y'] = dic_nt['C'] + dic_nt['T']
        dic_nt['K'] = dic_nt['G'] + dic_nt['T']
        for type_nt in dic_nt.keys():
            if dic_nt[type_nt] > 1:
                res_dic[type_nt][0] += 1
    for type_nt in res_dic.keys():
        res_dic[type_nt][1] = round((float(res_dic[type_nt][0]) / len(list_of_codon)) * 100, 3)
    return res_dic


def sequence_poor_analyser(list_of_codon, list_type_nt):
    """
    :param list_of_codon: (list of string) the list of codon of an exon/ a metaexon ( i.e. sequence formed by the
        concatenation of every exon sequence from the user/control set of exons)
    :param list_type_nt: a list of the possible nucleotide or groups of to nucleotides we want to study i.e A, T, G, C,
    AC, AG, AT, CG, CT, GT
    :return: (dictionary of lists [ int, float ]), the number of codon in the user/control set of exons that
        are poor in A, T, G, C, AC, AG, AT, CG, CT, GT (key of the dictionary)
    (i.e codon without A, C, T, G, A and C, A and G, A and T, C and G, C and T, G and T)
    """
    res_dic = {}
    nt2res = {'A': 'B', 'C': 'D', 'G': 'H', 'T': 'V'}
    for type_nt in list_type_nt:
        res_dic[type_nt] = [0, 0]
    for codon in list_of_codon:
        dic_nt = {'A': 0, "C": 0, "G": 0, 'T': 0}
        for letters in codon:
            dic_nt[letters] += 1
        for type_nt in dic_nt.keys():
            if dic_nt[type_nt] == 0:
                res_dic[nt2res[type_nt]][0] += 1
    for type_nt in res_dic.keys():
        res_dic[type_nt][1] = round((float(res_dic[type_nt][0]) / len(list_of_codon)) * 100, 3)
    return res_dic


def sequence_plus_analyser(list_of_codon):
    """
    :param list_of_codon: (list of string) the list of codon of an exon/ a metaexon ( i.e. sequence formed by the
        concatenation of every exon sequence from the user/control set of exons)
    :return: (dictionary of lists [ int, float ]), the number of codon in the user/control set of exons that
        only contains Y, R, W, S K, M nucleotides (key of the dictionary)
    """
    res_dic = {}
    nt2res = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['C', 'G'], 'K': ['G', 'T'], 'M': ['C', 'A']}
    for type_nt in nt2res.keys():
        res_dic[type_nt] = [0, 0]
    for codon in list_of_codon:
        dic_nt = {'A': 0, "C": 0, "G": 0, 'T': 0}
        for letters in codon:
            dic_nt[letters] += 1
        for type_nt in nt2res.keys():
            nt_tot = dic_nt[nt2res[type_nt][0]] + dic_nt[nt2res[type_nt][1]]
            if nt_tot == 3:
                res_dic[type_nt][0] += 1
    for type_nt in res_dic.keys():
        res_dic[type_nt][1] = round((float(res_dic[type_nt][0]) / len(list_of_codon)) * 100, 3)
    return res_dic

def sequence_all_analyser(list_of_codon, list_type_nt):
    """
    :param list_of_codon: (list of string) the list of codon of an exon/ a metaexon ( i.e. sequence formed by the
        concatenation of every exon sequence from the user/control set of exons)
    :param list_type_nt: a list of the possible nucleotide or groups of to nucleotides we want to study i.e A, T, G, C,
    AC, AG, AT, CG, CT, GT
    :return: (dictionary of lists [ int, float ]), the number of codon in the user/control set of exons that
        are rich or poor in A, T, G, C, AC, AG, AT, CG, CT, GT (key of the dictionary)
    (codons rich : codons without A, C, T, G, A and C, A and G, A and T, C and G, C and T, G and T)
    (codons poor : codons with more than one A, C, T, G, AT, GC, purine or pyrimidine)
    """
    res_dic = {}
    for type_nt in list_type_nt:
        res_dic[type_nt] = [0, 0]
    for codon in list_of_codon:
        dic_nt = {'A': 0, "C": 0, "G": 0, 'T': 0}
        for letters in codon:
            dic_nt[letters] += 1
        dic_nt['AC'] = dic_nt['A'] + dic_nt['C']
        dic_nt['AG'] = dic_nt['A'] + dic_nt['G']
        dic_nt['AT'] = dic_nt['A'] + dic_nt['T']
        dic_nt['CG'] = dic_nt['C'] + dic_nt['G']
        dic_nt['CT'] = dic_nt['C'] + dic_nt['T']
        dic_nt['GT'] = dic_nt['G'] + dic_nt['T']
        for type_nt in dic_nt.keys():
            if dic_nt[type_nt] > 1:
                res_dic[type_nt + ">=2"][0] += 1
            if dic_nt[type_nt] == 0:
                res_dic[type_nt + "=0"][0] += 1
    for type_nt in res_dic.keys():
        res_dic[type_nt][1] = round((float(res_dic[type_nt][0]) / len(list_of_codon)) * 100, 3)
    return res_dic

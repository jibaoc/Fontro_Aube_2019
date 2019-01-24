#!/usr/bin/python3


import subprocess
import rpy2.robjects
import rpy2.robjects.vectors as v
import pandas as pd
import os


stretch_i = {"SRSF1": ["G", 6, 7], "SRSF2": ["S", 9, 10], "SRSF3": ["C", 6, 7], "TRA2A-B": ["A", 6, 7]}
folder = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/data_fig1c/")
output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/fig1c_result/")

def retrieve_file_per_sf(folder):
    """
    :param folder: (string) the path to a file containing the cds sequence of a set of exons
    :return:
        sf_files : (dictionary of list of string) each SF (key of the dictionary) is linked with up n down count file
    """
    all_file = subprocess.check_output(["find", folder, "-name", "query_results*.xlsx", "-type", "f"],
                                       stderr=subprocess.STDOUT).decode("ascii").split("\n")[:-1]

    sf_files = {}

    for mfile in all_file:
        c = len(folder.split("/")) -1
        sf = mfile.split("/")[c].split("_")[2]
        if sf in stretch_i:
            if sf not in sf_files:
                sf_files[sf] = [mfile]
            else:
                sf_files[sf].append(mfile)

    return sf_files


def read_sequence(excel_file, stretch_len):
    """
    :param stretch_len: (int) the len of the stretch of interest
    :param excel_file: (string) path to a query file given by the tRNA program
    :return: (list of strings) the list of cds sequence
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
    for row in df.itertuples():
        if not isinstance(row.CDS_genomic_sequence, float):
            if len(row.CDS_genomic_sequence) > stretch_len-1:
                cds.append(row.CDS_genomic_sequence)
    return cds


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


def get_stretch(list_of_sequence, nt, stretch_len, stretch_content):
    """

    :param list_of_sequence: (list of string) list of amino acid sequences
    :param nt: (string) the name of the nt of interest
    :param stretch_len: (int) the length of the stretch of interest
    :param unit_type: (string) the name of the unit of interest (aa, feature, nt)
    :param stretch_content: (int) the number of amino acids participating to the feature
    "feature" that needs to be present in the subsequence of length "stretch_len" to
    say that there ise a stretch in the sub-sequence
    :return: 2 int,
        1 - the number of sequence in ``list_of_sequence`` having at least one stretch of ``nt``
        1 - the number of sequence in ``list_of_sequence`` having no stretch of ``nt``
    """
    stretches_count = 0
    no_stretch_count = 0
    for sequence in list_of_sequence:
        nb_stretch = stretch_finder_nt(sequence, nt, stretch_len, stretch_content)
        if nb_stretch >= 1:
            stretches_count += 1
        else:
            no_stretch_count +=1
    return stretches_count, no_stretch_count


def create_vector(list_of_files, nt, stretch_len):
    """

    :param list_of_files: (list of string) list of files containing exons regulated by the same splicing factor
    :param nt: (string) the nt for which we want to find a stretch
    :param stretch_len: (the length of the stretch we want to find
    :return:
        - val_full (list of int) : one if the exon contains at least one stretch of ``nt`` of size ``stretch_len`` 0 else
        - cell_full (list o string) : cell line of the exon
        - reg_full (list of string) the regulation of the exons
    """
    val_full = []
    cell_full = []
    reg_full = []
    for mfile in list_of_files:
        seqs = read_sequence(mfile, stretch_len)
        stretches_count, no_stretch_count = get_stretch(seqs, nt, stretch_len, stretch_len-1)
        val = [1] * stretches_count + [0] * no_stretch_count
        c = len(folder.split("/")) - 1
        cell = mfile.split("/")[c].split("_")[4]
        r = mfile.split("/")[c+1].split("_")[2].replace(".xlsx","")
        cell_line = [cell] * len(val)
        reg = [r] * len(val)
        val_full += val
        cell_full += cell_line
        reg_full += reg
    return val_full, cell_full, reg_full



def pval_getter(val, cell, reg):

    glm = rpy2.robjects.r("""

    function(val, cell, reg){

    data <- as.data.frame(cbind(val, cell, reg))
    data$val <- as.factor(data$val)
    data$reg <- as.factor(data$reg)
    data$cell <- as.factor(data$cell)
    print(data)
    md0 <-glm(val ~ cell, family=binomial("logit"),data=data)
    md1 <-glm(val ~ reg+cell, family=binomial("logit"),data=data)
    print(summary(md0))
    print(summary(md1))
    a <- anova(md1, md0, test="Chisq")
    print(a)
    return(as.numeric(a$"Pr(>Chi)"[2]))
    }
    """)
    return(glm(v.IntVector(val), v.StrVector(cell), v.StrVector(reg)))


def pval_getter2(val, cell, reg):

    glm = rpy2.robjects.r("""

    function(val, cell, reg){

    data <- as.data.frame(cbind(val, cell, reg))
    data$val <- as.factor(data$val)
    data$reg <- as.factor(data$reg)
    data$cell <- as.factor(data$cell)
    print(data)
    md0 <-glm(val ~ 1, family=binomial("logit"),data=data)
    md1 <-glm(val ~ reg, family=binomial("logit"),data=data)
    print(summary(md0))
    print(summary(md1))
    a <- anova(md1, md0, test="Chisq")
    print(a)
    return(as.numeric(a$"Pr(>Chi)"[2]))
    }
    """)
    return(glm(v.IntVector(val), v.StrVector(cell), v.StrVector(reg)))


def main():
    dic_files = retrieve_file_per_sf(folder)
    if not os.path.isdir(output):
        os.mkdir(output)
    with open(output + "glmstat_fig1c.txt", "w") as outfile:
        for key in dic_files:
            for i in range(1, len(stretch_i[key]), 1):
                line = str(key) + "\t" + str(stretch_i[key][i]-1) + "/" + str(stretch_i[key][i]) + "\tnt : " + stretch_i[key][0]
                print(line)
                val_full, cell_full, reg_full = create_vector(dic_files[key], stretch_i[key][0], stretch_i[key][i])
                if key== "TRA2A-B":
                    pval = pval_getter2(val_full, cell_full, reg_full)
                else:
                    pval = pval_getter(val_full, cell_full, reg_full)
                outfile.write(line + "\tpvall : " + str(pval) + "\n")



if __name__ == "__main__":
    main()

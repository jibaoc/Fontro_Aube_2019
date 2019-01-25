#!/usr/bin/python3


import subprocess
import rpy2.robjects
import os
import rpy2.robjects.vectors as v

folder = os.path.realpath(os.path.dirname(__file__)).replace("src", "data/data_fig3c/")
output = os.path.realpath(os.path.dirname(__file__)).replace("src", "result/fig3c_result/")

aa_i = {"SRSF1": "G", "SRSF2":"A", "SRSF3":"P", "TRA2A-B":"K"}

def retrieve_file_per_sf(folder):
    """
    :param folder: (string) the path to a file containing the number of exons having 0, 1, 2, or 3+ given amino acid in their sequence.
    :return:
        sf_files : (dictionary of list of string) each SF (key of the dictionary) is linked with up n down count file
    """
    all_file = subprocess.check_output(["find", folder, "-name", "*_count_aa", "-type", "f"],
                                       stderr=subprocess.STDOUT).decode("ascii").split("\n")[:-1]

    sf_files = {}

    for mfile in all_file:
        c = len(folder.split("/")) -1
        sf = mfile.split("/")[c].split("_")[2]
        if sf in aa_i:
            if sf not in sf_files:
                sf_files[sf] = [mfile]
            else:
                sf_files[sf].append(mfile)

    return sf_files


def file_reader(a_file, pos):
    """

    :param a_file: (string) a file containing the number of exons having 0, 1, 2, or 3+ given amino acid in their sequence.
    :pos: (int) the number of exons having `pos` in the file
    :return:
        - val (list of int) : 1 for each exons having at least `pos` aa of interest 0 else
        - reg = the regulation times the number of exon in the file
    """


    c = len(folder.split("/")) - 1
    sf = a_file.split("/")[c].split("_")[2]
    with open(a_file, "r") as infile:
        line = infile.readline()
        while line.split("\t")[0] != aa_i[sf]:
            line = infile.readline()
        print(line)
        ival = int(line.split("\t")[pos + 1])
        tot = int(line.split("\t")[5]) - ival
        val = [1] * ival + [0] * tot
        cell = [a_file.split("/")[c].split("_")[4]] * len(val)
        creg = a_file.split("/")[-1].split("_")[0]
        reg = [creg] * len(val)
    return val, cell, reg


def calcul_all_vectors(sf_files, pos):
    """

    :param sf_files: (dictionary of list of string) each SF (key of the dictionary) is linked with up n down count file
    :param pos: (int) the number of exons having `pos` in the file
    :return: (dictionary of list of 3 list of int, string and string) each key is a splicing factor and
        - the 1st list of int says if this exons contains ``pos`` aa interest
        - the 2nd list of strings say the cell line of the exons
        - the 3nd list of strings say the regulation (up/down) of the exons
    """
    reg_dic = {}
    for key in sf_files:
        full_val = []
        full_cell = []
        full_reg = []
        for files in sf_files[key]:
            val, cell, reg = file_reader(files, pos)
            full_val += val
            full_cell += cell
            full_reg += reg
        reg_dic[key] = [full_val, full_cell, full_reg]
    return reg_dic


def pval_getter(val, cell, reg):

    glm = rpy2.robjects.r("""


    function(val, cell, reg){

    data <- as.data.frame(cbind(val, cell, reg))
    data$val <- as.factor(data$val)
    data$reg <- as.factor(data$reg)
    data$cell <- as.factor(data$cell)
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
    sf_files = retrieve_file_per_sf(folder)
    if not os.path.isdir(output):
        os.mkdir(output)
    with open(output + "glmstst_fig3c.txt", "w") as outfile:
        for pos in [0, 1, 2, 3]:
            dic_vector = calcul_all_vectors(sf_files, pos)
            for key in dic_vector:
                line = str(key) + "\tAA : " + str(aa_i[key]) + "\tpos : " + str(pos)
                print(line)
                val, cell, reg = dic_vector[key]
                if key == "TRA2A-B":
                    pval = pval_getter2(val, cell, reg)
                else:
                    pval = pval_getter(val, cell, reg)
                outfile.write(line + "\t pval : " + str(pval) + "\n")

if __name__ == "__main__":
    main()





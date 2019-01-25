import rpy2.robjects as robj
import rpy2.robjects.vectors as v
import warnings
from rpy2.rinterface import RRuntimeWarning
import pandas as pd
import argparse
import os
from matplotlib import pyplot as plt


def read_hexanucleotide(excel_file, fasta, p_cor):
    list_hexa = []
    # opening the excel file
    xl = pd.ExcelFile(excel_file)
    df = "NA"
    # opening the sheet of interest
    for sheet in xl.sheet_names:
        if "hexa_info" == sheet or "hexanucleotide" == sheet:
            df = xl.parse(sheet)
    # if the sheet "exon_skipping*" doesn't exist, the faRLine file cannot be used so we stop the program
    if str(df) == "NA":
        print("the sheet names exon_skipping wasn't found")
        print("exiting...")
        exit(1)
    # creating the inputs for the tRNA program

    if p_cor is True:
        if not fasta:
            for row in df.itertuples():
                if row[10].strip() == "+":
                    list_hexa.append(row[1])
        else:
            for row in df.itertuples():
                if row.reg_p_cor.strip() == "+":
                    list_hexa.append(row[1])
    elif p_cor is False:
        if not fasta:
            for row in df.itertuples():
                if row[9].strip() == "+":
                    list_hexa.append(row[1])
        else:
            for row in df.itertuples():
                if row.reg_pval.strip() == "+":
                    list_hexa.append(row[1])
    elif p_cor == "ten":
        dic_hexa = {}
        if not fasta:
            for row in df.itertuples():
                if row[9].strip() == "+":
                    dic_hexa[row[1]] = float(row[2] - row[5]) / row[5]
            sor_dic = sorted(dic_hexa.items(), key=lambda l: l[1], reverse=True)
        else:
            for row in df.itertuples():
                if row.reg_pval.strip() == "+":
                    dic_hexa[row[1]] = float(row[3] - row[4]) / row[4]
            sor_dic = sorted(dic_hexa.items(), key=lambda l: l[1], reverse=True)
        for i in range(min(len(sor_dic), 10)):
            list_hexa.append(sor_dic[i][0])
    else:
        list_hexa = []
        p_val = []
        if not fasta:
            for row in df.itertuples():
                list_hexa.append(row[1])
                p_val.append(float(row[7]))
                if len(p_val) > 2:
                    if p_val[-1] != p_val[-2] and len(list_hexa) > 9:
                        break
        else:
            for row in df.itertuples():
                list_hexa.append(row[1])
                p_val.append(float(row.pvalues))
                if len(p_val) > 2:
                    if p_val[-1] != p_val[-2] and len(list_hexa) > 9:
                        break

    print(list_hexa)
    return list_hexa


def web_logo_creator(sequences, name_file, output):
    """
    :param sequence_list: (tuple of  8 list of strings) - each list in the tuple corresponds to a list of sequence
    :param sequence_name: (list of string) each string identifies on list of sequence in sequence_list
    :param output: (string) the folder where the results will be created
    """
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    weblogo_maker = robj.r("""
    library("ggplot2")
    library("ggseqlogo")

    function(mys_seq, name_file, mytitle){
        s1 = 15
        cs1 = make_col_scheme(chars=c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M'), groups=c('g1','g2','g3','g4','g5', 'g6', 'g7', 'g8', 'g9', 'g10'),cols=c('limegreen','brown1','gold','dodgerblue3','darkorange', "brown1", "limegreen", "dodgerblue3", "darkorchid3", "dodgerblue3"), name='custom1')
        p1 = ggseqlogo(mys_seq,  method = "probability", col_scheme=cs1, namespace = c('A','T','G','C', 'R', 'Y', 'W', 'S', 'K', 'M')) + theme_logo() + theme(axis.title.y=element_text(size=s1+25), legend.position="none")
        p1 = p1 + ggtitle(mytitle) +  theme(plot.title = element_text(hjust = 0.5))

        p1 = p1 + theme(axis.text=element_text(size=s1 + 50), plot.title = element_text(size=s1 + 60))
        #p1 = p1 + ylim(0,1)
        png(file=name_file,height=300 * 2,width=400 * 2 )
        print(p1)
        dev.off()
    }
    """)
    name = output + name_file
    weblogo_maker(v.StrVector(sequences), name, "")


def counting_nt(list_hexa):
    """

    :param list_hexa: (list of string) list of hexanucleotides of interest
    :return: freq : a list of frequencies,
             labels :  the name associated to the frequencies in "freq"
    """
    seq = "".join(list_hexa)
    freq = [float(seq.count("C")) / len(seq),
            float(seq.count("G")) / len(seq),
            float(seq.count("A")) / len(seq),
            float(seq.count("T")) / len(seq)]
    labels = ["C", "G", "A", "T"]
    return freq, labels


def pie_chart(label_nt, freq_nt, name_fig, output):
    """

    :param label_nt: (list of float) the name of each nucleotides
    :param freq_nt: (list of float) list of frequencies of each nucleotides
    :param name_fig: (string) the name of the figure to create
    :param output: (string) the name of the figure of interest
    """
    fig, ax = plt.subplots(figsize=(2,2))
    ax.pie(freq_nt, labels=label_nt, autopct='%1.1f%%', shadow=True)
    plt.savefig(output + str(name_fig) + "_pie.png")
    plt.clf()
    plt.cla()


def pie_chart_maker(excel_file, fasta, name_fig, output, p_cor):
    """

    :param excel_file:  (string) an excel enrichment file
    :param fasta: (boolean) True if the enrichment file comes from a fasta file containing random sequences
    faslse else.
    :param name_fig: (string) the name of the weblogo to create
    :param output: (string) the path where the logo will be created
    :param p_cor: (boolean/string) True if we want to take only the hexanucleotide having
    a corrected p_value below 0.05 False if we want to take only the enriched hexanucleotides
    having a p_value below 0.05. "ten" if we want to take the hexanucleotides having the
    most difference in frequencies
    """
    sequences = read_hexanucleotide(excel_file, fasta, p_cor)
    freq_nt, label_nt = counting_nt(sequences)
    pie_chart(label_nt, freq_nt, name_fig, output)


def web_logo_maker(excel_file, fasta, name_file, output, p_cor):
    """
    Creation of the logo figure
    :param excel_file: (string) an excel enrichment file
    :param fasta: (boolean) True if the enrichment file comes from a fasta file containing random sequences
    faslse else.
    :param output: (string) the path where the logo will be created
    :param name_file: (string) the name of the weblogo to create
    :param p_cor: (boolean/string) True if we want to take only the hexanucleotide having
    a corrected p_value below 0.05 False if we want to take only the enriched hexanucleotides
    having a p_value below 0.05. "ten" if we want to take the hexanucleotides having the
    most difference in frequencies
    """
    sequences = read_hexanucleotide(excel_file, fasta, p_cor)
    web_logo_creator(sequences, name_file, output)


def launcher():
    """
    function that contains a parser to launch the program
    """
    # description on how to use the program
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""create a weblogo of the hexanucleotide enriched in the excel file
    """,
                                     usage='%(prog)s --input input_file.txt [--output an output folder] ')
    # Arguments for the parser

    req_group = parser.add_argument_group("required arguments")
    parser.add_argument('--output', dest='output', help="An output folder",
                        default=".")
    parser.add_argument('--fasta', dest='fasta', help="True is the enrichment file was generated from a random fasta false else",
                        default=False)
    req_group.add_argument('--excel_file', dest='excel_file', help="An excel folder",
                        required=True)
    req_group.add_argument('--name', dest='name', help="the name of the result file",
                        required=True)
    parser.add_argument('--p_cor', dest='p_cor', help="True if we want the corrected p_val False else or 10 if we want"
                                                      "to create a weblogo with 10 enriched hexanucleotides having "
                                                      "the most different frequencies between the control and the interest freq",
                        default=True)
    parser.add_argument('--pie', dest='pie', help="True if you want to create a pie_chart, false else",
                        default=False)

    args = parser.parse_args()  # parsing arguments


    if not os.path.isdir(args.output):
        print("The given path in 'output' doesn't exist !")
        print("fasta file will be created in your current working directory")
        args.output = "./"

    if args.output[-1] != "/":
        args.output += "/"


    if args.fasta == "False":
        args.fasta = False
    if args.fasta == "True":
        args.fasta = True

    if args.p_cor == "False":
        args.p_cor = False
    if args.p_cor == "True":
        args.p_cor = True

    if args.pie == "False":
        args.pie = False
    if args.pie == "True":
        args.pie = True

    if args.p_cor not in [True, False, "ten", "ten_sig"]:
        print("wrong value for the p_cor argument")
        print("Exiting")
        exit(1)

    if args.pie:
        pie_chart_maker(args.excel_file, args.fasta, args.name, args.output, args.p_cor)
    else:
        web_logo_maker(args.excel_file, args.fasta, args.name, args.output, args.p_cor)

if __name__ == "__main__":
    launcher()
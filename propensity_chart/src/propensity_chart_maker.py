"""Summary: this script wil create propensity chart of a set of exons."""

# imports
import pandas as pd

# functions


def exons_reader(excel_file):
    """Description:read a query_result file produce by the tRNA program.

    :param excel_file: (string) path to an excel file.
        it must be produced by the tRNA program
    :return: (list of string) the list of each sequence in the mapping file.
    """
    sequences = []
    xl = pd.ExcelFile(excel_file)
    df = "NA"
    for sheet in xl.sheet_names:
        if "sequence" == sheet:
            df = xl.parse(sheet)
    # if the sheet nt info doesn't exist, we end the program
    if str(df) == "NA":
        print("The sheet sequence was not found in " + str(excel_file))
        print("Terminating...")
        exit(1)
    for row in df.itertuple():
        if isinstance(row.CDS_peptide_sequence, str):
            seq = row.CDS_peptide_sequence.replace("*", "")
            if len(seq) > 0:
                sequences.append(seq)
    return sequences

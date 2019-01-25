# Goal

## `analysis_maker.py` program

From a folder with the following structure :

```sh
folder/
 |
 folder1/
 |  |
 |  query_result_up.xlsx
 |  query_result_down.xlsx
 folder2/
 |  |
 |  query_result_up.xlsx
 |  query_result_down.xlsx
 ...
 foldern/
```

create in a given output_folder all the sub-folders in your given folder (i.e.
folder1, folder2..., foldern here) and create in those sub_folders the propensity charts
for the 2 sets of peptides (contained in query_result_up.xlsx, query_result_down.xlsx)
given in the sequence sheet of query_result_up/down.xlsx.

The files ``query_result_up.xlsx`` or ``query_result_down.xlsx`` are produced by the tRNA program

This program is based on the program `propensity_chart_maker.py`

## `propensity_chart_maker.py` program

From 2 query results files produced by the randomization-enrichment produce figures
of different features for those 2 sets of peptides within the files

# Prerequisites

Those program use python **version 2.7 or 3** and the following modules:

`propensity_chart_maker.py` program:
* `pandas` to read xlsx files
* `numpy` for array manipulation
* `copy` for creating deep copies
* `matplotlib` for figure creation
* `argparse` for parsing argument
* `scipy` for using statistical test
* `math`


 `analysis_maker.py` program
 * `argparse` for parsing argument
 * `os` to interroged the file system
 * `subprocess` to call `propensity_chart_maker.py`

 # Example

```sh
# command executed for the last version of the program
mkdir result/resultv0.3/
python3 src/analysis_maker.py --data_folder data/ --result_folder result/resultv0.3/

# to get some help :
python2 src/analysis_maker.py --help
python2 src/propensity_chart_maker.py --help
```

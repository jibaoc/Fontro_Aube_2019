# Description

The goal of this script is to produce the statistical analysis for the figure 2E.

# Prerequisites

The script `FigE5_stat_calculator.py` works with **python 3** and uses the following modules :
* pandas
* rpy2
* subprocess
* os

To launch the script `FigE5_stat_calculator.py` you must put in the `data` folder the 8 enrichment folders obtained after launching the tRNA program with the inputs corresponding to the repressed exon list of hnRNPL and PTBP1 in the four cell lines studied.
Each enrichment folder must contained the name of the splicing factor that represses the exons analyzed within that folder.

# Launching the program

```sh
python3 src/Fig5E_stat_calculator.py
```
The above command line will create a file named `Fig5E.stat.txt` in your the `result` folder.

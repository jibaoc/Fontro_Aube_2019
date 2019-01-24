# Description

The goal of this script is to perfom the logistic regression analysis of the figure 1C and 3C.

# Prerequisites

To launch the program you must create two subdirectories in your folder ``data/``:
1. A folder ``data_fig1c`` that must contains the ``query_results`` files produced by the tRNA program corresponding to activated and repressed exons by SRSF1, SRSF2, SRSF3 and TRA2 factors.
2. A folder ``data_fig3c`` that must contains folders produced by the tRNA program. Those corresponds to folders named ``recap_graphics`` created by the tRNA program when launched with an input corresponding to the actived or repressed exons by SRSF1, SRSF2, SRSF3 and TRA2.


# Execution

```sh
python3 src/glm_calculator_fig1c.py
python3 src/glm_calculator_fig3C.py
```

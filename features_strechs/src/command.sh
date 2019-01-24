#!/bin/bash

################################################################################
#                        Stretches evalutions
################################################################################

# the stretch are set to [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12]

data_folder="data/SF_regulated_exons"
folder="result/exploration_strethes"
subfolder="${folder}/exploration_strethes_sf_regulated_exons"
project=($(ls ${data_folder}))
for i in ${!project[@]} ; do
  eval "mkdir ${subfolder}/${project[$i]}"
  echo "python src/stretch_evaluator.py --up ${data_folder}/${project[$i]}/query_results_up.xlsx --down ${data_folder}/${project[$i]}/query_results_down.xlsx --fasta False --output ${subfolder}/${project[$i]} --project \"${project[$i]}\""
  eval "python src/stretch_evaluator.py --up ${data_folder}/${project[$i]}/query_results_up.xlsx --down ${data_folder}/${project[$i]}/query_results_down.xlsx --fasta False --output ${subfolder}/${project[$i]} --project \"${project[$i]}\""
done

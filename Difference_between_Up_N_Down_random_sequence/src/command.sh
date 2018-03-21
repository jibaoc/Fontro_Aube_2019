#!/bin/sh


python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit A --freq_high 0.345 --freq_low 0.23 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit C --freq_high 0.29 --freq_low 0.21 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit G --freq_high 0.273 --freq_low 0.238 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit S --freq_high 0.54 --freq_low 0.47 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit Y --freq_high 0.493 --freq_low 0.455 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit M --freq_high 0.53 --freq_low 0.51 --output result/ --iteration 100 --iscub True

python3 src/random_fasta_comparsion_high_N_low.py --type_unit dnt --unit TC --freq_high 0.07 --freq_low 0.05 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit dnt --unit AC --freq_high 0.063 --freq_low 0.051 --output result/ --iteration 100 --iscub True

python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit A --freq_high 0.345 --freq_low 0.23 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit C --freq_high 0.29 --freq_low 0.21 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit G --freq_high 0.273 --freq_low 0.238 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit S --freq_high 0.54 --freq_low 0.47 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit Y --freq_high 0.493 --freq_low 0.455 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit M --freq_high 0.53 --freq_low 0.51 --output result/ --iteration 100 --iscub False

python3 src/random_fasta_comparsion_high_N_low.py --type_unit dnt --unit TC --freq_high 0.07 --freq_low 0.05 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit dnt --unit AC --freq_high 0.063 --freq_low 0.051 --output result/ --iteration 100 --iscub False

python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Hydrophilic#1 --freq_high 0.4 --freq_low 0.26 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Small#2 --freq_high 0.44 --freq_low 0.41 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Disorder_promoting#1 --freq_high 0.52 --freq_low 0.49 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Polar_uncharged#2 --freq_high 0.29 --freq_low 0.25 --output result/ --iteration 100 --iscub False
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Hydroxylic --freq_high 0.19 --freq_low 0.17 --output result/ --iteration 100 --iscub False

python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Hydrophilic#1 --freq_high 0.4 --freq_low 0.26 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Small#2 --freq_high 0.44 --freq_low 0.41 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Disorder_promoting#1 --freq_high 0.52 --freq_low 0.49 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Polar_uncharged#2 --freq_high 0.29 --freq_low 0.25 --output result/ --iteration 100 --iscub True
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Hydroxylic --freq_high 0.19 --freq_low 0.17 --output result/ --iteration 100 --iscub True







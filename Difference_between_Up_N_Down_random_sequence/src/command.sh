

mkdir result/frequency_explorer_feature_4enrichment_in_nt.dnt
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit A --freq_high 0.345 --freq_low 0.23 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub True --type_unit_interest feature,feature --unit_interest Hydrophilic#1,Hydrophobic#1
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit C --freq_high 0.29 --freq_low 0.21 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub True --type_unit_interest feature,feature,feature --unit_interest Polar-uncharged#2,Neutral,Charged#2
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit C --freq_high 0.29 --freq_low 0.21 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub True --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged

python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit S --freq_high 0.53 --freq_low 0.47 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub True --type_unit_interest feature,feature,feature --unit_interest Very-small,Small#2,Large
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit Y --freq_high 0.491 --freq_low 0.46 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub True --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit M --freq_high 0.53 --freq_low 0.51 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub True --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged

python3 src/random_fasta_comparsion_high_N_low.py --type_unit dnt --unit TC --freq_high 0.065 --freq_low 0.055 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub True --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged
python3 src/random_fasta_comparsion_high_N_low.py --type_unit dnt --unit AC --freq_high 0.061 --freq_low 0.051 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub True --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged



python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit A --freq_high 0.345 --freq_low 0.23 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub False --type_unit_interest feature,feature --unit_interest Hydrophilic#1,Hydrophobic#1
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit C --freq_high 0.29 --freq_low 0.21 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub False --type_unit_interest feature,feature,feature --unit_interest Polar-uncharged#2,Neutral,Charged#2
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit C --freq_high 0.29 --freq_low 0.21 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub False --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged

python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit S --freq_high 0.53 --freq_low 0.47 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub False --type_unit_interest feature,feature,feature --unit_interest Very-small,Small#2,Large
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit Y --freq_high 0.491 --freq_low 0.46 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub False --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged
python3 src/random_fasta_comparsion_high_N_low.py --type_unit nt --unit M --freq_high 0.53 --freq_low 0.51 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub False --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged

python3 src/random_fasta_comparsion_high_N_low.py --type_unit dnt --unit TC --freq_high 0.065 --freq_low 0.055 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub False --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged
python3 src/random_fasta_comparsion_high_N_low.py --type_unit dnt --unit AC --freq_high 0.061 --freq_low 0.051 --output result/frequency_explorer_feature_4enrichment_in_nt.dnt/ --iteration 100 --iscub False --type_unit_interest feature,feature --unit_interest Hydroxylic,Negatively-charged



mkdir result/frequency_explorer_nt_4enrichment_in_feature
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Hydrophilic#1 --freq_high 0.4 --freq_low 0.26 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub False  --type_unit_interest nt,dnt --unit_interest A,AA
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Small#2 --freq_high 0.44 --freq_low 0.41 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub False  --type_unit_interest nt,dnt --unit_interest S,GC
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Very-small --freq_high 0.27 --freq_low 0.21 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub False  --type_unit_interest nt,dnt --unit_interest S,GC
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Polar-uncharged#2 --freq_high 0.29 --freq_low 0.25 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub False  --type_unit_interest nt,dnt --unit_interest C,CC
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Neutral --freq_high 0.38 --freq_low 0.31 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub False  --type_unit_interest nt,dnt --unit_interest C,CC
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Hydroxylic --freq_high 0.19 --freq_low 0.17 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub False --type_unit_interest nt,dnt --unit_interest C,CC

python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Hydrophilic#1 --freq_high 0.4 --freq_low 0.26 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub True  --type_unit_interest nt,dnt --unit_interest A,AA
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Small#2 --freq_high 0.44 --freq_low 0.41 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub True  --type_unit_interest nt,dnt --unit_interest S,GC
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Very-small --freq_high 0.27 --freq_low 0.21 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub True  --type_unit_interest nt,dnt --unit_interest S,GC
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Polar-uncharged#2 --freq_high 0.29 --freq_low 0.25 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub True  --type_unit_interest nt,dnt --unit_interest C,CC
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Neutral --freq_high 0.38 --freq_low 0.31 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub True  --type_unit_interest nt,dnt --unit_interest C,CC
python3 src/random_fasta_comparsion_high_N_low.py --type_unit feature --unit Hydroxylic --freq_high 0.19 --freq_low 0.17 --output result/frequency_explorer_nt_4enrichment_in_feature/ --iteration 100 --iscub True --type_unit_interest nt,dnt --unit_interest C,CC



mkdir result/frequency_explorer_nt_4enrichment_in_2features
python3 src/random_fasta_dependant_feature_high_N_low.py --type_unit feature --unit Small#2,Large --freq_high 0.44,0.34 --freq_low 0.41,0.38 --output result/frequency_explorer_nt_4enrichment_in_2features --iteration 100 --type_unit_interest nt,dnt --unit_interest S,GC
python3 src/random_fasta_dependant_feature_high_N_low.py --type_unit feature --unit Very-small,Large --freq_high 0.27,0.34 --freq_low 0.21,0.38 --output result/frequency_explorer_nt_4enrichment_in_2features --iteration 100 --type_unit_interest nt,dnt --unit_interest S,GC
python3 src/random_fasta_dependant_feature_high_N_low.py --type_unit feature --unit Polar-uncharged#2,Charged#2 --freq_high 0.29,0.17 --freq_low 0.25,0.26 --output result/frequency_explorer_nt_4enrichment_in_2features --iteration 100 --type_unit_interest nt,dnt --unit_interest C,CC
python3 src/random_fasta_dependant_feature_high_N_low.py --type_unit feature --unit Neutral,Charged#2 --freq_high 0.38,0.17 --freq_low 0.31,0.26 --output result/frequency_explorer_nt_4enrichment_in_2features --iteration 100 --type_unit_interest nt,dnt --unit_interest C,CC
python3 src/random_fasta_dependant_feature_high_N_low.py --type_unit feature --unit Hydrophilic#1,Hydrophobic#1  --freq_high 0.4,0.33 --freq_low 0.26,0.39 --output result/frequency_explorer_nt_4enrichment_in_2features --iteration 100 --type_unit_interest nt,dnt --unit_interest A,AA

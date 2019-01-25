# Command executed to populate the result folder

> Note :
>
> The `old_version_command.md` contains the result obtained with an old version of the project.
>
> Make sur to be at the good version of the project by executing :
> `git co v0.2`

# Command line executed

## SRSF2

```sh
mkdir result/with_feature_frequency
mkdir result/with_feature_frequency/SRSF2

# for SRSF2 creation of 2 fasta sequences
# 1. the feature choosen to create is small amino acid

# in 02_1_SRSF2_GSE65349_K562 project for the down of exon s, the small amino acid have a frequency of of 0.52

python src/fasta_reverse_generator.py --feature Small --prob 60 --ctrl CCE --output result/with_feature_frequency/SRSF2
# frequency is 0.594866666667

mkdir result/with_feature_frequency/SRSF2/result_small/

python src/make_comparison.py --output result/with_feature_frequency/SRSF2/result_small/ --ctrl CCE --fasta result/with_feature_frequency/SRSF2/CCE_Small_60.fasta --motif [GC][GC][ATGC][G]

# 1. the frequency choosen to create is disorder amino acid
# in 02_1_SRSF2_GSE65349_K562 project for the down of exon s, the small amino acid have a frequency of of 0.55

python src/fasta_reverse_generator.py --feature Disorder_promoting --prob 60 --ctrl CCE --output result/with_feature_frequency/SRSF2
# disorder freq in the file :  0.602566666667

mkdir result/with_feature_frequency/SRSF2/result_disorder/

python src/make_comparison.py --output result/with_feature_frequency/SRSF2/result_disorder/ --ctrl CCE --fasta result/with_feature_frequency/SRSF2/CCE_Disorder_promoting_60.fasta --motif [GC][GC][ATGC][G]

```

## SRSF3

```sh
mkdir result/with_feature_frequency/SRSF3

# for SRSF2 creation of 2 fasta sequences
# 1. the feature choosen to create is uncharged2 amino acid

# in 03_1_SRSF3_GSE52834_GM19238 project for the down of exon s, the uncharged amino acid have a frequency of of 0.377

python src/fasta_reverse_generator.py --feature Polar_uncharged2 --prob 45 --ctrl CCE --output result/with_feature_frequency/SRSF3
# frequency is 0.45

mkdir result/with_feature_frequency/SRSF3/result_uncharged2/

python src/make_comparison.py --output result/with_feature_frequency/SRSF3/result_uncharged2/ --ctrl CCE --fasta result/with_feature_frequency/SRSF3/CCE_Polar_uncharged2_45.fasta --motif X

# 1. the frequency choosen to create is charged amino acid
# in 03_1_SRSF3_GSE52834_GM19238 project for the down of exon s, the charged amino acid have a frequency of  0.20

python src/fasta_reverse_generator.py --feature Charged --prob 15 --ctrl CCE --output result/with_feature_frequency/SRSF3
# charged freq in the file :  0.15

mkdir result/with_feature_frequency/SRSF3/result_charged/

python src/make_comparison.py --output result/with_feature_frequency/SRSF3/result_charged/ --ctrl CCE --fasta result/with_feature_frequency/SRSF3/CCE_Charged_15.fasta --motif X
```

## TRA2

```sh
mkdir result/with_feature_frequency/TRA2

 # for TRA2 creation of 1 fasta sequence
# the feature choosen to create is charged amino acid
# in 08_1_TRA2A-B_GSE59335_MDA-MB-231 project for the down of exon s, the charged amino acid have a frequency of  0.33

python src/fasta_reverse_generator.py --feature Charged --prob 40 --ctrl CCE --output result/with_feature_frequency/TRA2
# charged freq in the file :  0.40

mkdir result/with_feature_frequency/TRA2/result_charged/

python src/make_comparison.py --output result/with_feature_frequency/TRA2/result_charged/ --ctrl CCE --fasta result/with_feature_frequency/TRA2/CCE_Charged_40.fasta --motif AGAA

```

## hnRNPK

```sh
mkdir result/with_feature_frequency/hnRNPK

 # for hnRNPK creation of 1 fasta sequence
# the feature choosen to create is uncharged amino acid
# in 18_1_hnRNPK_GSE52834_GM19238 project for the up of exon s, the uncharged2 amino acid have a frequency of  0.37

python src/fasta_reverse_generator.py --feature Polar_uncharged2 --prob 45 --ctrl CCE --output result/with_feature_frequency/hnRNPK
# uncharged freq in the file :  0.44

mkdir result/with_feature_frequency/hnRNPK/result_uncharged/

python src/make_comparison.py --output result/with_feature_frequency/hnRNPK/result_uncharged/ --ctrl CCE --fasta result/with_feature_frequency/hnRNPK/CCE_Polar_uncharged2_45.fasta --motif CC

```

## PTBP1

```sh
mkdir result/with_feature_frequency/PTBP1

# for PTBP1 creation of 1 fasta sequence
# the feature choosen to create is uncharged amino acid
# in 23_1_PTBP1_ENCSR064DXG_HepG2 project for the up of exon s, the hydroxilic amino acid have a frequency of  0.20

python src/fasta_reverse_generator.py --feature Hydroxylic --prob 25 --ctrl CCE --output result/with_feature_frequency/PTBP1
# Hydroxylic freq in the file :  0.25

mkdir result/with_feature_frequency/PTBP1/result_hydroxilic/

python src/make_comparison.py --output result/with_feature_frequency/PTBP1/result_hydroxilic/ --ctrl CCE --fasta result/with_feature_frequency/PTBP1/CCE_Hydroxylic_25.fasta --motif CT

```

## hnRNPL

```sh
mkdir result/with_feature_frequency/hnRNPL

# for hnRNPL creation of 1 fasta sequence
# the feature choosen to create is uncharged amino acid
# in  project 19_1_hnRNPL_ENCSR563YIS_K562 for the up of exon s, the hydroxilic amino acid have a frequency of  0.18

python src/fasta_reverse_generator.py --feature Hydroxylic --prob 25 --ctrl CCE --output result/with_feature_frequency/hnRNPL
# Hydroxylic freq in the file :  0.25

mkdir result/with_feature_frequency/hnRNPL/result_hydroxilic/

python src/make_comparison.py --output result/with_feature_frequency/hnRNPL/result_hydroxilic/ --ctrl CCE --fasta result/with_feature_frequency/hnRNPL/CCE_Hydroxylic_25.fasta --motif CA

```

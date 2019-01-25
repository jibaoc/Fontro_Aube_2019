
# Executed commands to populate the result/old_analysis folder

**to reproduce those old result make sur to be at the right version of the script by tapping `git co v0.1`**

1 - Creation of a `control_negatif` folder where random sequences CCE generated with the codons frequency of CCE exons will be compared to a control with the codon frequency in CCE exons. We expect here to find no enrichement !
```sh
mkdir result/old_analysis
mkdir result/old_analysis/control_negatif
mkdir result/old_analysis/control_negatif/CCE_comparison

# generating the fasta sequences with the same codon frequency as the one in CCE exons
python2 src/fasta_generator.py --filename "random_CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/old_analysis/control_negatif --ctrl CCE

# Making the compaison with CCE codon frequencies as control :
python src/make_comparison.py --fasta result/old_analysis/old_analysis/control_negatif/random_CCE_sequences.fasta --motif None --ctrl CCE --output result/old_analysis/control_negatif/CCE_comparison/
```
Conclusion : Nothing here is significant **except the hexanucleotides because it is only the frequency of codons that was used to generate the random sequences***.


2 - In the `control_negatif` folder where completely random sequences (generated with the nt frequency of A:0.25 C:0.25 T:0.25 G:0.25) will be compared to a control with completely random frequency.


```sh
mkdir result/old_analysis/control_negatif/random_comparison

# generating the fasta sequences with the same codon frequency as the one in CCE exons
python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/control_negatif

# Making the compaison with RD codon frequencies as control :
python src/make_comparison.py --fasta result/old_analysis/control_negatif/random_sequences.fasta --motif None --ctrl RD --output result/old_analysis/control_negatif/random_comparison/
```

## Analysing SRSF1

Creation of an `SRSF1` folder and an `01_1_SRSF1_ENCSR066VOO_K562` subfolder.

Creation of a random sequence (no control) and fixing the frequency of the di-nucleotide **GG**. We then make an enrichment comparison with a random control.

The frequency of 0.077 for the di-nucleotide GG was choosed because it's the same as the one found in the enricheùment report of the 01_1_SRSF1_ENCSR066VOO_K562 down exons

```sh
mkdir result/old_analysis/SRSF1
mkdir result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562
mkdir result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562 --GG 0.085 # for having 0.08 % of GG

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/random_sequences.fasta --motif GGAGGA --ctrl RD --output result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/random_comparison/
```

Same thing with CCE control :


```sh
mkdir result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562 --dnt GG --freq 0.072 --ctrl CCE # for having 0.08% of GG

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/CCE_sequences.fasta --motif GGAGGA --ctrl CCE --output result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/CCE_comparison/
```

Random sequence VS CCE control

```sh
mkdir  result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/random_sequences.fasta --motif GGAGGA --ctrl CCE --output result/old_analysis/SRSF1/01_1_SRSF1_ENCSR066VOO_K562/random_vs_CCE/
```

Same thing for the 3 other replicat of SRSF1


```sh

# Project 01_2_SRSF1_ENCSR094KBY_HepG2 of SRSF1

mkdir result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2
mkdir result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2 --GG 0.077

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/random_sequences.fasta --motif GGAGGA --ctrl RD --output result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/random_comparison/


mkdir result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2 --dnt GG --freq 0.077 --ctrl CCE

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/CCE_sequences.fasta --motif GGAGGA --ctrl CCE --output result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/CCE_comparison/

mkdir  result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/random_sequences.fasta --motif GGAGGA --ctrl CCE --output result/old_analysis/SRSF1/01_2_SRSF1_ENCSR094KBY_HepG2/random_vs_CCE/


# project 01_3_SRSF1_GSE52834_GM19238 of SRSF1

mkdir result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238
mkdir result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238 --GG 0.074

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/random_sequences.fasta --motif GGAGGA --ctrl RD --output result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/random_comparison/

mkdir result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238 --dnt GG --freq 0.074 --ctrl CCE

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/CCE_sequences.fasta --motif GGAGGA --ctrl CCE --output result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/CCE_comparison/

mkdir  result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/random_sequences.fasta --motif GGAGGA --ctrl CCE --output result/old_analysis/SRSF1/01_3_SRSF1_GSE52834_GM19238/random_vs_CCE/



mkdir result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA
mkdir result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA --GG 0.087

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/random_sequences.fasta --motif GGAGGA --ctrl RD --output result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/random_comparison/


mkdir result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA --dnt GG --freq 0.087 --ctrl CCE

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/CCE_sequences.fasta --motif GGAGGA --ctrl CCE --output result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/CCE_comparison/

mkdir  result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/random_sequences.fasta --motif GGAGGA --ctrl CCE --output result/old_analysis/SRSF1/01_4_SRSF1_GSE26463_HeLA/random_vs_CCE/
```

## Analysing SRSF2

The frequency of 0.26 for the di-nucleotide SS was choosed because it's the same as the one found in the enricheùment report of the 02_1_SRSF2_GSE65349_K562 down exons

```sh
# Project 02_1_SRSF2_GSE65349_K562 of SRSF2

mkdir result/old_analysis/SRSF2
mkdir result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562
mkdir result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562 --GG 0.065 --GC 0.065 --CG 0.065 --CC 0.065

python src/make_comparison.py --fasta result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/random_sequences.fasta --motif [CG][CG][ATGC]G --ctrl RD --output result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/random_comparison/


mkdir result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562 --dnt SS --freq 0.26 --ctrl CCE

python src/make_comparison.py --fasta result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/CCE_sequences.fasta --motif [CG][CG][ATGC]G --ctrl CCE --output result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/CCE_comparison/

mkdir  result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/random_sequences.fasta --motif [CG][CG][ATGC]G --ctrl CCE --output result/old_analysis/SRSF2/02_1_SRSF2_GSE65349_K562/random_vs_CCE/
```

## Analysing SRSF3
The frequency of 0.091 for the di-nucleotide CC was choosed because it's the same as the one found in the enricheùment report of the 03_1_SRSF3_GSE52834_GM19238 down exons

```sh
# Project 03_1_SRSF3_GSE52834_GM19238 of SRSF3

mkdir result/old_analysis/SRSF3
mkdir result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238
mkdir result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238 --CC 0.098 # for having 91 in the sequence

python src/make_comparison.py --fasta result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/random_sequences.fasta --motif NNN --ctrl RD --output result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/random_comparison/


mkdir result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238 --dnt CC --freq 0.085 --ctrl CCE # for having 0.091 in the sequence

python src/make_comparison.py --fasta result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/CCE_sequences.fasta --motif NNN --ctrl CCE --output result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/CCE_comparison/

mkdir  result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/random_sequences.fasta --motif NNN --ctrl CCE --output result/old_analysis/SRSF3/03_1_SRSF3_GSE52834_GM19238/random_vs_CCE/

```

## Analysing TRA2
The frequency of 0.13 for the di-nucleotide AA was choosed because it's the same as the one found in the enricheùment report of the 08_1_TRA2A-B_GSE59335_MDA down exons

```sh
# Project 08_1_TRA2A-B_GSE59335_MDA-MB-231 of TRA2

mkdir result/old_analysis/TRA2
mkdir result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231
mkdir result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231 --AA 0.13

python src/make_comparison.py --fasta result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/random_sequences.fasta --motif AGAA --ctrl RD --output result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/random_comparison/


mkdir result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231 --dnt AA --freq 0.13 --ctrl CCE # for having 0.091 in the sequence

python src/make_comparison.py --fasta result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/CCE_sequences.fasta --motif AGAA --ctrl CCE --output result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/CCE_comparison/

mkdir  result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/random_sequences.fasta --motif AGAA --ctrl CCE --output result/old_analysis/TRA2/08_1_TRA2A-B_GSE59335_MDA-MB-231/random_vs_CCE/
```

## Analysing hnRNPK

The frequency of 0.1 for the di-nucleotide CC was choosed because it's the same as the one found in the enricheùment report of the 18_1_hnRNPK_GSE52834_GM19238 u exons
```sh
# Project 18_1_hnRNPK_GSE52834_GM19238 of hnRNPK

mkdir result/old_analysis/hnRNPK
mkdir result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238
mkdir result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238 --CC 0.11 # for having 0.1 in the sequence

python src/make_comparison.py --fasta result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/random_sequences.fasta --motif CC --ctrl RD --output result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/random_comparison/


mkdir result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238 --dnt CC --freq 0.1 --ctrl CCE

python src/make_comparison.py --fasta result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/CCE_sequences.fasta --motif CC --ctrl CCE --output result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/CCE_comparison/

mkdir  result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/random_sequences.fasta --motif CC --ctrl CCE --output result/old_analysis/hnRNPK/18_1_hnRNPK_GSE52834_GM19238/random_vs_CCE/
```

## Analysing hnRNPH



The frequency of 0.087 for the di-nucleotide GG was choosed because it's the same as the one found in the enricheùment report of the 17_1_hnRNPH1_GSE34996_293T up exons


```sh
# Project 17_1_hnRNPH1_GSE34996_293T of hnRNPH1

mkdir result/old_analysis/hnRNPH1
mkdir result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T
mkdir result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T --GG 0.092 # for having 0.087 in the sequence

python src/make_comparison.py --fasta result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/random_sequences.fasta --motif GG --ctrl RD --output result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/random_comparison/


mkdir result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T --dnt GG --freq 0.080 --ctrl CCE # for having 0.087 in the sequence

python src/make_comparison.py --fasta result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/CCE_sequences.fasta --motif GG --ctrl CCE --output result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/CCE_comparison/

mkdir  result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/random_sequences.fasta --motif GG --ctrl CCE --output result/old_analysis/hnRNPH1/17_1_hnRNPH1_GSE34996_293T/random_vs_CCE/
```

## Analysing hnRNPL

The frequency of 0.088 for the di-nucleotide CA was choosed because it's the same as the one found in the enricheùment report of the 19_1_hnRNPL_ENCSR563YIS_K562 up exons


```sh
# Project 19_1_hnRNPL_ENCSR563YIS_K562 of hnRNPL

mkdir result/old_analysis/hnRNPL
mkdir result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562
mkdir result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562 --CA 0.092 # for having 0.088 in the sequence

python src/make_comparison.py --fasta result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/random_sequences.fasta --motif CA --ctrl RD --output result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/random_comparison/


mkdir result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562 --dnt CA --freq 0.082 --ctrl CCE # for having 0.088 in the sequence

python src/make_comparison.py --fasta result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/CCE_sequences.fasta --motif CA --ctrl CCE --output result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/CCE_comparison/

mkdir  result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/random_sequences.fasta --motif CA --ctrl CCE --output result/old_analysis/hnRNPL/19_1_hnRNPL_ENCSR563YIS_K562/random_vs_CCE/
```

## Analysing PTBP1

The frequency of 0.083 for the di-nucleotide CT was choosed because it's the same as the one found in the enricheùment report of the 19_1_hnRNPL_ENCSR563YIS_K562 up exons

```sh
# Project 23_1_PTBP1_ENCSR064DXG_HepG2 of PTBP1

mkdir result/old_analysis/PTBP1
mkdir result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2
mkdir result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/random_comparison

python2 src/fasta_generator.py --filename "random_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2 --CT 0.086 # for having 0.83 in the sequence

python src/make_comparison.py --fasta result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/random_sequences.fasta --motif CT --ctrl RD --output result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/random_comparison/


mkdir result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/CCE_comparison

python2 src/fasta_generator.py --filename "CCE_sequences" --nbr_seq 300 --size_inf 50 --size_max 300 --output result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2 --dnt CT --freq 0.078 --ctrl CCE # for having 0.083 in the sequence

python src/make_comparison.py --fasta result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/CCE_sequences.fasta --motif CT --ctrl CCE --output result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/CCE_comparison/

mkdir  result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/random_vs_CCE/

python src/make_comparison.py --fasta result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/random_sequences.fasta --motif CT --ctrl CCE --output result/old_analysis/PTBP1/23_1_PTBP1_ENCSR064DXG_HepG2/random_vs_CCE/
```

--- 

####Finally a recap was made manually from those result in a file named `description_des_resultats.odt`

# Population of the folder /result/version_0.4/

Make sure you're useing the good version of the program
`git co v0.4`

>**The first thing you have to do is to create the file CCE_reading_frame.csv**
>
> To do so:
>1. Open the file `/home/nicolas/Documents/22_distribution_codons_decil_exon/CCE/query_results.xlsx`
>2. Go to the sequence sheet and supresse the genomic_sequence (B), the cds_genomic_sequence (C),  the CDS_peptide_sequence column (D) and only keep the Name(A) and reading_frame column (E).
>3. Eliminate every empty reading frame sequence by sorting the columne 'reading frame' and deleting the line where the reading frame column is empty
>4. Sort by the column name in alphabetical order then save the current sheet as csv with the name 'CCE_reading_frame.csv' in the folder srf of this directory
>5. Delete the header

**Do exactly the same thing with the file `/home/nicolas/Documents/22_distribution_codons_decil_exon/ACE/query_results.xlsx`**

*NB : you can create those files if you want by creating a file with two collumn the first one containing a name : and the second one containing a nucleotide sequence. Those files will be the reference sequence for the program `src/src/fasta_generator_from_real_exons.py`*

> Note : the frequency of feature that will be used for the creation of fasta files is given in the folder /data of this current directory


## control


```sh
mkdir result/version_0.4
mkdir result/version_0.4/control_feature

python src/fasta_reverse_generator.py --output result/version_0.4/control_feature --feature Hydroxylic --prob 16 --filename CCERD_Hydroxilic_16
python src/fasta_reverse_generator.py --output result/version_0.4/control_feature --feature "Hydroxylic,Serine" --prob "16,8" --filename CCERD_Hydroxilic_serin_16_8
python src/fasta_generator_from_real_exons.py --output result/version_0.4/control_feature --feature Hydroxylic --prop 16
python src/fasta_generator_from_real_exons.py --output result/version_0.4/control_feature --feature "Hydroxylic,Serine" --prop "16,7"

python src/fasta_generator_dinucleotide_from_real_exon.py --dnt CC --freq 0.06807 --output result/version_0.4/control_feature --filename CC_0.06807


 mkdir result/version_0.4/control_feature/comparison_rd_exon
 mkdir result/version_0.4/control_feature/comparison_rd_exon_hydro_n_serine
 mkdir result/version_0.4/control_feature/comparison_real_exon
 mkdir result/version_0.4/control_feature/comparison_real_exon_hydro_n_serine
mkdir result/version_0.4/control_feature/comparison_real_exon_dnt


python src/make_comparison.py --output result/version_0.4/control_feature/comparison_rd_exon/ --fasta result/version_0.4/control_feature/CCERD_Hydroxilic_16.fasta --ctrl CCE --motif X

python src/make_comparison.py --output result/version_0.4/control_feature/comparison_real_exon/ --fasta result/version_0.4/control_feature/CCE_Hydroxylic_16.fasta --ctrl CCE --motif X

python src/make_comparison.py --output result/version_0.4/control_feature/comparison_rd_exon_hydro_n_serine/ --fasta result/version_0.4/control_feature/CCERD_Hydroxilic_serin_16_8.fasta --ctrl CCE --motif X

python src/make_comparison.py --output result/version_0.4/control_feature/comparison_real_exon_hydro_n_serine/ --fasta result/version_0.4/control_feature/CCE_Hydroxylic,Serine_16,8.fasta --ctrl CCE --motif X

python src/make_comparison.py --output result/version_0.4/control_feature/comparison_real_exon_hydro_n_serine/ --fasta result/version_0.4/control_feature/CCE_Hydroxylic,Serine_16,8.fasta --ctrl CCE --motif X

python src/make_comparison.py --output result/version_0.4/control_feature/comparison_real_exon_dnt/ --fasta result/version_0.4/control_feature/AA_0.078.fasta --ctrl CCE --motif X

```


## I A - Creation of fasta files containing totaly random sequence.

The totaly random sequence means that the sequence where generated using the codon frequency of CCE but the chain of codon is totally random

```sh
mkdir result/version_0.4/fasta_random_sequence_feature_freq

feature=(Small Small Small,Charged Small,Charged Tiny Tiny Tiny,Charged Tiny,Charged Hydrophilic Hydrophilic Hydrophilic,Negatively_charged Hydrophilic,Negatively_charged Polar_uncharged2 Polar_uncharged2 Polar_uncharged2,Hydroxylic Polar_uncharged2,Hydroxylic Hydroxylic Hydroxylic Hydroxylic,Serine Hydroxylic,Serine Hydroxylic Hydroxylic Hydroxylic,Threonine Hydroxylic,Threonine)

proportion=(53 48 53,25 48,27 32 27 32,25 27,27 54 44 54,16 44,10 36 29 36,20 29,16 19 15 19,11 15,7 19 15 19,7 19,5)

let end=${#feature[@]}-1
for i in `seq 0 $end`; do
  cmd="python src/fasta_reverse_generator.py --output result/version_0.4/fasta_random_sequence_feature_freq/ --ctrl CCE --feature ${feature[$i]} --prob ${proportion[$i]}"
  eval $cmd
  echo "----"
done
# frequence of Small amino acids in the file : 0.526333333333
# ----
# frequence of Small amino acids in the file : 0.482133333333
# ----
# frequence of Small amino acids in the file : 0.527083333333
# frequence of Charged amino acids in the file : 0.25142
# ----
# frequence of Small amino acids in the file : 0.477276666667
# frequence of Charged amino acids in the file : 0.270416666667
# ----
# frequence of Tiny amino acids in the file : 0.322733333333
# ----
# frequence of Tiny amino acids in the file : 0.272633333333
# ----
# frequence of Tiny amino acids in the file : 0.323786666667
# frequence of Charged amino acids in the file : 0.25297
# ----
# frequence of Tiny amino acids in the file : 0.269433333333
# frequence of Charged amino acids in the file : 0.26779
# ----
# frequence of Hydrophilic amino acids in the file : 0.5398
# ----
# frequence of Hydrophilic amino acids in the file : 0.4425
# ----
# frequence of Hydrophilic amino acids in the file : 0.53998
# frequence of Negatively_charged amino acids in the file : 0.161196666667
# ----
# frequence of Hydrophilic amino acids in the file : 0.43001
# frequence of Negatively_charged amino acids in the file : 0.101176666667
# ----
# frequence of Polar_uncharged2 amino acids in the file : 0.359366666667
# ----
# frequence of Polar_uncharged2 amino acids in the file : 0.2969
# ----
# frequence of Polar_uncharged2 amino acids in the file : 0.356686666667
# frequence of Hydroxylic amino acids in the file : 0.200793333333
# ----
# frequence of Polar_uncharged2 amino acids in the file : 0.287313333333
# frequence of Hydroxylic amino acids in the file : 0.158276666667
# ----
# frequence of Hydroxylic amino acids in the file : 0.191633333333
# ----
# frequence of Hydroxylic amino acids in the file : 0.146533333333
# ----
# frequence of Hydroxylic amino acids in the file : 0.195986666667
# frequence of Serine amino acids in the file : 0.113203333333
# ----
# frequence of Hydroxylic amino acids in the file : 0.151763333333
# frequence of Serine amino acids in the file : 0.0706366666667
# ----
# frequence of Hydroxylic amino acids in the file : 0.1934
# ----
# frequence of Hydroxylic amino acids in the file : 0.150333333333
# ----
# frequence of Hydroxylic amino acids in the file : 0.190033333333
# frequence of Threonine amino acids in the file : 0.0688133333333
# ----
# frequence of Hydroxylic amino acids in the file : 0.199186666667
# frequence of Threonine amino acids in the file : 0.0516533333333
# ----
```

## I B - Creation of fasta files containing mutated CCE sequence (to mach a wanted feature frequency).

```sh

mkdir result/version_0.4/fasta_CCEmutated_sequence_feature_freq

feature=(Small Small Small,Charged Small,Charged Tiny Tiny Tiny,Charged Tiny,Charged Hydrophilic Hydrophilic Hydrophilic,Negatively_charged Hydrophilic,Negatively_charged Polar_uncharged2 Polar_uncharged2 Polar_uncharged2,Hydroxylic Polar_uncharged2,Hydroxylic Hydroxylic Hydroxylic Hydroxylic,Serine Hydroxylic,Serine Hydroxylic Hydroxylic Hydroxylic,Threonine Hydroxylic,Threonine)

proportion=(53 48 53,25 48,27 32 27 32,25 27,27 54 44 54,16 44,10 36 29 36,20 29,16 19 15 19,11 15,7 19 15 19,7 19,5)

let end=${#feature[@]}-1
for i in `seq 0 $end`; do
  cmd="python src/fasta_generator_from_real_exons.py --output result/version_0.4/fasta_CCEmutated_sequence_feature_freq/ --ctrl CCE --feature ${feature[$i]} --prop ${proportion[$i]}"
  echo $cmd
  echo "----"
done

# frequence of Small amino acids in the file : 0.535349514284
# ----
# frequence of Small amino acids in the file : 0.482460622509
# ----
# frequence of Smallamino acids in the file : 0.525800425151
# frequence of Chargedamino acids in the file : 0.250656482815
# ----
# frequence of Smallamino acids in the file : 0.476822617375
# frequence of Chargedamino acids in the file : 0.270287921038
# ----
# frequence of Tiny amino acids in the file : 0.327340516439
# ----
# frequence of Tiny amino acids in the file : 0.270911271558
# ----
# frequence of Tinyamino acids in the file : 0.318408759775
# frequence of Chargedamino acids in the file : 0.252447539529
# ----
# frequence of Tinyamino acids in the file : 0.270684812671
# frequence of Chargedamino acids in the file : 0.27032751542
# ----
# frequence of Hydrophilic amino acids in the file : 0.5451024955
# ----
# frequence of Hydrophilic amino acids in the file : 0.437037556396
# ----
# frequence of Hydrophilicamino acids in the file : 0.536050128442
# frequence of Negatively_chargedamino acids in the file : 0.159725766013
# ----
# frequence of Hydrophilicamino acids in the file : 0.441203809625
# frequence of Negatively_chargedamino acids in the file : 0.101253584301
# ----
# frequence of Polar_uncharged2 amino acids in the file : 0.368493603519
# ----
# frequence of Polar_uncharged2 amino acids in the file : 0.290990988054
# ----
# frequence of Polar_uncharged2amino acids in the file : 0.355723497982
# frequence of Hydroxylicamino acids in the file : 0.204813622522
# ----
# frequence of Polar_uncharged2amino acids in the file : 0.289826499252
# frequence of Hydroxylicamino acids in the file : 0.167632735071
# ----
# frequence of Hydroxylic amino acids in the file : 0.195781311838
# ----
# frequence of Hydroxylic amino acids in the file : 0.150389231012
# ----
# frequence of Hydroxylicamino acids in the file : 0.186784317389
# frequence of Serineamino acids in the file : 0.108554771484
# ----
# frequence of Hydroxylicamino acids in the file : 0.150648603379
# frequence of Serineamino acids in the file : 0.0702786461171
# ----
# frequence of Hydroxylic amino acids in the file : 0.196556496796
# ----
# frequence of Hydroxylic amino acids in the file : 0.150427903424
# ----
# frequence of Hydroxylicamino acids in the file : 0.187446868054
# frequence of Threonineamino acids in the file : 0.0712322756285
# ----
# frequence of Hydroxylicamino acids in the file : 0.186256111143
# frequence of Threonineamino acids in the file : 0.0518043315365
# ----
```


## II  - Creation of comparison_result

### II A- For random fasta

```sh
mkdir mkdir result/version_0.4/comparison_random_sequence_feature_freq

list_file=($(ls result/version_0.4/fasta_random_sequence_feature_freq))
motif=(AGAA AGAA AGAA AGAA [TC][CT] [TC][CT] [TC][CT] [CT][TC] [CA][AC] [CA][AC] CC CC CC CC [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G)
echo ${list_file[1]:0:-6}
let c=0
for i in ${list_file[@]}; do
  eval "mkdir result/version_0.4/comparison_random_sequence_feature_freq/${i:0:-6}"
  eval "python src/make_comparison.py --output result/version_0.4/comparison_random_sequence_feature_freq/${i:0:-6}/ --fasta result/version_0.4/fasta_random_sequence_feature_freq/$i --ctrl CCE --motif ${motif[$c]}"
  let c=$c+1
done
```

### II B- For CCE mutated fasta

```sh
mkdir mkdir result/version_0.4/comparison_CCEmutated_sequence_feature_freq

list_file=($(ls result/version_0.4/fasta_CCEmutated_sequence_feature_freq))
motif=(AGAA AGAA AGAA AGAA [TC][CT] [TC][CT] [TC][CT] [CT][TC] [CA][AC] [CA][AC] CC CC CC CC [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G [GC][GC][TACG]G)
echo ${list_file[1]:0:-6}
let c=0
for i in ${list_file[@]}; do
  #eval "mkdir result/version_0.4/comparison_CCEmutated_sequence_feature_freq/${i:0:-6}"
  eval "python src/make_comparison.py --output result/version_0.4/comparison_CCEmutated_sequence_feature_freq/${i:0:-6}/ --fasta result/version_0.4/fasta_CCEmutated_sequence_feature_freq/$i --ctrl CCE --motif ${motif[$c]}"
  let c=$c+1
done
```



# III - Creation of fasta files having a particular dinucleotide frequency


## III A - For CCE mutated fasta

```sh

mkdir result/version_0.4/fasta_CCE_mutated_sequence_dinucleotide_freq

dnt=(GC CC CA CT AA)

# real proportion
proportion=(0.08 0.095 0.09 0.083 0.13)



let end=${#dnt[@]}-1
for i in `seq 0 $end`; do
  cmd="python src/fasta_generator_dinucleotide_from_real_exon.py --output result/version_0.4/fasta_CCE_mutated_sequence_dinucleotide_freq/ --dnt ${dnt[$i]} --freq ${proportion[$i]} --filename CCE_mutated_${dnt[$i]}_${proportion[$i]}"
  eval $cmd
done
# AA : 0.0703132604413 | AC : 0.0523766471446 | AG : 0.0828372743101 | AT : 0.0547994185285 | CA : 0.0824730224829 | CC : 0.0683991647823 | CG : 0.0388748386706 | CT : 0.0675903634985 | GA : 0.0741309149572 | GC : 0.0822768118047 | GG : 0.0668121598595 | GT : 0.0429892819965 | TA : 0.0330524245311 | TC : 0.0542682871177 | TT : 0.0506288867248 | TG : 0.0781772431496 |
# proportion in the file :
# AA : 0.0724259356547 | AC : 0.0564161289101 | AG : 0.0709677243718 | AT : 0.0536120285335 | CA : 0.0773159977961 | CC : 0.0998658069519 | CG : 0.0327686687788 | CT : 0.071989841871 | GA : 0.0718992266505 | GC : 0.0636138899411 | GG : 0.0600247922187 | GT : 0.0443050287151 | TA : 0.0322572510519 | TC : 0.0611073011782 | TT : 0.0547341135616 | TG : 0.076696263815 |
# proportion in the file :
# AA : 0.077733329744 | AC : 0.063245188346 | AG : 0.0803969162573 | AT : 0.0566375724645 | CA : 0.0916725959862 | CC : 0.0707905927952 | CG : 0.0272662849312 | CT : 0.0645605605467 | GA : 0.0756322298824 | GC : 0.0645194074423 | GG : 0.0659324981699 | GT : 0.0440544048103 | TA : 0.0329151858778 | TC : 0.0551595925856 | TT : 0.0522960146994 | TG : 0.0771876254612 |
# proportion in the file :
# AA : 0.0663078986555 | AC : 0.0602363498644 | AG : 0.0728810417811 | AT : 0.0558066675053 | CA : 0.075771657555 | CC : 0.068551568958 | CG : 0.0263489858196 | CT : 0.0851450736958 | GA : 0.0720452703709 | GC : 0.0628164710259 | GG : 0.059426543846 | GT : 0.046135455283 | TA : 0.0406710210475 | TC : 0.0640933574091 | TT : 0.0612120977599 | TG : 0.0825505394229 |
# proportion in the file :
# AA : 0.134705334015 | AC : 0.061160351211 | AG : 0.077655186969 | AT : 0.0631862564252 | CA : 0.0815511942939 | CC : 0.0571987799919 | CG : 0.0210073711422 | CT : 0.0611695252406 | GA : 0.0787092385176 | GC : 0.0523835072748 | GG : 0.0562162254178 | GT : 0.0396653668273 | TA : 0.0420042673667 | TC : 0.04998878342 | TT : 0.0508330955203 | TG : 0.0725655163665 |
```


# IV  Creation of comparison result

## For CCE mutated fasta

```sh

mkdir result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq

dnt=(GC CC CA CT AA)

# real proportion
proportion=(0.08 0.095 0.09 0.083 0.13)

list_file=($(ls result/version_0.4/fasta_CCE_mutated_sequence_dinucleotide_freq | sort))
motif=(AGAA CA CC CT [GC][GC][TAGC]G)
echo ${list_file[1]:0:-6}
let c=0
for i in ${list_file[@]}; do
  #eval "mkdir result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/${i:0:(-6)}"
  eval "python src/make_comparison.py --output result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/${i:0:-6}/ --fasta result/version_0.4/fasta_CCE_mutated_sequence_dinucleotide_freq/$i --ctrl CCE --motif ${motif[$c]}"
  let c=$c+1
done
```

## Hexanucleotides - Weblogo

### pval not corrected

Creation of the weblogos related to the hexanucleotide enriched in the random sequence with specific enrichement
here all the enriched hexanucleotide are taken if their pvalue **not corrected is below 0.05**


```sh

mkdir result/version_0.4/weblogo_hexanucleotides/
mkdir result/version_0.4/weblogo_hexanucleotides/weblogo_feature_random
mkdir result/version_0.4/weblogo_hexanucleotides/weblogo_feature_CCEmutated
mkdir result/version_0.4/weblogo_hexanucleotides/weblogo_CCE_mutated_dnt


file_CCE_mutated=($(find result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_mutated=($(ls result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE_mutated[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_hexanucleotides/weblogo_feature_CCEmutated --name ${name_CCE_mutated[$c]}.png --p_cor False"
  let c=$c+1
done


file_CCE=($(find result/version_0.4/comparison_random_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE=($(ls result/version_0.4/comparison_random_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_hexanucleotides/weblogo_feature_random --name ${name_CCE[$c]}.png --p_cor False"
  let c=$c+1
done


file_CCE_dnt=($(find result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_dnt=($(ls result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_CCE_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_hexanucleotides/weblogo_CCE_mutated_dnt --name ${name_CCE_dnt[$c]}.png --p_cor False
  let c=$c+1
done
```

### pval corrected

Creation of the weblogos related to the hexanucleotide enriched in the random sequence with specific enrichement
here all the enriched hexanucleotide are taken if their pvalue **not corrected is below 0.05**


```sh

mkdir result/version_0.4/weblogo_pcor/
mkdir result/version_0.4/weblogo_pcor/weblogo_feature_random
mkdir result/version_0.4/weblogo_pcor/weblogo_feature_CCEmutated
mkdir result/version_0.4/weblogo_pcor/weblogo_CCE_mutated_dnt


file_CCE_mutated=($(find result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_mutated=($(ls result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE_mutated[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_pcor/weblogo_feature_CCEmutated --name ${name_CCE_mutated[$c]}.png --p_cor True"
  let c=$c+1
done


file_CCE=($(find result/version_0.4/comparison_random_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE=($(ls result/version_0.4/comparison_random_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_pcor/weblogo_feature_random --name ${name_CCE[$c]}.png --p_cor True"
  let c=$c+1
done


file_CCE_dnt=($(find result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_dnt=($(ls result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_CCE_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_pcor/weblogo_CCE_mutated_dnt --name ${name_CCE_dnt[$c]}.png --p_cor True
  let c=$c+1
done
```


### 10 most enriched hexnucleotides



```sh
mkdir result/version_0.4/weblogo_hexanucleotides/10_hexnt_sig/
mkdir result/version_0.4/weblogo_hexanucleotides/10_hexnt_sig/weblogo_feature_random
mkdir result/version_0.4/weblogo_hexanucleotides/10_hexnt_sig/weblogo_feature_CCEmutated
mkdir result/version_0.4/weblogo_hexanucleotides/10_hexnt_sig/weblogo_random_dnt
mkdir result/version_0.4/weblogo_hexanucleotides/10_hexnt_sig/weblogo_CCE_mutated_dnt


file_CCE_mutated=($(find result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_mutated=($(ls result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE_mutated[@]}; do
  echo "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_hexanucleotides/10_hexnt_sig/weblogo_feature_CCEmutated --name ${name_CCE_mutated[$c]}.png --p_cor ten_sig"
  let c=$c+1
done


file_CCE=($(find result/version_0.4/comparison_random_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE=($(ls result/version_0.4/comparison_random_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_hexanucleotides/10_hexnt_sig/weblogo_feature_random --name ${name_CCE[$c]}.png --p_cor ten_sig
  let c=$c+1
done

file_CCE_dnt=($(find result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_dnt=($(ls result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_CCE_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/weblogo_hexanucleotides/10_hexnt_sig/weblogo_CCE_mutated_dnt --name ${name_CCE_dnt[$c]}.png --p_cor ten_sig
  let c=$c+1
done

## Hexanucleotides - Pie_Chart

You must execute this command becouse of a bad commit order...
git co v0.5
 

### exons set of interest

```sh
mkdir result/version_0.4/pie_chart/
mkdir result/version_0.4/pie_chart/weblogo_exon_regulated_by_sf

file=(down/01_1_SRSF1_ENCSR066VOO_K562_result_CCE_down_enrichment_report.xlsx down/01_2_SRSF1_ENCSR094KBY_HepG2_result_CCE_down_enrichment_report.xlsx down/01_3_SRSF1_GSE52834_GM19238_result_CCE_down_enrichment_report.xlsx down/01_4_SRSF1_GSE26463_HeLA_result_CCE_down_enrichment_report.xlsx down/02_1_SRSF2_GSE65349_K562_result_CCE_down_enrichment_report.xlsx down/02_2_SRSF2_GSE78705_Huh7_result_CCE_down_enrichment_report.xlsx down/03_1_SRSF3_GSE52834_GM19238_result_CCE_down_enrichment_report.xlsx down/03_2_SRSF3_ENCSR376FGR_HepG2_result_CCE_down_enrichment_report.xlsx down/05_1_SRSF7_ENCSR464ADT_K562_result_CCE_down_enrichment_report.xlsx down/05_2_SRSF7_ENCSR017PRS_HepG2_result_CCE_down_enrichment_report.xlsx down/08_1_TRA2A-B_GSE59335_MDA-MB-231_result_CCE_down_enrichment_report.xlsx up/17_1_hnRNPH1_GSE34996_293T_result_CCE_up_enrichment_report.xlsx  up/18_1_hnRNPK_GSE52834_GM19238_result_CCE_up_enrichment_report.xlsx up/18_2_hnRNPK_ENCSR529JNJ_K562_result_CCE_up_enrichment_report.xlsx up/18_4_hnRNPK_ENCSR853ZJS_HepG2_result_CCE_up_enrichment_report.xlsx up/19_1_hnRNPL_ENCSR563YIS_K562_result_CCE_up_enrichment_report.xlsx up/19_2_hnRNPL_ENCSR155BMF_HepG2_result_CCE_up_enrichment_report.xlsx up/19_3_hnRNPL_GSE72842_LNCaP_result_CCE_up_enrichment_report.xlsx up/19_4_hnRNPL_GSE52834_GM19238_result_CCE_up_enrichment_report.xlsx up/23_1_PTBP1_ENCSR064DXG_HepG2_result_CCE_up_enrichment_report.xlsx up/23_2_PTBP1_GSE59884_293T_result_CCE_up_enrichment_report.xlsx up/23_3_PTBP1_GSE42701_HeLA_result_CCE_up_enrichment_report.xlsx up/23_4_PTBP1_ENCSR239BCO_K562_result_CCE_up_enrichment_report.xlsx)


name=(SRSF1_down_K562 SRSF1_down_HepG2 SRSF1_down_GM19238 SRSF1_down_Hela SRSF2_down_K562 SRSF2_down_Huh7 SRSF3_down_GM19238 SRSF3_down_HepG2 SRSF7_down_K562 SRSF7_down_HepG2 TRA2_down_MDA hnRNPH1_up_293T hnRNPK_up_GM19238 hnRNPK_up_K562 hnRNPK_up_HepG2 hnRNPL_up_K562 hnRNPL_up_HepG2 hnRNPL_up_LNCap hnRNPL_up_GM19238 PTBP1_up_HepG2 PTBP1_up_293T PTBP1_up_Hela PTBP1_up_K562)

let c=0
for i in ${file[@]}; do
  echo "python src/weblogo_maker.py --excel_file /media/nicolas/DD_2/Projects/splicing_factor_analysis_group_enrichement/new_analysis_January_12/exons_regulated_by_sf/CCE/$i --fasta False --output result/version_0.4/pie_chart/weblogo_exon_regulated_by_sf --name ${name[$c]} --pie True"
  let c=$c+1
done
```

# Fabien hexant

Creating weblogo for fabien enriched di-nucleotide. They are in the folder data.

```sh
mkdir result/version_0.4/pie_chart/pie_fab

file_fab=($(find data/Fabien_hexant_enrichement/ -name "*.xls*" -type f | sort))
name=(hnRNPH1 hnRNPK hnRNPL PTBP1 SRSF1 SRSF2 SRSF3 SRSF7 TRA2)

let c=0
for i in ${file_fab[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta False --output result/version_0.4/pie_chart/pie_fab --name ${name[$c]}.png --pie True
  let c=$c+1
done
```



### pval not corrected

Creation of the weblogos related to the hexanucleotide enriched in the random sequence with specific enrichement
here all the enriched hexanucleotide are taken if their pvalue **not corrected is below 0.05**


```sh

mkdir result/version_0.4/pie_chart/uncor_pval/
mkdir result/version_0.4/pie_chart/uncor_pval/weblogo_feature_random
mkdir result/version_0.4/pie_chart/uncor_pval/weblogo_feature_CCEmutated
mkdir result/version_0.4/pie_chart/uncor_pval/weblogo_CCE_mutated_dnt


file_CCE_mutated=($(find result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_mutated=($(ls result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE_mutated[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/uncor_pval/weblogo_feature_CCEmutated --name ${name_CCE_mutated[$c]}.png --p_cor False --pie True"
  let c=$c+1
done


file_CCE=($(find result/version_0.4/comparison_random_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE=($(ls result/version_0.4/comparison_random_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/uncor_pval/weblogo_feature_random --name ${name_CCE[$c]}.png --p_cor False --pie True"
  let c=$c+1
done


file_CCE_dnt=($(find result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_dnt=($(ls result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_CCE_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/uncor_pval/weblogo_CCE_mutated_dnt --name ${name_CCE_dnt[$c]}.png --p_cor False --pie True
  let c=$c+1
done
```

### pval corrected

Creation of the weblogos related to the hexanucleotide enriched in the random sequence with specific enrichement
here all the enriched hexanucleotide are taken if their pvalue **not corrected is below 0.05**


```sh

mkdir result/version_0.4/pie_chart/cor_pval/
mkdir result/version_0.4/pie_chart/cor_pval/weblogo_feature_random
mkdir result/version_0.4/pie_chart/cor_pval/weblogo_feature_CCEmutated
mkdir result/version_0.4/pie_chart/cor_pval/weblogo_CCE_mutated_dnt


file_CCE_mutated=($(find result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_mutated=($(ls result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE_mutated[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/cor_pval/weblogo_feature_CCEmutated --name ${name_CCE_mutated[$c]}.png --p_cor True --pie True"
  let c=$c+1
done


file_CCE=($(find result/version_0.4/comparison_random_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE=($(ls result/version_0.4/comparison_random_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/cor_pval/weblogo_feature_random --name ${name_CCE[$c]}.png --p_cor True --pie True"
  let c=$c+1
done


file_CCE_dnt=($(find result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_dnt=($(ls result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_CCE_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/cor_pval/weblogo_CCE_mutated_dnt --name ${name_CCE_dnt[$c]}.png --p_cor True --pie True
  let c=$c+1
done

```


### 10 most enriched hexnucleotides



```sh
mkdir result/version_0.4/pie_chart/10_hexnt_sig/
mkdir result/version_0.4/pie_chart/10_hexnt_sig/weblogo_feature_random
mkdir result/version_0.4/pie_chart/10_hexnt_sig/weblogo_feature_CCEmutated
mkdir result/version_0.4/pie_chart/10_hexnt_sig/weblogo_random_dnt
mkdir result/version_0.4/pie_chart/10_hexnt_sig/weblogo_CCE_mutated_dnt


file_CCE_mutated=($(find result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_mutated=($(ls result/version_0.4/comparison_CCEmutated_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE_mutated[@]}; do
  eval "python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/10_hexnt_sig/weblogo_feature_CCEmutated --name ${name_CCE_mutated[$c]}.png --p_cor ten_sig --pie True"
  let c=$c+1
done


file_CCE=($(find result/version_0.4/comparison_random_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE=($(ls result/version_0.4/comparison_random_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/10_hexnt_sig/weblogo_feature_random --name ${name_CCE[$c]}.png --p_cor ten_sig --pie True
  let c=$c+1
done

file_CCE_dnt=($(find result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_dnt=($(ls result/version_0.4/comparison_CCE_mutated_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_CCE_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.4/pie_chart/10_hexnt_sig/weblogo_CCE_mutated_dnt --name ${name_CCE_dnt[$c]}.png --p_cor ten_sig --pie True
  let c=$c+1
done

# Population of the folder /result/version_0.3/

Make sure you're useing the good version of the program
`git co v0.3`

>**The first thing you have to do is to create the file CCE_reading_frame.csv**
>
> To do so:
>1. Open the file `/home/nicolas/Documents/22_distribution_codons_decil_exon/CCE/query_results.xlsx`
>2. Go to the sequence sheet and supresse the genomic_sequence (B), the cds_genomic_sequence (C),  the CDS_peptide_sequence column (D) and only keep the Name(A) and reading_frame column (E).
>3. Eliminate every empty reading frame sequence by sorting the columne 'reading frame' and deleting the line where the reading frame column is empty
>4. Sort by the column name in alphabetical order then save the current sheet as csv with the name 'CCE_reading_frame.csv' in the folder srf of this directory
>5. Delete the header

**Do exactly the sanme thing with the file `/home/nicolas/Documents/22_distribution_codons_decil_exon/ACE/query_results.xlsx`**

*NB : you can create those file if you want by creating a file with two collumn the first one containing a name : and the second one containing a nucleotide sequence. Those file will be the reference sequence for the program `src/src/fasta_generator_from_real_exons.py`*

> **Note : the frequency of feature that will be used for the creation of fasta files is given in the folder /data of this current directory


# Command line used

## I  - Creation of fasta files having a particular feature frequency

### I A - Creation of fasta files containing totaly random sequence.

The totaly random sequence means that the sequence where generated using the codon frequency of CCE fasta_generator_from_real_exons but that the codon where chosen one by one according to their frequency in CCE exons.

```sh
mkdir result/version_0.3
mkdir result/version_0.3/fasta_random_sequence_feature_freq

feature=(Small Small Disorder_promoting Disorder_promoting Order_promoting Order_promoting Hydrophobic Hydrophobic Hydrophilic Hydrophilic Polar_uncharged2 Polar_uncharged2 Charged Charged Negatively_charged Negatively_charged Hydroxylic Hydroxylic Serine Serine Threonine Threonine)

proportion=(53 43 54 44 39 32 52 39 55 40 38 20 34 19 16 9 20 11 10 5 7 4)

let end=${#feature[@]}-1
for i in `seq 0 $end`; do
  cmd="python src/fasta_reverse_generator.py --output result/version_0.3/fasta_random_sequence_feature_freq/ --ctrl CCE --feature ${feature[$i]} --prob ${proportion[$i]}"
  eval $cmd
done
# frequence of Small amino acids in the file : 0.5286
# frequence of Small amino acids in the file : 0.431233333333
# frequence of Disorder_promoting amino acids in the file : 0.5395
# frequence of Disorder_promoting amino acids in the file : 0.442166666667
# frequence of Order_promoting amino acids in the file : 0.3881
# frequence of Order_promoting amino acids in the file : 0.318633333333
# frequence of Hydrophobic amino acids in the file : 0.5259
# frequence of Hydrophobic amino acids in the file : 0.3868
# frequence of Hydrophilic amino acids in the file : 0.547933333333
# frequence of Hydrophilic amino acids in the file : 0.4027
# frequence of Polar_uncharged2 amino acids in the file : 0.381966666667
# frequence of Polar_uncharged2 amino acids in the file : 0.204666666667
# frequence of Charged amino acids in the file : 0.340533333333
# frequence of Charged amino acids in the file : 0.191466666667
# frequence of Negatively_charged amino acids in the file : 0.1611
# frequence of Negatively_charged amino acids in the file : 0.0948666666667
# frequence of Hydroxylic amino acids in the file : 0.1981
# frequence of Hydroxylic amino acids in the file : 0.1114
# frequence of Serine amino acids in the file : 0.101666666667
# frequence of Serine amino acids in the file : 0.0475666666667
# frequence of Threonine amino acids in the file : 0.0686
# frequence of Threonine amino acids in the file : 0.0418
```

### I B - Creation of fasta files containing mutated CCE sequence (to mach a wanted feature frequency).
```sh

mkdir result/version_0.3/fasta_CCEmutated_sequence_feature_freq

feature=(Small Small Disorder_promoting Disorder_promoting Order_promoting Order_promoting Hydrophobic Hydrophobic Hydrophilic Hydrophilic Polar_uncharged2 Polar_uncharged2 Charged Charged Negatively_charged Negatively_charged Hydroxylic Hydroxylic Serine Serine Threonine Threonine)

proportion=(0.5302 0.4338 0.5409 0.4426 0.3928 0.3213 0.5288 0.3909 0.5514 0.4076 0.3791 0.2041 0.3445 0.1855 0.1619 0.0872 0.2026 0.1091 0.09674 0.0521 0.0655 0.0393)

let end=${#feature[@]}-1
for i in `seq 0 $end`; do
  cmd="python src/fasta_generator_from_real_exons.py --output result/version_0.3/fasta_CCEmutated_sequence_feature_freq/ --ctrl CCE --feature ${feature[$i]} --prop ${proportion[$i]}"
  eval $cmd
done
# frequence of Small amino acids in the file : 0.535599898622
# frequence of Small amino acids in the file : 0.430081908632
# frequence of Disorder_promoting amino acids in the file : 0.545610389327
# frequence of Disorder_promoting amino acids in the file : 0.436762058137
# frequence of Order_promoting amino acids in the file : 0.394477767959
# frequence of Order_promoting amino acids in the file : 0.314901048271
# frequence of Hydrophobic amino acids in the file : 0.535868929161
# frequence of Hydrophobic amino acids in the file : 0.383945805126
# frequence of Hydrophilic amino acids in the file : 0.559055052037
# frequence of Hydrophilic amino acids in the file : 0.400420626813
# frequence of Polar_uncharged2 amino acids in the file : 0.389265580178
# frequence of Polar_uncharged2 amino acids in the file : 0.195300689537
# frequence of Charged amino acids in the file : 0.352515831896
# frequence of Charged amino acids in the file : 0.177423748629
# frequence of Negatively_charged amino acids in the file : 0.16915936134
# frequence of Negatively_charged amino acids in the file : 0.0830715631585
# frequence of Hydroxylic amino acids in the file : 0.209572686017
# frequence of Hydroxylic amino acids in the file : 0.103380391995
# frequence of Serine amino acids in the file : 0.100294746554
# frequence of Serine amino acids in the file : 0.048264320168
# frequence of Threonine amino acids in the file : 0.0713770057469
# frequence of Threonine amino acids in the file : 0.034804902143
```

## II  - Creation of comparison_result

### II A- For random fasta

```sh
mkdir mkdir result/version_0.3/comparison_random_sequence_feature_freq

list_file=($(ls result/version_0.3/fasta_random_sequence_feature_freq))
motif=(AGAA AGAA [GC][GC][ATGC]G [GC][GC][ATGC]G AGAA AGAA AGAA AGAA [CT][TC] [CT][TC] [CT][TC] [CT][TC] [GC][GC][ATGC]G [GC][GC][ATGC]G [CT][ACT][TA][CAT][ACT] [CT][ACT][TA][CAT][ACT] [CT][TC] [CT][TC] [GC][GC][ATGC]G [GC][GC][ATGC]G [CA][AC] [CA][AC])
echo ${list_file[1]:0:-6}
let c=0
for i in ${list_file[@]}; do
  #echo "mkdir result/version_0.3/comparison_random_sequence_feature_freq/${i:0:-6}"
  eval "python src/make_comparison.py --output result/version_0.3/comparison_random_sequence_feature_freq/${i:0:-6}/ --fasta result/version_0.3/fasta_random_sequence_feature_freq/$i --ctrl CCE --motif ${motif[$c]}"
  let c=$c+1
done
```

### II B- For CCE mutated fasta

```sh
mkdir mkdir result/version_0.3/comparison_CCEmutated_sequence_feature_freq

list_file=($(ls result/version_0.3/fasta_CCEmutated_sequence_feature_freq))
motif=(AGAA AGAA [GC][GC][ATGC]G [GC][GC][ATGC]G AGAA AGAA AGAA AGAA [CT][TC] [CT][TC] [CT][TC] [CT][TC] [GC][GC][ATGC]G [GC][GC][ATGC]G [CT][ACT][TA][CAT][ACT] [CT][ACT][TA][CAT][ACT] [CT][TC] [CT][TC] [GC][GC][ATGC]G [GC][GC][ATGC]G [CA][AC] [CA][AC])
echo ${list_file[1]:0:-6}
let c=0
for i in ${list_file[@]}; do
  #echo "mkdir result/version_0.3/comparison_random_sequence_feature_freq/${i:0:-6}"
  eval "python src/make_comparison.py --output result/version_0.3/comparison_CCEmutated_sequence_feature_freq/${i:0:-6}/ --fasta result/version_0.3/fasta_CCEmutated_sequence_feature_freq/$i --ctrl CCE --motif ${motif[$c]}"
  let c=$c+1
done
```

# III - Creation of fasta files having a particular dinucleotide frequency


## III A - For random fasta

```sh

mkdir result/version_0.3/fasta_random_sequence_dinucleotide_freq

dnt=(GC CC CA CT AA)

# real proportion
proportion_r=(0.08 0.095 0.09 0.083 0.13)

proportion=(0.085 0.105 0.1 0.087 0.14)

let end=${#dnt[@]}-1
for i in `seq 0 $end`; do
  cmd="python src/fasta_generator.py --output result/version_0.3/fasta_random_sequence_dinucleotide_freq/ --${dnt[$i]} ${proportion[$i]} --filename RD_${dnt[$i]}_${proportion_r[$i]}"
  eval $cmd
done
# di-nucleotides proportions
# AA : 0.061 | AC : 0.061 | AG : 0.061 | AT : 0.061 | CA : 0.061 | CC : 0.061 | CG : 0.061 | CT : 0.061 | GA : 0.061 | GC : 0.085 | GG : 0.061 | GT : 0.061 | TA : 0.061 | TC : 0.061 | TT : 0.061 | TG : 0.061 |
# proportion in the file :
# AA : 0.057501633791 | AC : 0.0614595579389 | AG : 0.0660679827712 | AT : 0.0608667045076 | CA : 0.0655731554877 | CC : 0.0542469838672 | CG : 0.06940081423 | CT : 0.0658368052808 | GA : 0.0607866167631 | GC : 0.0813450246928 | GG : 0.0574461531215 | GT : 0.0584842765538 | TA : 0.0619555129688 | TC : 0.0578826972326 | TT : 0.056112439107 | TG : 0.065033641686 |
# di-nucleotides proportions
# AA : 0.0596666666667 | AC : 0.0596666666667 | AG : 0.0596666666667 | AT : 0.0596666666667 | CA : 0.0596666666667 | CC : 0.105 | CG : 0.0596666666667 | CT : 0.0596666666667 | GA : 0.0596666666667 | GC : 0.0596666666667 | GG : 0.0596666666667 | GT : 0.0596666666667 | TA : 0.0596666666667 | TC : 0.0596666666667 | TT : 0.0596666666667 | TG : 0.0596666666667 |
# proportion in the file :
# AA : 0.056081651056 | AC : 0.0614027292039 | AG : 0.0592098623198 | AT : 0.0629443869152 | CA : 0.0626698104061 | CC : 0.0963565486176 | CG : 0.0618721118962 | CT : 0.0616285168151 | GA : 0.0617017930055 | GC : 0.0633930557677 | GG : 0.0542982255189 | GT : 0.0589996820331 | TA : 0.059006290236 | TC : 0.0615837455343 | TT : 0.0557861071279 | TG : 0.0630654835467 |
# di-nucleotides proportions
# AA : 0.06 | AC : 0.06 | AG : 0.06 | AT : 0.06 | CA : 0.1 | CC : 0.06 | CG : 0.06 | CT : 0.06 | GA : 0.06 | GC : 0.06 | GG : 0.06 | GT : 0.06 | TA : 0.06 | TC : 0.06 | TT : 0.06 | TG : 0.06 |
# proportion in the file :
# AA : 0.0544980688884 | AC : 0.0746980124959 | AG : 0.06450746096 | AT : 0.0658321298972 | CA : 0.0917105116164 | CC : 0.0555763615788 | CG : 0.0562456092652 | CT : 0.0585070015935 | GA : 0.0570938799459 | GC : 0.0656650484134 | GG : 0.0552071866579 | GT : 0.0592657182456 | TA : 0.0564626377748 | TC : 0.0657822682101 | TT : 0.0576287220179 | TG : 0.0613193824391 |
# di-nucleotides proportions
# AA : 0.0608666666667 | AC : 0.0608666666667 | AG : 0.0608666666667 | AT : 0.0608666666667 | CA : 0.0608666666667 | CC : 0.0608666666667 | CG : 0.0608666666667 | CT : 0.087 | GA : 0.0608666666667 | GC : 0.0608666666667 | GG : 0.0608666666667 | GT : 0.0608666666667 | TA : 0.0608666666667 | TC : 0.0608666666667 | TT : 0.0608666666667 | TG : 0.0608666666667 |
# proportion in the file :
# AA : 0.0556380381727 | AC : 0.0664894221104 | AG : 0.0608070628319 | AT : 0.0589064234874 | CA : 0.0587112582072 | CC : 0.056890282659 | CG : 0.059643096215 | CT : 0.0850372200161 | GA : 0.0614370139079 | GC : 0.0647977881382 | GG : 0.0545135909655 | GT : 0.0588706746759 | TA : 0.0664735219352 | TC : 0.0719642421913 | TT : 0.0555585159776 | TG : 0.0642618485087 |
# di-nucleotides proportions
# AA : 0.14 | AC : 0.0573333333333 | AG : 0.0573333333333 | AT : 0.0573333333333 | CA : 0.0573333333333 | CC : 0.0573333333333 | CG : 0.0573333333333 | CT : 0.0573333333333 | GA : 0.0573333333333 | GC : 0.0573333333333 | GG : 0.0573333333333 | GT : 0.0573333333333 | TA : 0.0573333333333 | TC : 0.0573333333333 | TT : 0.0573333333333 | TG : 0.0573333333333 |
# proportion in the file :
# AA : 0.130918904505 | AC : 0.0599645267323 | AG : 0.0604563485038 | AT : 0.0612089414797 | CA : 0.0605448352299 | CC : 0.0522358841665 | CG : 0.0575075873954 | CT : 0.0586213589806 | GA : 0.0596090849982 | GC : 0.0595783743659 | GG : 0.0528377128454 | GT : 0.0566356212867 | TA : 0.061196526202 | TC : 0.0575958603391 | TT : 0.0532106437211 | TG : 0.057877789249 |
```


## III B - For CCE fasta

```sh

mkdir result/version_0.3/fasta_CCE_sequence_dinucleotide_freq

dnt=(GC CC CA CT AA)

# real proportion
proportion_r=(0.08 0.095 0.09 0.083 0.13)

proportion=(0.075 0.087 0.083 0.078 0.123)

let end=${#dnt[@]}-1
for i in `seq 0 $end`; do
  cmd="python src/fasta_generator.py --output result/version_0.3/fasta_CCE_sequence_dinucleotide_freq/ --dnt ${dnt[$i]} --freq ${proportion[$i]} --filename RD_${dnt[$i]}_${proportion_r[$i]} --ctrl CCE"
  echo $cmd
done
  # proportion in the file :
  # AA : 0.0712029336115 | AC : 0.0544392819394 | AG : 0.0739719412183 | AT : 0.0540630124559 | CA : 0.0725431899183 | CC : 0.0693247624796 | CG : 0.054418383133 | CT : 0.065767894752 | GA : 0.0717734135182 | GC : 0.0817505010694 | GG : 0.0667182154635 | GT : 0.0468398014839 | TA : 0.0373564303081 | TC : 0.0572358838506 | TT : 0.0506975974649 | TG : 0.0718967573334 |
  # proportion in the file :
  # AA : 0.0722081760988 | AC : 0.0616691773637 | AG : 0.0708183543657 | AT : 0.0521660863447 | CA : 0.0729649612287 | CC : 0.0961988548519 | CG : 0.0485963002379 | CT : 0.0648653547605 | GA : 0.0720365812787 | GC : 0.065058163317 | GG : 0.0618286558181 | GT : 0.0468679531148 | TA : 0.0389457268891 | TC : 0.0597613467211 | TT : 0.051838179747 | TG : 0.0641761278625 |
  # proportion in the file :
  # AA : 0.0738281032095 | AC : 0.0674569821107 | AG : 0.0754122535892 | AT : 0.0581787193615 | CA : 0.089072855154 | CC : 0.0708148484696 | CG : 0.0443576393031 | CT : 0.0613578965044 | GA : 0.0730913766176 | GC : 0.0680770147839 | GG : 0.0628171388665 | GT : 0.0436643848125 | TA : 0.0387599973907 | TC : 0.0597760751418 | TT : 0.0487949605663 | TG : 0.0645397541188 |
  # proportion in the file :
  # AA : 0.0718841486041 | AC : 0.0615412749696 | AG : 0.068824560828 | AT : 0.0528489975599 | CA : 0.0670182323012 | CC : 0.0679842951492 | CG : 0.0442279425877 | CT : 0.0841686902733 | GA : 0.0730084054194 | GC : 0.0666587602453 | GG : 0.0581526291523 | GT : 0.0451708985217 | TA : 0.0435115359722 | TC : 0.0669863127623 | TT : 0.057169317753 | TG : 0.0708439979009 |
  # proportion in the file :
  # AA : 0.129321828885 | AC : 0.0627079621457 | AG : 0.0740682554413 | AT : 0.0594892897575 | CA : 0.072483539515 | CC : 0.0615213525259 | CG : 0.0415086709545 | CT : 0.0579565778324 | GA : 0.0772974287649 | GC : 0.0567726707187 | GG : 0.056139867995 | GT : 0.0428530356666 | TA : 0.045794900443 | TC : 0.052841310665 | TT : 0.0477662955291 | TG : 0.0614770131604 |
```


# IV  Creation of comparison result

## For Random fasta

```sh

mkdir result/version_0.3/comparison_random_sequence_dinucleotide_freq

dnt=(GC CC CA CT AA)

# real proportion
proportion=(0.08 0.095 0.09 0.083 0.13)

list_file=($(ls result/version_0.3/fasta_random_sequence_dinucleotide_freq))
motif=([GC][GC][ATGC]G  CC [CA][AC] [CT][TC] AGAA)
echo ${list_file[1]:0:-6}
let c=0
for i in ${list_file[@]}; do
  #eval "mkdir result/version_0.3/comparison_random_sequence_dinucleotide_freq/${i:0:(-6)}"
  eval "python src/make_comparison.py --output result/version_0.3/comparison_random_sequence_dinucleotide_freq/${i:0:-6}/ --fasta result/version_0.3/fasta_random_sequence_dinucleotide_freq/$i --ctrl RD --motif ${motif[$c]}"
  let c=$c+1
done
```


## For CCE fasta

```sh

mkdir result/version_0.3/comparison_CCE_sequence_dinucleotide_freq

list_file=($(ls result/version_0.3/fasta_CCE_sequence_dinucleotide_freq))
motif=([GC][GC][ATGC]G  CC [CA][AC] [CT][TC] AGAA)
echo ${list_file[1]:0:-6}
let c=0
for i in ${list_file[@]}; do
  #eval "mkdir result/version_0.3/comparison_CCE_sequence_dinucleotide_freq/${i:0:(-6)}"
  eval "python src/make_comparison.py --output result/version_0.3/comparison_CCE_sequence_dinucleotide_freq/${i:0:-6}/ --fasta result/version_0.3/fasta_CCE_sequence_dinucleotide_freq/$i --ctrl CCE --motif ${motif[$c]}"
  let c=$c+1
done
```

## control


```sh
mkdir result/version_0.3/control_feature

eval "python src/fasta_reverse_generator.py --output result/version_0.3/control_feature --feature Serine --prob 8 --filename CCE_Serine_08"
eval "python src/fasta_generator_from_real_exons.py --output result/version_0.3/control_feature --feature Serine --prop 0.08"

#mkdir result/version_0.3/control_feature/comparison_rd_exon
#mkdir result/version_0.3/control_feature/comparison_real_exon


python src/make_comparison.py --output result/version_0.3/control_feature/comparison_rd_exon/ --fasta result/version_0.3/control_feature/CCE_Serine_08.fasta --ctrl CCE --motif X

python src/make_comparison.py --output result/version_0.3/control_feature/comparison_real_exon/ --fasta result/version_0.3/control_feature/CCE_Serine_0.08.fasta --ctrl CCE --motif X

```

# Hexanucleotides


Creation of the weblogos related to the hexanucleotide enriched in the exons regulated by splicing factors of interest

```sh
mkdir result/version_0.3/weblogo_hexanucleotides
mkdir result/version_0.3/weblogo_hexanucleotides/weblogo_exon_regulated_by_sf

file=(down/01_1_SRSF1_ENCSR066VOO_K562_result_CCE_down_enrichment_report.xlsx down/01_2_SRSF1_ENCSR094KBY_HepG2_result_CCE_down_enrichment_report.xlsx down/01_3_SRSF1_GSE52834_GM19238_result_CCE_down_enrichment_report.xlsx down/01_4_SRSF1_GSE26463_HeLA_result_CCE_down_enrichment_report.xlsx down/02_1_SRSF2_GSE65349_K562_result_CCE_down_enrichment_report.xlsx down/02_2_SRSF2_GSE78705_Huh7_result_CCE_down_enrichment_report.xlsx down/03_1_SRSF3_GSE52834_GM19238_result_CCE_down_enrichment_report.xlsx down/03_2_SRSF3_ENCSR376FGR_HepG2_result_CCE_down_enrichment_report.xlsx down/05_1_SRSF7_ENCSR464ADT_K562_result_CCE_down_enrichment_report.xlsx down/05_2_SRSF7_ENCSR017PRS_HepG2_result_CCE_down_enrichment_report.xlsx down/08_1_TRA2A-B_GSE59335_MDA-MB-231_result_CCE_down_enrichment_report.xlsx up/17_1_hnRNPH1_GSE34996_293T_result_CCE_up_enrichment_report.xlsx  up/18_1_hnRNPK_GSE52834_GM19238_result_CCE_up_enrichment_report.xlsx up/18_2_hnRNPK_ENCSR529JNJ_K562_result_CCE_up_enrichment_report.xlsx up/18_4_hnRNPK_ENCSR853ZJS_HepG2_result_CCE_up_enrichment_report.xlsx up/19_1_hnRNPL_ENCSR563YIS_K562_result_CCE_up_enrichment_report.xlsx up/19_2_hnRNPL_ENCSR155BMF_HepG2_result_CCE_up_enrichment_report.xlsx up/19_3_hnRNPL_GSE72842_LNCaP_result_CCE_up_enrichment_report.xlsx up/19_4_hnRNPL_GSE52834_GM19238_result_CCE_up_enrichment_report.xlsx up/23_1_PTBP1_ENCSR064DXG_HepG2_result_CCE_up_enrichment_report.xlsx up/23_2_PTBP1_GSE59884_293T_result_CCE_up_enrichment_report.xlsx up/23_3_PTBP1_GSE42701_HeLA_result_CCE_up_enrichment_report.xlsx up/23_4_PTBP1_ENCSR239BCO_K562_result_CCE_up_enrichment_report.xlsx)


name=(SRSF1_down_K562 SRSF1_down_HepG2 SRSF1_down_GM19238 SRSF1_down_Hela SRSF2_down_K562 SRSF2_down_Huh7 SRSF3_down_GM19238 SRSF3_down_HepG2 SRSF7_down_K562 SRSF7_down_HepG2 TRA2_down_MDA hnRNPH1_up_293T hnRNPK_up_GM19238 hnRNPK_up_K562 hnRNPK_up_HepG2 hnRNPL_up_K562 hnRNPL_up_HepG2 hnRNPL_up_LNCap hnRNPL_up_GM19238 PTBP1_up_HepG2 PTBP1_up_293T PTBP1_up_Hela PTBP1_up_K562)

let c=0
for i in ${file[@]}; do
  echo "python src/weblogo_maker.py --excel_file /media/nicolas/DD_2/Projects/splicing_factor_analysis_group_enrichement/new_analysis_January_12/exons_regulated_by_sf/CCE/$i --fasta False --output result/version_0.3/weblogo_hexanucleotides/weblogo_exon_regulated_by_sf/ --name ${name[$c]}"
  let c=$c+1
done
```
Creation of the weblogos related to the hexanucleotide enriched in the random sequence with specific enrichement
here all the enriched hexanucleotide are taken if their pvalue **not corrected is below 0.05**

```sh
mkdir result/version_0.3/weblogo_hexanucleotides/weblogo_feature_random
mkdir result/version_0.3/weblogo_hexanucleotides/weblogo_feature_CCEmutated
mkdir result/version_0.3/weblogo_hexanucleotides/weblogo_random_dnt
mkdir result/version_0.3/weblogo_hexanucleotides/weblogo_CCE_dnt


file_CCE_mutated=($(find result/version_0.3/comparison_CCEmutated_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_mutated=($(ls result/version_0.3/comparison_CCEmutated_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE_mutated[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.3/weblogo_hexanucleotides/weblogo_feature_CCEmutated --name ${name_CCE_mutated[$c]}.png --p_cor False
  let c=$c+1
done


file_CCE=($(find result/version_0.3/comparison_random_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE=($(ls result/version_0.3/comparison_random_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.3/weblogo_hexanucleotides/weblogo_feature_random --name ${name_CCE[$c]}.png --p_cor False
  let c=$c+1
done

file_rd_dnt=($(find result/version_0.3/comparison_random_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_rd_dnt=($(ls result/version_0.3/comparison_random_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_rd_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.3/weblogo_hexanucleotides/weblogo_random_dnt --name ${name_rd_dnt[$c]}.png --p_cor False
  let c=$c+1
done

file_CCE_dnt=($(find result/version_0.3/comparison_CCE_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_dnt=($(ls result/version_0.3/comparison_CCE_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_CCE_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.3/weblogo_hexanucleotides/weblogo_CCE_dnt --name ${name_CCE_dnt[$c]}.png --p_cor False
  let c=$c+1
done
```

Creating weblogo for fabien enriched di-nucleotide. They are in the folder data.

```sh
mkdir result/version_0.3/weblogo_hexanucleotides/weblogo_Fabien

file_fab=($(find data/Fabien_hexant_enrichement/ -name "*.xls*" -type f | sort))
name=(hnRNPH1 hnRNPK hnRNPL PTBP1 SRSF1 SRSF2 SRSF3 SRSF7 TRA2)

let c=0
for i in ${file_fab[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta False --output result/version_0.3/weblogo_hexanucleotides/weblogo_Fabien/ --name ${name[$c]}.png
  let c=$c+1
done
```

Creation of the weblogos related to the hexanucleotide enriched in the random sequence with specific enrichement
here all the enriched hexanucleotide are taken if their pvalue ** corrected is below 0.05**

```sh
mkdir result/version_0.3/weblogo_hexanucleotides/p_cor/
mkdir result/version_0.3/weblogo_hexanucleotides/p_cor/weblogo_feature_random
mkdir result/version_0.3/weblogo_hexanucleotides/p_cor/weblogo_feature_CCEmutated
mkdir result/version_0.3/weblogo_hexanucleotides/p_cor/weblogo_random_dnt
mkdir result/version_0.3/weblogo_hexanucleotides/p_cor/weblogo_CCE_dnt


file_CCE_mutated=($(find result/version_0.3/comparison_CCEmutated_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_mutated=($(ls result/version_0.3/comparison_CCEmutated_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE_mutated[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.3/weblogo_hexanucleotides/p_cor/weblogo_feature_CCEmutated --name ${name_CCE_mutated[$c]}.png --p_cor True
  let c=$c+1
done


file_CCE=($(find result/version_0.3/comparison_random_sequence_feature_freq/ -name "enrichment_report*" -type f | sort))
name_CCE=($(ls result/version_0.3/comparison_random_sequence_feature_freq/ | sort))

let c=0
for i in ${file_CCE[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.3/weblogo_hexanucleotides/p_cor/weblogo_feature_random --name ${name_CCE[$c]}.png --p_cor True
  let c=$c+1
done

file_rd_dnt=($(find result/version_0.3/comparison_random_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_rd_dnt=($(ls result/version_0.3/comparison_random_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_rd_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.3/weblogo_hexanucleotides/p_cor/weblogo_random_dnt --name ${name_rd_dnt[$c]}.png --p_cor True
  let c=$c+1
done

file_CCE_dnt=($(find result/version_0.3/comparison_CCE_sequence_dinucleotide_freq/ -name "enrichment_report*" -type f | sort))
name_CCE_dnt=($(ls result/version_0.3/comparison_CCE_sequence_dinucleotide_freq/ | sort))

let c=0
for i in ${file_CCE_dnt[@]}; do
  python src/weblogo_maker.py --excel_file $i --fasta True --output result/version_0.3/weblogo_hexanucleotides/p_cor/weblogo_CCE_dnt --name ${name_CCE_dnt[$c]}.png --p_cor True
  let c=$c+1
done
```

#2018 04 06
# JBC
Figures last version
w. Panalyse_NT_AA_freqs.pl
~/Bureau/Didier/ARTICLE_AA_nt/201804_new_new_groups/CSV

for i in *csv; do echo work on $i ; Panalyse_NT_AA_freqs.pl -minlength 0 -in $i |& tee $i.freq.log ; done

#compute Density graphs for "n" UP vs Down
for up in *up.n_freq.txt ; do down=`echo $up | sed s/_up/_down/`; short=`echo $up | sed s/_up//`;  Panalyse_NT_AA_graphs.pl -mode freq -in $up,$down -ks -p 6 -cond UP,DOWN -out PROP_FREQ.$short | tee PROP_FREQ.$short.log ; done

#compute Density graphs for aa properties UP vs Down
for up in *up.prop_freq.txt ; do down=`echo $up | sed s/_up/_down/`; short=`echo $up | sed s/_up//`;  Panalyse_NT_AA_graphs.pl -mode freq -in $up,$down -ks -p 6 -cond UP,DOWN -out PROP_FREQ.$short | tee PROP_FREQ.$short.log ; done

#compute Density graphs  for aa properties TRA2 MDA vs SRSF2 K562
# $Rpalette_color='c("'.$PALETTE{'TRA2'}.'","'.$PALETTE{'SRSF2'}.'","#000000")';

title=TRA2_MDA_DOWN_vs_SRSF2_K562_DOWN ;  Panalyse_NT_AA_graphs.pl -mode freq -in 08_1_TRA2A-B_GSE59335_MDA-MB-231_query_results_down.prop_freq.txt,02_1_SRSF2_GSE65349_K562_query_results_down.prop_freq.txt  -ks -p 7 -cond TRA2,SRSF2 -out PROP_FREQ.$title | tee PROP_FREQ.$title.graph.log

title=TRA2_MDA_DOWN_vs_SRSF2_all_DOWN ;  Panalyse_NT_AA_graphs.pl -mode freq -in 08_1_TRA2A-B_GSE59335_MDA-MB-231_query_results_down.prop_freq.txt,2_Union_SRSF2_down.prop_freq.txt  -ks -p 7 -cond TRA2,SRSF2 -out PROP_FREQ.$title | tee PROP_FREQ.$title.graph.log


#compute Density graphs  for aa properties SRSF3 HepG2 vs SRSF1 HepG2
# $Rpalette_color='c("'.$PALETTE{'SRSF3'}.'","'.$PALETTE{'SRSF1'}.'","#000000")';

title=SRSF3_HepG2_DOWN_vs_SRSF1_HepG2_DOWN ;  Panalyse_NT_AA_graphs.pl -mode freq -in 03_2_SRSF3_ENCSR376FGR_HepG2_query_results_down.prop_freq.txt,01_2_SRSF1_ENCSR094KBY_HepG2_query_results_down.prop_freq.txt  -ks -p 8 -cond  SRSF3,SRSF1 -out PROP_FREQ.$title | tee PROP_FREQ.$title.graph.log

title=SRSF3_all_DOWN_vs_SRSF1_all_DOWN ;  Panalyse_NT_AA_graphs.pl -mode freq -in 3_Union_SRSF3_down.prop_freq.txt,1_Union_SRSF1_down.prop_freq.txt  -ks -p 8 -cond  SRSF3,SRSF1 -out PROP_FREQ.$title | tee PROP_FREQ.$title.graph.log

# Dispatch images
for dir in aa_features n ; do mkdir -p Images/$dir ; done
for dir in hnRNPH1 hnRNPK hnRNPL SRSF1 SRSF2 SRSF3 PTBP1 TRA2 SRSF3_vs_SRSF1 TRA2_vs_SRSF2 ; do  mkdir -p Images/aa_features/$dir ; done
for dir in hnRNPH1 hnRNPK hnRNPL SRSF1 SRSF2 SRSF3 PTBP1 TRA2 ; do  mkdir -p Images/n/$dir ; done

mkdir -p R TXT PDF
mv *pdf PDF/
mv *R R/
mv *figures.txt TXT/

mv *prop_freq*png Images/aa_features/
mv *n_freq*png Images/n/
mv PROP_FREQ.TRA2*png Images/aa_features/TRA2_vs_SRSF2/
mv PROP_FREQ.SRSF*png Images/aa_features/SRSF3_vs_SRSF1/
cd Images/aa_features/
for i in  hnRNPH1 hnRNPK hnRNPL SRSF1 SRSF2 SRSF3 PTBP1 TRA2 ; do mv *$i*png $i/ ; done
cd ../n/
for i in  hnRNPH1 hnRNPK hnRNPL SRSF1 SRSF2 SRSF3 PTBP1 TRA2 ; do mv *$i*png $i/  ; done

my %PROPERTIES;
#$PROPERTIES{'tiny'}="A,C,G,S,T";
$PROPERTIES{'very_small'}="A,C,G,S";
$PROPERTIES{'smallII'}="A,C,D,G,N,P,S,T";
#$PROPERTIES{'smallI'}="A,C,D,G,N,P,S,T,V";
$PROPERTIES{'large'}="F,I,K,L,M,R,W,Y";

$PROPERTIES{'disorder_promII'}="A,E,G,K,P,Q,S";
$PROPERTIES{'disorder_promI'}="A,E,G,K,P,Q,R,S";
$PROPERTIES{'order_promI'}="C,F,I,L,N,V,W,Y";
#$PROPERTIES{'order_promII'}="C,F,H,I,L,M,N,V,W,Y";

#$PROPERTIES{'polar_unchargedII'}="N,P,Q,S,T";
$PROPERTIES{'polar_unchargedII'}="N,Q,S,T,Y";
#$PROPERTIES{'polar_unchargedI'}="C,N,P,Q,S,T";
#$PROPERTIES{'polar_charged'}="D,E,H,K,R";
$PROPERTIES{'polar_chargedII'}="D,E,K,R";

#$PROPERTIES{'hydrophobicII'}="A,C,F,I,L,M,P,V,W,Y";
$PROPERTIES{'hydrophobicI'}="A,C,F,I,L,M,V";
$PROPERTIES{'hydrophilicI'}="D,E,K,N,Q,R";
#$PROPERTIES{'hydrophilicI'}="D,E,H,K,N,Q,R";
#$PROPERTIES{'hydrophilicII'}="D,E,H,K,N,Q,R,S,T";

$PROPERTIES{'neutral'}="G,H,P,S,T,Y";
#$PROPERTIES{'chargedI'}="D,E,H,K,R";
#$PROPERTIES{'chargedII'}="D,E,K,R";

$PROPERTIES{'positively_chargedII'}="K,R";
#$PROPERTIES{'positively_charged'}="H,K,R";
$PROPERTIES{'negatively_charged'}="D,E";
$PROPERTIES{'hydroxyl'}="S,T,Y";

my @PAIRED_PROP=("very_small::large",
                 "polar_unchargedII::polar_chargedII",
                 "neutral::polar_chargedII",
                 "disorder_promI::order_promI",
                 "hydrophilicI::hydrophobicI",
                 "negatively_charged::hydroxyl","negatively_charged::positively_chargedII");

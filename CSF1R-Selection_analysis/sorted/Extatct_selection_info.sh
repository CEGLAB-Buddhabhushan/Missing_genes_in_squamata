echo -e "Species_name\tFEL_positive_sites\tFEL_nagative_sites\tMEME_Episodic_div_site\trelax_pvalue\trelax_kvalue\tabsrel_pvalue\tbusted_pvalue" > FEL_MEME_RELAX_aBSREL_BUSTED.out.tsv
for i in `cat Species.lst`
do
fel_pos=`grep "pervasive positive diversifying" ../CSF1R/FEL/"$i"_treeoutput_fel | grep -oP '_\K\d+(?=_)'|head -1`
fel_neg=`grep "pervasive positive diversifying" ../CSF1R/FEL/"$i"_treeoutput_fel | grep -oP '_\K\d+(?=_)'|tail -1`
meme_EDS=`grep "episodic diversifying positive selection" ../CSF1R/MEME/"$i"_treeoutput_meme | sed -n 's/.*_\(.*\)_.*/\1/p'`
relax_pvalue=`grep -w "$i" ../CSF1R/RELAX/HYPHY_RELAX.Results.sorted.txt|cut -f2`
relax_kvalue=`grep -w "$i" ../CSF1R/RELAX/HYPHY_RELAX.Results.sorted.txt|cut -f3`
absrel_pvalue=`grep -w "$i" ../CSF1R/aBSREL/HYPHY_aBSREL.Results.txt|cut -f2`
busted_pvalue=`grep -w "$i" ../CSF1R/BUSTED/HYPHY_BUSTED.Results.txt|cut -f2`
echo -e "$i\t$fel_pos\t$fel_neg\t$meme_EDS\t$relax_pvalue\t$relax_kvalue\t$absrel_pvalue\t$busted_pvalue" >> FEL_MEME_RELAX_aBSREL_BUSTED.out.tsv
echo $i
done


awk 'BEGIN{OFS="\t"} NR==1{print $1,$2,$3,$4,$5,$6,"Selection"$7,$8; next} {print $1,$2,$3,$4,$5,$6,"Episodic diversifying selection",$7,$8}' HYPHY_MEME_merged_data.tsv > HYPHY_MEME_merged_data_with_selection.tsv
cat HYPHY_FEL_merged_data.tsv <(tail -n+2 HYPHY_MEME_merged_data_with_selection.tsv)> HYPHY_FEL-MEME_merged_data.tsv


#relax
for i in SAMD9 VWA7 SAMD9L HSD17B13 19spCUZD1 CUZD1
do
cd $i/RELAX
cut -f1,3,4 HYPHY_RELAX.Results.txt > HYPHY_RELAX.Results.sorted.txt
Rscript ../../adj.p-value_relax.R
Rscript ../../adj.p-value_relax_logp.R
cut -f2- HYPHY_RELAX.Results.p_val_with_log10.tsv > "$i".Final_HYPHY_RELAX.Results.p_val_with_log10.tsv
cp "$i".Final_HYPHY_RELAX.Results.p_val_with_log10.tsv /media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/Thylacine_selection/RELAX
cd ../..
echo $i
done

##do this for  aBSREL BUSTED FEL MEME
for i in SAMD9 VWA7 SAMD9L HSD17B13 19spCUZD1 CUZD1
do
cd $i/MEME
cp HYPHY_MEME_merged_data.tsv /media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/Thylacine_selection/MEME/"$i".HYPHY_MEME_merged_data.tsv
cp HYPHY_MEME_merged_data.tsv.png /media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/Thylacine_selection/MEME/"$i".HYPHY_MEME_merged_data.tsv.png
cd ../..
echo $i
done


##codeml
for i in SAMD9 VWA7 SAMD9L HSD17B13 19spCUZD1 CUZD1
do
cd $i/codeml
echo  "clade mdl fgspecies bgspecies om0 lnlm0 npm0 ombnfg ombnbg lnlbn npbn ombffg ombfbg lnlbf npbf pvalm0_bfree pavalbne_bfree" |sed 's/ /\t/g' > Final_PAML-codeml.Results.txt
cat Final_result.txt >> Final_PAML-codeml.Results.txt
Rscript ../../adj.p-value_codeml.R
cut -f2- Final_PAML-codeml.Results.p_adj.tsv |sed 's/"//g' > "$i".Final_PAML-codeml.Results.p_adj.tsv
cp "$i".Final_PAML-codeml.Results.p_adj.tsv /media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/Thylacine_selection/codeml
cd ../..
echo $i
done

##gBGC
for i in SAMD9 VWA7 SAMD9L HSD17B13 19spCUZD1 CUZD1
do
cd $i/gBGC
bash ../../gBGC_higher_than_0.5.sh "$i"_phastout > "$i".Species_with_gBGC_higher_than_0.5.list
cp "$i".Species_with_gBGC_higher_than_0.5.list /media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/Thylacine_selection/gBGC
cp "$i".jpeg /media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/Thylacine_selection/gBGC
cd ../..
echo $i
done

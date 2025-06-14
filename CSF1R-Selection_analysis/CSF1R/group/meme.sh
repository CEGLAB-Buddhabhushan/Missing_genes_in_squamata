##generate MEME output using HYPHY
clade=$1
fg=$2
cd $clade
mkdir MEME
cp "$clade".aln "$clade".nwk MEME
cd MEME
tree="$clade".nwk
for j in `cat ../$fg`
do
sed -i "s/$j/$j{fg}/g" $tree 
done

HYPHYMPI meme --alignment "$clade".aln --tree $tree --branches fg  > treeoutput_meme


##makking table
mkdir MEME_output
echo -e "Codon\tPartition\talpha\tbeta_plus\tp_plus\tLRT\tpValue\tSpecies" > HYPHY_MEME_merged_data.tsv
for i in treeoutput_meme
do
species=`echo $i|sed 's/_treeoutput_meme//g'`
echo $species
awk '/### For partition/,0' $i|awk 'BEGIN { FS="|"; OFS="\t" } NR>2'|grep -v "###"|sed 's/:----------://g;s/|/\t/g;s/ //g;s/:--------------------------://g;s/Yes,p=//g;s/\t//1;s/Episodicselectiondetected?/pValue/g;s/beta+/beta_plus/g;s/p+/p_plus/g'|grep "\S" |cut -f1-7 > "$species".tsv

awk -v species="$species" 'BEGIN{OFS="\t"; print "Codon\tPartition\talpha\tbeta_plus\tp_plus\tLRT\tpValue\tSpecies"} NR>1{$8=species; print}' "$species".tsv > "$species"_with_species.tsv

cat "$species"_with_species.tsv|awk 'NR>2' >> HYPHY_MEME_merged_data.tsv
rm "$species".tsv
mv "$species"_with_species.tsv MEME_output 
done
Rscript ../../MEME_plot.R HYPHY_MEME_merged_data.tsv
cd ..
mv MEME ${fg%.lst}_MEME
cd ..

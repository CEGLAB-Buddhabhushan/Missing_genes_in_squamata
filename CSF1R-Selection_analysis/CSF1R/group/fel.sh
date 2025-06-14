##generate FEL output using HYPHY
clade=$1
fg=$2
cd $clade
mkdir FEL
cp "$clade".aln "$clade".nwk FEL
cd FEL
tree="$clade".nwk
for j in `cat ../$fg`
do
sed -i "s/$j/$j{fg}/g" $tree 
done

HYPHYMPI fel --alignment "$clade".aln --tree $tree --branches fg > treeoutput_fel



##makking table
mkdir FEL_output
echo -e "Codon\tPartition\talpha\tbeta\tLRT\tSelection\tpValue\tSpecies" > HYPHY_FEL_merged_data.tsv
for i in treeoutput_fel
do
species=`echo $i|sed 's/_treeoutput_fel//g'`
echo $species
awk '/### For partition/,0' $i|awk 'BEGIN { FS="|"; OFS="\t" } NR>2'|grep -v "###"|sed 's/:--------------://g;s/|/\t/g;s/ //g;s/:-----------------://g;s/p = //g;s/\t//1;s/Selectiondetected?/Selection/g;s/Neg./Negative/g;s/Pos./Positive/g;s/p=/\t/g'|grep "\S" > "$species".tsv

awk -v species="$species" 'BEGIN{OFS="\t"; print "Codon\tPartition\talpha\tbeta\tLRT\tSelection\tpValue\tSpecies"} NR>1{$8=species; print}' "$species".tsv > "$species"_with_species.tsv

cat "$species"_with_species.tsv|awk 'NR>1'|head -n -1|sed 's/.*[1-9]m//g' >> HYPHY_FEL_merged_data.tsv
rm "$species".tsv
mv "$species"_with_species.tsv FEL_output 
done
Rscript ../../FEL_plot.R HYPHY_FEL_merged_data.tsv
cd ..
mv FEL ${fg%.lst}_FEL
cd ..


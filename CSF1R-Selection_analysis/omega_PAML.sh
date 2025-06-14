clade=$1
cd $clade
mkdir omega_PAML
cp "$clade".aln "$clade".nwk omega_PAML
cp ../scripts/codeml1 omega_PAML
cd omega_PAML

##calculating the dn and ds values by codeml
j=$1
tree="$j".nwk

for i in "$j".aln
do
j=`echo $i|cut -d. -f1`
sed  "s/aln/$i/g" codeml1 > codeml.ctl
sed  -i "s/bestTree/$tree/g" codeml.ctl
sed -i "s/output_file/$i.out/g" codeml.ctl
sed -i "s/model = 1/model = 0/g" codeml.ctl
codeml codeml.ctl
grep -A1 "dN tree" $i.out > dn_paml_tree
grep -A1 "dS tree" $i.out > ds_paml_tree

sed '1d' dn_paml_tree | tr "," "\n" |sed 's/(//g' |sed 's/).*//g' |sed 's/:/\t/g' |awk '{print $1,$2}' |sed 's/ /\t/g' > dn_paml_tree.1
sed '1d' ds_paml_tree | tr "," "\n" |sed 's/(//g' |sed 's/).*//g' |sed 's/:/\t/g' |awk '{print $1,$2}' |sed 's/ /\t/g' > ds_paml_tree.1

paste dn_paml_tree.1 ds_paml_tree.1 > dnds_paml_tree.result
sort -k1,1 dnds_paml_tree.result > dnds_paml_tree.result.final
done
echo "Species_name dN dS omega"|sed 's/ /\t/g' > Final_dnds_paml_tree.tsv
awk '{print $1, $2, $4, $2/$4}' dnds_paml_tree.result.final|sed 's/ /\t/g' >> Final_dnds_paml_tree.tsv
Rscript /media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/plot_tree_with_omega.R "$clade".nwk Final_dnds_paml_tree.tsv
cd ..


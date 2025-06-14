##generate aBSREL output using HYPHY
clade=$1
fg=$2
cd $clade
mkdir aBSREL
cp "$clade".aln "$clade".nwk aBSREL
cd aBSREL
tree="$clade".nwk
for j in `cat ../$fg`
do
sed -i "s/$j/$j{fg}/g" $tree
done
HYPHYMPI aBSREL --alignment "$clade".aln --tree $tree --branches fg > treeoutput_aBSREL

##makking table

echo -e "test pval"|sed 's/ /\t/g' > HYPHY_aBSREL.Results.txt
for d in treeoutput_aBSREL
do
pval=$(grep ", p-value =  " "$d" | cut -f2 -d'=' | sed 's/ //g')
test=$(grep "Selected 1 branches for testing:" "$d" | cut -f2 -d':' | sed 's/ //g; s/`//g')
if [ -z "$pval" ]
then pval="NA"
fi

echo -e "$test\t$pval" >> HYPHY_aBSREL.Results.txt
done

cd ..
mv aBSREL ${fg%.lst}_aBSREL
cd ..

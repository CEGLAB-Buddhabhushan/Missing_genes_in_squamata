##generate BUSTED output using HYPHY
clade=$1
fg=$2
cd $clade
mkdir BUSTED
cp "$clade".aln "$clade".nwk BUSTED
cd BUSTED
tree="$clade".nwk

for j in `cat ../$fg`
do
sed -i "s/$j/$j{fg}/g" $tree
done
HYPHYMPI busted --alignment "$clade".aln --tree $tree --branches fg > treeoutput_BUSTED


##makking table

echo -e "test pval"|sed 's/ /\t/g' > HYPHY_BUSTED.Results.txt
for d in treeoutput_BUSTED
do
pval=$(grep "Likelihood ratio test for episodic diversifying positive selection" "$d" | cut -f2 -d'=' | sed 's/ //g;s/[*]//g;s/\.//2')
test=$(grep "Selected 1 branches to test in the BUSTED analysis:" "$d" | cut -f2 -d':' | sed 's/ //g; s/`//g')

if [ -z "$pval" ]
then pval="NA"
fi

echo -e "$test\t$pval" >> HYPHY_BUSTED.Results.txt
done

cd ..
mv BUSTED ${fg%.lst}_BUSTED
cd ..

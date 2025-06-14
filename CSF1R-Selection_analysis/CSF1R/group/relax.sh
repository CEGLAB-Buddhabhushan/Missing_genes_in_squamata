clade=$1
cd $clade
fg=$2
bg=$3
mkdir RELAX
cp "$clade".aln "$clade".nwk RELAX
cd RELAX
tree="$clade".nwk

for j in `cat ../$fg`
do
sed -i "s/$j/$j{fg}/g" $tree
done
for j in `cat ../$bg`
do
sed -i "s/$j/$j{bg}/g" $tree 
done

HYPHYMPI relax --alignment "$clade".aln --tree $tree --test fg --reference bg > treeoutput_relax


### taking HYPHYMP output result in one file
echo -e "test back pval kval"|sed 's/ /\t/g' > HYPHY_RELAX.Results.txt
for d in `ls -1 treeoutput_relax`
do
pval=`grep "^Like" $d|awk '{print $6}'|sed 's/\*\*\.//g'|head -1`
kval1=`grep "Relaxation/intensification" $d|awk '{print $6}'|head -1`
test=`grep "_Test_ set:" $d|awk '{print $9}'|sed 's/\`//g'`
back=`grep "_Reference_ set:" $d|cut -f2 -d ":"|sed 's/ //g'`

echo -e "$test $back $pval $kval1"|sed 's/ /\t/g' >> HYPHY_RELAX.Results.txt
done
cut -f1,3,4 HYPHY_RELAX.Results.txt > HYPHY_RELAX.Results.sorted.txt
cd ..
mv RELAX ${fg%.lst}-${bg%.lst}_RELAX
cd ..

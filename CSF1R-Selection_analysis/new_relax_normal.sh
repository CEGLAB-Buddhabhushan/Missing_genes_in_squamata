#!/bin/bash
#conda activate hyphy_env
clade=$1
cd $clade
mkdir new_RELAX_normal
cp "$clade".aln "$clade".nwk new_RELAX_normal
cd new_RELAX_normal
tree="$clade".nwk
for i in "$clade".aln
do
t=`grep ">" $i|wc -l`
grep ">" $i|sed 's/>//g' > taxlist.txt
for j in `cat taxlist.txt`
do
sed "s/$j/$j{fg}/g" $tree > "$j"_treeLabled.txt
hyphy relax --alignment $i --tree "$j"_treeLabled.txt --test fg --output "$j"_treeoutput_relax
done
done


### taking HYPHYMP output result in one file
echo -e "Test\tLRT\tp-value\tK" > HYPHY_RELAX.Results.txt
for d in `ls -1 *_relax`
do
test=`grep "Test" $d |tail -1|cut -f1 -d':'|sed 's/"//g;s/ //g'`
lrt=`grep "LRT" $d |cut -f2 -d':'|sed 's/,//g'`
pval=`grep "p-value" $d |cut -f2 -d':'|sed 's/,//g'`
kval=`grep "relaxation or intensification parameter" $d |cut -f2 -d':'|sed 's/,//g'`
echo -e "$test\t$lrt\t$pval\t$kval" >> HYPHY_RELAX.Results.txt
done
cd ../..

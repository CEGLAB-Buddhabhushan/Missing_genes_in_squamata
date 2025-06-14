for gene in CSF1R
do
for fg in Alligator_mississippiensis.lst Amphibia.lst Chelonioidea.lst Clupeocephala.lst Episquamata.lst Latimeria_chalumnae.lst Mammalia.lst Neognathae.lst Sphenodon_punctatus.lst
do
./aBSREL.sh $gene $fg
./BUSTED.sh $gene $fg
./meme.sh $gene $fg
./fel.sh $gene $fg
done
done




for gene in CSF1R
do
species_groups=(Alligator_mississippiensis.lst Amphibia.lst Chelonioidea.lst Clupeocephala.lst Episquamata.lst Latimeria_chalumnae.lst Mammalia.lst Neognathae.lst Sphenodon_punctatus.lst)
for ((i = 0; i < ${#species_groups[@]}; i++))
do
fg=${species_groups[i]}
for ((j = i + 1; j < ${#species_groups[@]}; j++))
do
bg=${species_groups[j]}
./relax.sh $gene $fg $bg
./relax.sh $gene $bg $fg
done
done
done

for gene in CSF1R
do
species_groups=(Mammalia.lst Sauria.lst)
for ((i = 0; i < ${#species_groups[@]}; i++))
do
fg=${species_groups[i]}
for ((j = i + 1; j < ${#species_groups[@]}; j++))
do
bg=${species_groups[j]}
./relax.sh $gene $fg $bg
./relax.sh $gene $bg $fg
done
done
done


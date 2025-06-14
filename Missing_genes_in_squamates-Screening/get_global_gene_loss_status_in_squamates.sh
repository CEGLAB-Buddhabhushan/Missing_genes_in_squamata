echo -e "Gene\tPantherophis_guttatus\tPituophis_catenifer_annectens\tThamnophis_elegans\tThamnophis_sirtalis\tAhaetulla_prasina\tCrotalus_tigris\tProtobothrops_mucrosquamatus\tNotechis_scutatus\tPseudonaja_textilis\tErythrolamprus_reginae\tPython_bivittatus\tCandoia_aspera\tAnolis_sagrei\tAnolis_carolinensis\tSceloporus_undulatus\tPogona_vitticeps\tVaranus_komodoensis\tElgaria_multicarinata_webbii\tPodarcis_muralis\tPodarcis_raffonei\tZootoca_vivipara\tLacerta_agilis\tRhineura_floridana\tTiliqua_scincoides\tHemicordylus_capensis\tSphaerodactylus_townsendi\tEuleptes_europaea\tHeteronotia_binoei\tGekko_japonicus\tEublepharis_macularius\tAlligator_mississippiensis\tGallus_gallus\tChelonia_mydas\tHomo_sapiens\tXenopus_tropicalis" |sed 's/\t/,/g' > Gene_events_squamates.csv

for gene_id in $(cat gene_id-name.lst)
do
i=`echo $gene_id|cut -f1 -d':'`
gene_name=`echo $gene_id|cut -f2 -d':'`
cd $i/ncbi_dataset/data/ || exit
echo -e "Species_name\t$i" > "$i".gene_status.tsv
/media/disk1/BUDDHA/BMC-BIO/dataformat tsv gene-product --inputfile product_report.jsonl > squamates_product_report.dataformat.tsv
head -1 squamates_product_report.dataformat.tsv|cut -f1,3,7-11,16,17,18,33 > squamates_LOW_QUALITY_PROTEIN.tsv
grep "LOW QUALITY PROTEIN" squamates_product_report.dataformat.tsv |cut -f1,3,7-11,16,17,18,33|sort -u >> squamates_LOW_QUALITY_PROTEIN.tsv

for acc in `cut -f6 squamates_LOW_QUALITY_PROTEIN.tsv|tail -n+2`
do
echo $acc
rm "$acc".gb
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acc}&rettype=gb&retmode=txt&api_key=61131f45442b718e60396eb45ce73b65bc08">> "$acc".gb
organism=`grep "/organism" "$acc".gb|sed 's/organism="//g;s/[/]//g;s/"//g'`
Note=`awk '/\/note="/ {flag=1; print; next} flag {if (/"/) {flag=0; print} else {print}}' "$acc".gb | sed 's/"//g;s/note=//g;s/[/]//g' | tr '\n' ' ' | sed 's/  */ /g;s/^ //;s/ $//'`
#check internal Stop
if grep -A3 "##RefSeq-Attributes-START##" "$acc".gb | grep -v "##" | sed 's/^ *//' | grep -q "internal stop codons"; then
    internal_stop="Yes"
else
    internal_stop="No"
fi
#Check for "frameshifts"
if grep -A3 "##RefSeq-Attributes-START##" "$acc".gb | grep -v "##" | sed 's/^ *//' | grep -q "frameshifts"; then
    frameshifts="Yes"
else
    frameshifts="No"
fi
#Check "support for all annotated introns"
if echo $Note|grep -q "support for all annotated introns"; then
    annotated_introns="Yes"
else
    annotated_introns="No"
fi

echo -e "$acc\t$organism\t$internal_stop\t$frameshifts\t$annotated_introns\t$coverage\t$GC_content\t$Note" >> squamates_LOW_QUALITY_PROTEIN.efetch.tsv
done

for sp in $(cat /media/disk1/BUDDHA/BMC-BIO/Bifurcata_and_outgroup_species.lst)
do
species_name=$(echo "$sp" | sed 's/_/ /g')
# Check if the species name exists in the file
if grep -q -w "$species_name" squamates_product_report.dataformat.tsv; then
# Check if it is associated with "LOW QUALITY PROTEIN"
if grep -w "$species_name" squamates_product_report.dataformat.tsv | grep -q "LOW QUALITY PROTEIN"; then
# Check if all annotated introns
if grep -w "$species_name" squamates_product_report.dataformat.tsv | cut -f5 | grep -q "Yes"; then
gene_status="LQP-Intact"
else
gene_status="LQP"
fi
else
gene_status="Intact"
fi
else
gene_status="Missing"
fi
echo -e "$sp\t$gene_status" >> "$i".gene_status.tsv
done
awk '{for (i=1; i<=NF; i++) a[i] = (a[i] ? a[i] FS : "") $i} END {for (i=1; i<=NF; i++) print a[i]}' "$i".gene_status.tsv |tr ' ' ',' > "$i".gene_status.transpose.csv
tail -1 "$i".gene_status.transpose.csv |sed "s/$i/$gene_name/g" >> /media/disk1/BUDDHA/Global_gene_loss/Human-gene_set/Gene_events_squamates.csv
cd /media/disk1/BUDDHA/Global_gene_loss/Human-gene_set/ || exit
echo $i
done

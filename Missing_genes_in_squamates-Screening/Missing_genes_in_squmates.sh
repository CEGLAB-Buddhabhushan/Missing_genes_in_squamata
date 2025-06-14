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


## GET THE GENES WHICH MISSING IN SQUAMATES
awk -F',' '                                        {                                                
    missing = 1;  # Assume all values in columns 2-31 are "Missing" initially
    for (i=2; i<=31; i++) {
        if ($i != "Missing" && $i != "LQP") {  # If any column is NOT "Missing", set flag to 0
            missing = 0;
            break;
        }
    }
    if (missing && $32=="Intact" && $33=="Intact" && $34=="Intact" && $35=="Intact" && $36=="Intact") {
        print $1;  # Print gene name (column 1)
    }
}' Gene_events_squamates.csv > missing_in_squamates.lst

for i in `cat missing_in_squamates.lst`; do grep -w "$i" gene_id-name.lst ; done > Complete_loss.lst
for i in $(cat missing_in_squamates.lst); do   awk -F',' -v gene="$i" '$1 == gene {print $0}' Gene_events_squamates.csv; done >> Missing_genes-Gene_events_squamates.csv


#### get cds for blast
sp1="Xenopus tropicalis"
sp2="Anolis carolinensis"
sp1sp2=`echo "${sp1}-${sp2}.gene_orthologs.tsv"|sed 's/ /_/g'`
echo -e "Gene\t${sp1}_LONGEST_acc\t${sp1}_LONGEST_length\tGene\t${sp2}_LONGEST_acc\t${sp2}_LONGEST_length" > $sp1sp2

# Ensure output FASTA files are empty before starting
> "${sp1}.cds.fasta"
> "${sp2}.cds.fasta"

while IFS=':' read -r id gene; do
    transcript_file="$id/ncbi_dataset/data/transcript_protein.tsv"
    
    if [[ ! -f "$transcript_file" ]]; then
        echo "Warning: File $transcript_file not found for gene $gene. Skipping..."
        continue
    fi
    
    # Extract longest transcript for species 1
    sp1_LONGEST_acc=$(grep "$sp1" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f4)
    sp1_LONGEST_length=$(grep "$sp1" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f5)
    
    if [[ -n "$sp1_LONGEST_acc" ]]; then
        efetch -db nuccore -id "$sp1_LONGEST_acc" -format fasta_cds_na 2>/dev/null | \
            sed "s/^>\(.*\) \[gene=.*$/>\1 $gene/" | \
            sed 's/^>\S*_cds_\S*_\S*_\S* />/' >> "${sp1}.cds.fasta"
    else
        echo "Warning: No longest transcript found for $gene in $sp1."
    fi
    
    # Extract longest transcript for species 2
    sp2_LONGEST_acc=$(grep "$sp2" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f4)
    sp2_LONGEST_length=$(grep "$sp2" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f5)
    
    if [[ -n "$sp2_LONGEST_acc" ]]; then
        efetch -db nuccore -id "$sp2_LONGEST_acc" -format fasta_cds_na 2>/dev/null | \
            sed "s/^>\(.*\) \[gene=.*$/>\1 $gene/" | \
            sed 's/^>\S*_cds_\S*_\S*_\S* />/' >> "${sp2}.cds.fasta"
    else
        echo "Warning: No longest transcript found for $gene in $sp2."
    fi
    
    echo -e "$gene\t$sp1_LONGEST_acc\t$sp1_LONGEST_length\t$gene\t$sp2_LONGEST_acc\t$sp2_LONGEST_length" >> $sp1sp2
done < gene_id-name.lst


## get syntenic genes
# Create header line
echo -e "Gene\tHomo_sapiens\tXenopus_tropicalis\tChelonia_mydas\tAlligator_mississippiensis\tGallus_gallus\tPantherophis_guttatus\tPituophis_catenifer_annectens\tThamnophis_elegans\tThamnophis_sirtalis\tAhaetulla_prasina\tCrotalus_tigris\tProtobothrops_mucrosquamatus\tNotechis_scutatus\tPseudonaja_textilis\tErythrolamprus_reginae\tPython_bivittatus\tCandoia_aspera\tAnolis_sagrei\tAnolis_carolinensis\tSceloporus_undulatus\tPogona_vitticeps\tVaranus_komodoensis\tElgaria_multicarinata_webbii\tPodarcis_muralis\tPodarcis_raffonei\tZootoca_vivipara\tLacerta_agilis\tRhineura_floridana\tTiliqua_scincoides\tHemicordylus_capensis\tSphaerodactylus_townsendi\tEuleptes_europaea\tHeteronotia_binoei\tGekko_japonicus\tEublepharis_macularius" > Squamates_all_genes_synteny.tsv

# Store species list in array
species_list=(Homo_sapiens Xenopus_tropicalis Chelonia_mydas Alligator_mississippiensis Gallus_gallus \
Pantherophis_guttatus Pituophis_catenifer_annectens Thamnophis_elegans Thamnophis_sirtalis \
Ahaetulla_prasina Crotalus_tigris Protobothrops_mucrosquamatus Notechis_scutatus Pseudonaja_textilis \
Erythrolamprus_reginae Python_bivittatus Candoia_aspera Anolis_sagrei Anolis_carolinensis \
Sceloporus_undulatus Pogona_vitticeps Varanus_komodoensis Elgaria_multicarinata_webbii \
Podarcis_muralis Podarcis_raffonei Zootoca_vivipara Lacerta_agilis Rhineura_floridana \
Tiliqua_scincoides Hemicordylus_capensis Sphaerodactylus_townsendi Euleptes_europaea \
Heteronotia_binoei Gekko_japonicus Eublepharis_macularius)

# Loop through each gene
while IFS=':' read -r id gene; do
    row="$gene"

    for species in "${species_list[@]}"; do
        sp=$(echo "$species" | sed 's/_/ /g')
        transcript_file="$id/ncbi_dataset/data/transcript_protein.tsv"
        gene_product_file="$id/ncbi_dataset/data/gene_product.tsv"

        if [[ ! -f "$transcript_file" || ! -f "$gene_product_file" ]]; then
            row+="\tNA"
            continue
        fi

        # Get longest transcript accession
        acc=$(grep -F "$sp" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f4)

        if [[ -z "$acc" ]]; then
            row+="\tNA"
            continue
        fi

        # Get info from gene_product file
        info=$(grep -F "$acc" "$gene_product_file" | cut -f7,18,24-26 | sort -u | \
            awk -F'\t' '{print $1" "$2" "$4" "$5" "$3}' | sed 's/  */,/g;s/minus/-/g;s/plus/+/g')

        if [[ -z "$info" ]]; then
            row+="\tNA"
        else
            row+="\t$info"
        fi
    done

    echo -e "$row" >> Squamates_all_genes_synteny.tsv
done < gene_id-name.lst

#### get the CDS of five outgroup
echo -e "Gene\tHomo_sapiens\tXenopus_tropicalis\tChelonia_mydas\tAlligator_mississippiensis\tGallus_gallus" > Outgroup_acc.tsv

# Store species list in array
species_list=(Homo_sapiens Xenopus_tropicalis Chelonia_mydas Alligator_mississippiensis Gallus_gallus)

# Loop through each gene
while IFS=':' read -r id gene; do
    row="$gene"

    for species in "${species_list[@]}"; do
        sp=$(echo "$species" | sed 's/_/ /g')
        transcript_file="$id/ncbi_dataset/data/transcript_protein.tsv"
        gene_product_file="$id/ncbi_dataset/data/gene_product.tsv"

        if [[ ! -f "$transcript_file" || ! -f "$gene_product_file" ]]; then
            row+="\tNA"
            continue
        fi

        # Get longest transcript accession
        acc=$(grep -F "$sp" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f4)

        if [[ -z "$acc" ]]; then
            row+="\tNA"
            continue
        fi

        if [[ -z "$acc" ]]; then
            row+="\tNA"
        else
            row+="\t$acc"
        fi
    done

    echo -e "$row" >> Outgroup_acc.tsv
done < Complete_loss.lst
# efetch cds

for i in `cut -f2 -d':' Complete_loss.lst`
do
for accession in `grep -w "$i" Outgroup_acc.tsv|cut -f2-|tr '\t' '\n'`
do
species_name=`efetch -db nuccore -id "$accession" -format summary|grep "ORGANISM"|awk -F' ' '{print $2"_"$3}'`
efetch -db nuccore -id "$accession" -format fasta_cds_na | sed "s/^>\(.*\) \[gene=.*$/>\1 $species_name:$accession-$i/"|sed 's/^>\S*_cds_\S*_\S*_\S* />/' >> "$i".fa
sed -i 's/:/_/g' "$i".fa
done
done


########################### Hybpiper
#conda activate hybpiper
for i in `cat cds_link.lst`
do
sp=`echo $i|cut -f1 -d'-'`
cds_link=`echo $i|cut -f2 -d'-'`
wget $cds_link
zcat "${cds_link##*/}" | sed 's/^.*gene=/>/' | cut -f1,4 -d']'|sed 's/\[protein_id=//g;s/\]//g;s/ /-/g' |sed "s/-/-$sp-/g" > "$sp".cds_from_genomic.fna
echo $sp
done
cat  *.cds_from_genomic.fna > 39_species.cds_from_genomic.fna
makeblastdb -in 39_species.cds_from_genomic.fna -out 39_species.cds_from_genomic.fna -dbtype nucl
grep ">" 39_species.cds_from_genomic.fna|sort -u |wc -l
#1867546

for i in *.fa
do
gene=$(echo $i | sed 's/\.fa//g')
rm -r $gene
mkdir $gene
cd $gene
cp ../$i .

### Naja naja Illumina HiSeq 2500
blastn -task blastn  -evalue 0.05 -db /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Naja_naja/DNA/SRR10428157.fa -out blastn.Naja_naja."$gene".1out -num_threads 64 -outfmt 1 -query $i -max_target_seqs 10000
#get read ids list
grep "^SRR" blastn.Naja_naja."$gene".1out|cut -f1 -d' '|sort -u > Naja_naja."$gene".SRR10428157.lst
#extract reads
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Naja_naja/DNA/SRR10428157_1.fastq.gz Naja_naja."$gene".SRR10428157.lst > Naja_naja."$gene".SRR10428157_1.fq
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Naja_naja/DNA/SRR10428157_2.fastq.gz Naja_naja."$gene".SRR10428157.lst > Naja_naja."$gene".SRR10428157_2.fq
#hybpiper
hybpiper check_targetfile --targetfile_dna $i
hybpiper assemble -r Naja_naja."$gene".SRR10428157_1.fq Naja_naja."$gene".SRR10428157_2.fq -t_dna $i --prefix Naja_naja.SRR10428157.blast --evalue 0.05 --cpu 64 --max_target_seqs 5000
#retrive sequence
hybpiper retrieve_sequences --targetfile_dna $i --single_sample_name Naja_naja.SRR10428157.blast --fasta_dir Naja_naja.SRR10428157.extracted_seq dna

### Pogona_vitticeps nanopore
blastn -task blastn  -evalue 0.05 -db /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Pogona_vitticeps_nanopore/SRR31361885.fa -out blastn.Pogona_vitticeps."$gene".1out -num_threads 64 -outfmt 1 -query $i -max_target_seqs 10000
#get read ids list
grep "^SRR" blastn.Pogona_vitticeps."$gene".1out|cut -f1 -d' '|sort -u > Pogona_vitticeps."$gene".SRR31361885.lst
#extract reads
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Pogona_vitticeps_nanopore/SRR31361885_1.fastq.gz Pogona_vitticeps."$gene".SRR31361885.lst > Pogona_vitticeps."$gene".SRR31361885.fq
#hybpiper
hybpiper check_targetfile --targetfile_dna $i
hybpiper assemble -r Pogona_vitticeps."$gene".SRR31361885.fq -t_dna $i --prefix Pogona_vitticeps.SRR31361885.blast --evalue 0.05 --cpu 64 --max_target_seqs 5000
#retrive sequence
hybpiper retrieve_sequences --targetfile_dna $i --single_sample_name Pogona_vitticeps.SRR31361885.blast --fasta_dir Pogona_vitticeps.SRR31361885.extracted_seq dna

### Pogona_vitticeps  Illumina HiSeq 2000
blastn -task blastn  -evalue 0.05 -db /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Pogona_vitticeps_nanopore/ERR409938.fa -out blastn.Pogona_vitticeps."$gene".ERR409938.1out -num_threads 64 -outfmt 1 -query $i -max_target_seqs 10000
#get read ids list
grep "^ERR" blastn.Pogona_vitticeps."$gene".ERR409938.1out|cut -f1 -d' '|sort -u > Pogona_vitticeps."$gene".ERR409938.lst
#extract reads
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Pogona_vitticeps_nanopore/ERR409938_1.fastq.gz Pogona_vitticeps."$gene".ERR409938.lst > Pogona_vitticeps."$gene".ERR409938_1.fq
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Pogona_vitticeps_nanopore/ERR409938_2.fastq.gz Pogona_vitticeps."$gene".ERR409938.lst > Pogona_vitticeps."$gene".ERR409938_2.fq
#hybpiper
hybpiper check_targetfile --targetfile_dna $i
hybpiper assemble -r Pogona_vitticeps."$gene".ERR409938_1.fq Pogona_vitticeps."$gene".ERR409938_2.fq -t_dna $i --prefix Pogona_vitticeps.ERR409938.blast --evalue 0.05 --cpu 64 --max_target_seqs 5000
#retrive sequence
hybpiper retrieve_sequences --targetfile_dna $i --single_sample_name Pogona_vitticeps.ERR409938.blast --fasta_dir Pogona_vitticeps.ERR409938.extracted_seq dna

### Rhabdophis_nuchalis PacBio-HiFi-seq-Revio
blastn -task blastn  -evalue 0.05 -db /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Rhabdophis_nuchalis/SRR28573917.fa -out blastn.Rhabdophis_nuchalis."$gene".1out -num_threads 64 -outfmt 1 -query $i -max_target_seqs 10000
#get read ids list
grep "^SRR" blastn.Rhabdophis_nuchalis."$gene".1out|cut -f1 -d' '|sort -u > Rhabdophis_nuchalis."$gene".SRR28573917.lst
#extract reads
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Rhabdophis_nuchalis/SRR28573917_subreads.fastq.gz Rhabdophis_nuchalis."$gene".SRR28573917.lst > Rhabdophis_nuchalis."$gene".SRR28573917.fq
#hybpiper
hybpiper check_targetfile --targetfile_dna $i
hybpiper assemble -r Rhabdophis_nuchalis."$gene".SRR28573917.fq -t_dna $i --prefix Rhabdophis_nuchalis.SRR28573917.blast --evalue 0.05 --cpu 64 --max_target_seqs 5000
#retrive sequence
hybpiper retrieve_sequences --targetfile_dna $i --single_sample_name Rhabdophis_nuchalis.SRR28573917.blast --fasta_dir Rhabdophis_nuchalis.SRR28573917.extracted_seq dna

### Zootoca_vivipara PacBio - HiFi-Revio
blastn -task blastn  -evalue 0.05 -db /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Zootoca_vivipara/ERR14708140.fa -out blastn.Zootoca_vivipara."$gene".1out -num_threads 64 -outfmt 1 -query $i -max_target_seqs 10000
#get read ids list
grep "^ERR" blastn.Zootoca_vivipara."$gene".1out|cut -f1 -d' '|sort -u > Zootoca_vivipara."$gene".ERR14708140.lst
#extract reads
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Zootoca_vivipara/ERR14708140.fastq.gz Zootoca_vivipara."$gene".ERR14708140.lst > Zootoca_vivipara."$gene".ERR14708140.fq
#hybpiper
hybpiper check_targetfile --targetfile_dna $i
hybpiper assemble -r Zootoca_vivipara."$gene".ERR14708140.fq -t_dna $i --prefix Zootoca_vivipara.ERR14708140.blast --evalue 0.05 --cpu 64 --max_target_seqs 5000
#retrive sequence
hybpiper retrieve_sequences --targetfile_dna $i --single_sample_name Zootoca_vivipara.ERR14708140.blast --fasta_dir Zootoca_vivipara.ERR14708140.extracted_seq dna

### Anolis_sagrei Hi-C sequencing
blastn -task blastn  -evalue 0.05 -db /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Anolis_sagrei/SRR30050337.fa -out blastn.Anolis_sagrei."$gene".1out -num_threads 64 -outfmt 1 -query $i -max_target_seqs 10000
#get read ids list
grep "^SRR" blastn.Anolis_sagrei."$gene".1out|cut -f1 -d' '|sort -u > Anolis_sagrei."$gene".SRR30050337.lst
#extract reads
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Anolis_sagrei/SRR30050337_1.fastq.gz Anolis_sagrei."$gene".SRR30050337.lst > Anolis_sagrei."$gene".SRR30050337_1.fq
seqtk subseq /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Anolis_sagrei/SRR30050337_2.fastq.gz Anolis_sagrei."$gene".SRR30050337.lst > Anolis_sagrei."$gene".SRR30050337_2.fq
#hybpiper
hybpiper check_targetfile --targetfile_dna $i
hybpiper assemble -r Anolis_sagrei."$gene".SRR30050337_1.fq Anolis_sagrei."$gene".SRR30050337_2.fq -t_dna $i --prefix Anolis_sagrei.SRR30050337.blast --evalue 0.05 --cpu 64 --max_target_seqs 5000
#retrive sequence
hybpiper retrieve_sequences --targetfile_dna $i --single_sample_name Anolis_sagrei.SRR30050337.blast --fasta_dir Anolis_sagrei.SRR30050337.extracted_seq dna

################# use hybpiper output to find orthologs #################
cat Naja_naja.SRR10428157.extracted_seq/*.fasta Pogona_vitticeps.SRR31361885.extracted_seq/*.fasta Pogona_vitticeps.ERR409938.extracted_seq/*.fasta Rhabdophis_nuchalis.SRR28573917.extracted_seq/*.fasta Zootoca_vivipara.ERR14708140.extracted_seq/*.fasta Anolis_sagrei.SRR30050337.extracted_seq/*.fasta |cut -f1 -d' ' > "$gene".hybpiper.fa

#Check if FASTA file is empty
if [ ! -s "${gene}.hybpiper.fa" ]; then
    echo "${gene} not assembled" >> "$gene".hybpiper.log
    exit 0
fi

blastn -task blastn  -evalue 0.001 -db /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/cds_of_representative_species/39_species.cds_from_genomic.aggregate.fna -out blastn.39_species.cds_from_genomic."$gene".1out -num_threads 64 -outfmt 1 -query "$gene".hybpiper.fa -max_target_seqs 10000
blastn -task blastn  -evalue 0.001 -db /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/cds_of_representative_species/39_species.cds_from_genomic.aggregate.fna -out blastn.39_species.cds_from_genomic."$gene".6out -num_threads 64 -outfmt '6 qseqid sseqid evalue qlen qcovs' -query "$gene".hybpiper.fa -max_target_seqs 10000
awk -F'\t' '$5>=50 {print $2}' blastn.39_species.cds_from_genomic."$gene".6out|sort -u > blastn.39_species.cds_from_genomic."$gene".lst
seqtk subseq  /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/cds_of_representative_species/39_species.cds_from_genomic.aggregate.fna blastn.39_species.cds_from_genomic."$gene".lst > blastn.39_species.cds_from_genomic."$gene".fa

cat "$gene".hybpiper.fa blastn.39_species.cds_from_genomic."$gene".fa "$gene".fa > "$gene".combined.fa
awk '/^>/ {print; next} {gsub(/-/, ""); gsub(/X/, "N"); print}' "$gene.combined.fa" |seqtk seq |seqkit rmdup -s -o "$gene".combined.cleaned.fa

mafft --auto --anysymbol --quiet "$gene".combined.cleaned.fa > "$gene".combined.mafft.aln
iqtree2 -T AUTO -s "$gene".combined.mafft.aln -m MFP --alrt 1000 -B 1000 --boot-trees

bash /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/HybPiper/get_tree_bootstrap_info.sh $gene | sed 's/\\n/\n/g' > "$gene".hybpiper.iqtree.out
sister_branch=$(grep "is tip" "$gene".hybpiper.iqtree.out|cut -f2 -d':'|cut -f1 -d'-' |sed 's/ //g'|sort -u)
if grep -q "No immediate sister branch to MRCA found." "$gene".hybpiper.iqtree.out; then
    BS_for_cluster="NA"
    BS_for_immediate_parent_clade="NA"
else
BS_for_cluster=$(grep "Bootstrap support for" "$gene".hybpiper.iqtree.out | cut -f2 -d':' | head -1)
BS_for_immediate_parent_clade=$(grep "Bootstrap support for" "$gene".hybpiper.iqtree.out | cut -f2 -d':' | tail -1)
fi

# Output the results
echo -e "$gene,$sister_branch,$BS_for_cluster,$BS_for_immediate_parent_clade"|sed 's/,/\t/g' > "$gene".hybpiper.iqtree.tsv
cd ..
done

#collect output 
for i in `ls -d */`; do   gene=`echo $i | sed 's/[/]//g'`;   cat $gene/"$gene".hybpiper.iqtree.tsv; done | sed 's/-e //g' | cut -f1-2 | awk '
{ a[NR]=$1; b[NR]=$2 }
END {
  for (i=1; i<=NR; i++) printf "%s%s", a[i], (i==NR ? "\n" : "\t")
  for (i=1; i<=NR; i++) printf "%s%s", b[i], (i==NR ? "\n" : "\t")
}'> Hybpiper.sorted.out.tsv

##################################################################################

echo -e "Gene\tHomo_sapiens\tXenopus_tropicalis\tChelonia_mydas\tAlligator_mississippiensis\tGallus_gallus" > All_gene.Outgroup_acc.tsv

# Store species list in array
species_list=(Homo_sapiens Xenopus_tropicalis Chelonia_mydas Alligator_mississippiensis Gallus_gallus)

# Loop through each gene
while IFS=':' read -r id gene; do
    row="$gene"

    for species in "${species_list[@]}"; do
        sp=$(echo "$species" | sed 's/_/ /g')
        transcript_file="$id/ncbi_dataset/data/transcript_protein.tsv"
        gene_product_file="$id/ncbi_dataset/data/gene_product.tsv"

        if [[ ! -f "$transcript_file" || ! -f "$gene_product_file" ]]; then
            row+="\tNA"
            continue
        fi

        # Get longest transcript accession
        acc=$(grep -F "$sp" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f4)

        if [[ -z "$acc" ]]; then
            row+="\tNA"
            continue
        fi

        if [[ -z "$acc" ]]; then
            row+="\tNA"
        else
            row+="\t$acc"
        fi
    done

    echo -e "$row" >> All_gene.Outgroup_acc.tsv
done < gene_id-name.lst


#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_feature_table.txt.gz

cut -f1,5,7-9,15-16 GCF_000001405.40_GRCh38.p14_feature_table.txt|awk '$1=="gene"'|grep -A10 -B10 "IL34"|cut -f6 > IL34_human.lst
for i in `cat IL34_human.lst`
do
grep -w "$i" All_gene.Outgroup_acc.tsv
done | cut -f1 |grep -A2 -B2 "IL34"|grep -v "IL34" > IL34_human_sorted.lst

for i in `cat IL34_human_sorted.lst`
do
awk -F'\t' -v var="$i" '$1=="gene" && $15==var {print $15"-"$16}' GCF_000001405.40_GRCh38.p14_feature_table.txt >> IL34_syntenic_genes.lst
done 



conda activate ncbi_datasets
rm -r IL34_deletion
mkdir  IL34_deletion
cd IL34_deletion
mkdir CDS_filtered
for i in `cat ../IL34_syntenic_genes.lst`
do
gene=`echo $i | cut -f1 -d'-'`
gene_id=`echo $i | cut -f2 -d'-'`

datasets download gene gene-id $gene_id --include product-report --ortholog all --no-progressbar --api-key c181c44f2dca87803a9c0f9a24a32f67a608 --filename "$gene".zip
rm -rf "$gene"
unzip -q -o "$gene".zip -d "$gene"
rm "$gene".zip
cd "$gene/ncbi_dataset/data"
dataformat tsv gene-product --inputfile product_report.jsonl  > gene_product.tsv
dataformat tsv gene-product --inputfile product_report.jsonl --fields gene-id,tax-name,symbol,transcript-accession,transcript-length,transcript-protein-accession,transcript-protein-length > transcript_protein.tsv

mkdir -p CDS_sorted

for species_name in `cut -f2 -d'-' /mnt/disk4/BUDDHA/BMC-Bio/Segmental_deleltion_global/Vertebrate-RefSeqAnno.444.lst|sort -u`
    do
        j=`echo $species_name | sed 's/_/ /g'`
        LONGEST=`grep "$j" transcript_protein.tsv | sort -t$'\t' -nr -k7 | head -n1 | cut -f4`

        # Proceed only if LONGEST is found (not empty)
        if [ -n "$LONGEST" ]; then
            accession=`grep "$species_name" /mnt/disk4/BUDDHA/BMC-Bio/Segmental_deleltion_global/Vertebrate-RefSeqAnno.444.lst | cut -f1 -d'-'`
            datasets download genome accession $accession --include gff3 --no-progressbar --api-key c181c44f2dca87803a9c0f9a24a32f67a608 --filename "$species_name".zip
            rm -rf $species_name
            unzip -q -o "$species_name".zip -d $species_name
            rm "$species_name".zip

            gff2bed < $species_name/ncbi_dataset/data/$accession/genomic.gff  > $species_name/ncbi_dataset/data/$accession/genomic.gff.bed
            grep -w "$LONGEST" $species_name/ncbi_dataset/data/$accession/genomic.gff.bed | grep -P "\tCDS\t" | cut -f1-6 | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 7500000 -s > CDS_sorted/"$species_name"."$gene".bed
            
            [ -s CDS_sorted/"$species_name"."$gene".bed ] || rm -f CDS_sorted/"$species_name"."$gene".bed
            rm -r $species_name
            echo $species_name
        else
            echo "Skipping $species_name for $gene (No longest transcript found)"
        fi
done

cd CDS_sorted
ls -1 *.bed|cut -f1 -d'.'|sort -u > /mnt/disk4/BUDDHA/BMC-Bio/Segmental_deleltion_global/IL34_deletion/CDS_filtered/"$gene".species.lst
cd /mnt/disk4/BUDDHA/BMC-Bio/Segmental_deleltion_global/IL34_deletion
done

# get common species which have annoated in all genes

cd CDS_filtered
gene1=`cat ../../IL34_syntenic_genes.lst|cut -f1 -d'-'|tr '\n' '\t' |awk '{print $1}'`
gene2=`cat ../../IL34_syntenic_genes.lst|cut -f1 -d'-'|tr '\n' '\t' |awk '{print $2}'`
gene3=`cat ../../IL34_syntenic_genes.lst|cut -f1 -d'-'|tr '\n' '\t' |awk '{print $3}'`
gene4=`cat ../../IL34_syntenic_genes.lst|cut -f1 -d'-'|tr '\n' '\t' |awk '{print $4}'`
echo $gene1  $gene2  $gene3  $gene4
comm -12 <(sort -u "$gene1".species.lst) <(sort -u "$gene2".species.lst) | comm -12 - <(sort -u "$gene3".species.lst) | comm -12 - <(sort -u "$gene4".species.lst) > common.species.lst 
for gene in `cat ../../IL34_syntenic_genes.lst|cut -f1 -d'-'|tr '\n' ' '`
do
mkdir $gene
for species_name in `cat common.species.lst`
do
cp /mnt/disk4/BUDDHA/BMC-Bio/Segmental_deleltion_global/IL34_deletion/$gene/ncbi_dataset/data/CDS_sorted/"$species_name"."$gene".bed $gene
done
echo $gene
done

##### Get gene distance table
genes=($(cat ../../IL34_syntenic_genes.lst|cut -f1 -d'-'|tr '\n' ' '))

# Create an output file
output_file="gene_distances.tsv"

# Print header with all gene pair combinations
echo -ne "Species\t" > $output_file
for ((i=0; i<${#genes[@]}; i++)); do
    for ((j=i+1; j<${#genes[@]}; j++)); do
        echo -ne "${genes[i]}-${genes[j]}\t" >> $output_file
    done
done
echo "" >> $output_file

# Read species list from common.species.lst
species_list=($(cat common.species.lst))

# Process each species
for species in "${species_list[@]}"; do
    echo -ne "$species\t" >> $output_file  # Start a new row with species name
    
    for ((i=0; i<${#genes[@]}; i++)); do
        for ((j=i+1; j<${#genes[@]}; j++)); do
            gene1=${genes[i]}
            gene2=${genes[j]}
            
            file1="$gene1/$species.$gene1.bed"
            file2="$gene2/$species.$gene2.bed"
            
            if [[ -f "$file1" && -f "$file2" ]]; then
                # Run bedtools closest
                result=$(bedtools closest -a "$file1" -b "$file2" -d | awk 'NR==1 {print $1, $NF}')  # Get chromosome and distance
                
                chr1=$(echo "$result" | awk '{print $1}')  # Chromosome of gene1
                distance=$(echo "$result" | awk '{print $2}')  # Distance
                
                # Check if chromosomes are different (EBR case)
                if [[ "$chr1" != "$(awk 'NR==1 {print $1}' $file2)" ]]; then
                    distance="EBR"
                fi
            else
                distance="NA"  # If one gene is missing
            fi
            
            echo -ne "$distance\t" >> $output_file
        done
    done
    echo "" >> $output_file  # New row for the next species
done

echo "Pairwise gene distances saved in $output_file."
grep -v "EBR" gene_distances.tsv > gene_distances.sorted.tsv
cut -f1 gene_distances.sorted.tsv | tail -n +2 > species_list.txt
for i in `cat species_list.txt`
do
species=`echo $i |sed 's/_/ /g'`
order=`esearch -db taxonomy -query "$species" | efetch -format xml | xmllint --xpath "//LineageEx/Taxon[Rank='order']/ScientificName/text()" -`
class=`esearch -db taxonomy -query "$species" | efetch -format xml | xmllint --xpath "//LineageEx/Taxon[Rank='class']/ScientificName/text()" -`
echo -e "$class\t$order\t$i" >> taxonomy_info.tsv
echo $i
done




### find gaps
echo -e "Species\tNum_Gaps" > gap_summary.tsv  # Header for TSV file

for i in $(cat species_list.txt); do
    # Merge BED regions
    merged_region=$(cat $gene1/"$i"."$gene1".bed $gene4/"$i"."$gene4".bed | \
                    sort -k1,1 -k2,2n | \
                    bedtools merge -i - -d 100000000)

    # Extract chromosome ID, start, and stop positions
    id=$(echo "$merged_region" | cut -f1)
    chr_start=$(echo "$merged_region" | cut -f2)
    chr_stop=$(echo "$merged_region" | cut -f3)

    # Fetch FASTA and find gaps
    efetch -db nuccore -id $id -format fasta -chr_start $chr_start -chr_stop $chr_stop | \
        gzip > "$i"."$gene1"_"$gene4".fasta.gz

    klumpy find_gaps --fasta "$i"."$gene1"_"$gene4".fasta.gz 

    # Check if the gaps file exists, else set num_gaps = 0
    if [[ -f "$i"."$gene1"_"$gene4"_gaps.tsv ]]; then
        num_gaps=$(tail -n +2 "$i"."$gene1"_"$gene4"_gaps.tsv | wc -l)
    else
        num_gaps=0
    fi

    echo -e "$i\t$num_gaps" >> gap_summary.tsv
    echo -e "$i\t$num_gaps"
    rm "$i"."$gene1"_"$gene4".fasta.gz  "$i"."$gene1"_"$gene4"_gaps.tsv
done


head -n 1 gene_distances.sorted.tsv | awk '{print "Group\tOrder\tSpecies\tSpecies\tNum_Gaps\t"$0}'| awk 'BEGIN{OFS="\t"} $3==$4 && $4==$6 {print $2, $3, $5, $7, $8, $9, $10, $11, $12}' |sed 's/Species/Species_name/g;s/-/_/g' > final_output.tsv
paste taxonomy_info.tsv <(tail -n +2 gap_summary.tsv) <(tail -n +2 gene_distances.sorted.tsv) | awk 'BEGIN{OFS="\t"} $3==$4 && $4==$6 {print $2, $3, $5, $7, $8, $9, $10, $11, $12}' >> final_output.tsv
awk 'BEGIN{OFS="\t"} NR==1 {print $0, "Gene_status"; next} $1=="Squamata" {print $0, "Loss"; next} {print $0, "Intact"}' final_output.tsv > final_output_with_status.tsv
awk 'BEGIN{OFS="\t"} $1 != ""' final_output_with_status.tsv > final_output_with_status.filtered.tsv

cut -f2 final_output_with_status.filtered.tsv|tail -n+2> species_tree.lst
cd ..
mv ../IL34_human.lst ../IL34_human_sorted.lst ../IL34_syntenic_genes.lst .
cd  /mnt/disk4/BUDDHA/BMC-Bio/Segmental_deleltion_global/

##################### get gc content, gc strenches and sequence divergence ./get_GC_content_seq_divergence_of_53_genes.sh
rm -r GC_content_seq_divergence
mkdir GC_content_seq_divergence
echo -e "Gene_name\tSample_size\tGC_content\tGC_Stretch\tSeq_iden" > GC_content_seq_divergence/Gene_sequence_properties.tsv
for i in `cat Complete_loss.lst`
do
gene_id=`echo $i|cut -f1 -d':'`
gene=`echo $i|cut -f2 -d':'`
cut -f2,4,5 $gene_id/ncbi_dataset/data/transcript_protein.tsv|sed 's/ /_/g;s/\t/-/g'|tail -n+2 > GC_content_seq_divergence/"$gene".transcript_info.lst
cd GC_content_seq_divergence/
for i in `cat "$gene".transcript_info.lst|sort -u`
do
species_name=`echo $i|cut -f1 -d'-'`
accession=`echo $i|cut -f2 -d'-'`
length=`echo $i|cut -f3 -d'-'`
efetch -db nuccore -id "$accession" -format fasta_cds_na | sed "s/^>\(.*\) \[gene=.*$/>\1 $species_name:$accession:$length/"|sed 's/^>\S*_cds_\S*_\S*_\S* />/' >> "$gene".fa
done
cdskit aggregate --seqfile "$gene".fa --outfile "$gene".longest.fa --expression ":.*"
perl ../mulfifasta.pl "$gene".longest.fa |cut -f1 -d':'|seqtk seq > "$gene".valid.longest.fa
Sample_size=`grep ">" "$gene".valid.longest.fa |wc -l`
GC_content_mean_sd=`seqkit fx2tab --name --gc "$gene".valid.longest.fa| awk '{sum+=$2; sumsq+=$2*$2; n+=1} END {mean=sum/n; stddev=sqrt((sumsq - (sum^2)/n)/n); print  mean, "+/-", stddev}'`
GC_Stretch_mean_sd=`perl ../GC_Stretch_finder.pl "$gene".valid.longest.fa | cut -f5 -d' ' | awk '{sum+=$1; sumsq+=$1*$1; n+=1} 
     END {
         mean=sum/n; 
         stddev=sqrt((sumsq - (sum^2)/n)/n); 
         print mean, "+/-", stddev
     }'`
mkdir "$gene"_query/
faSplit byname "$gene".valid.longest.fa "$gene"_query/
echo "species ident"|sed 's/ /\t/g' > "$gene".identity_wrt_human.out
cd  "$gene"_query/
for i in *.fa
do
n=`echo $i|sed 's/\.fa//g'`
needle -asequence Homo_sapiens.fa -bsequence $i -outfile Homo_sapiens-"$n".needle.out -gapopen 10.0 -gapextend 0.5
sim=`grep "Similarity" Homo_sapiens-"$n".needle.out |awk '{print $4}'|sed 's/(//g'|sed 's/)//g'|sed 's/%//g'`
ident=`grep "Identity" Homo_sapiens-"$n".needle.out |awk '{print $4}'|sed 's/(//g'|sed 's/)//g'|sed 's/%//g'`
gaps=`grep "Gaps" Homo_sapiens-"$n".needle.out |awk '{print $5}'|sed 's/(//g'|sed 's/)//g'|sed 's/%//g'`
echo $n $ident|sed 's/ /\t/g' >> ../"$gene".identity_wrt_human.out
done
cd ..
Seq_iden_mean_sd=`cat "$gene".identity_wrt_human.out| awk '{sum+=$2; sumsq+=$2*$2; n+=1} END {mean=sum/n; stddev=sqrt((sumsq - (sum^2)/n)/n); print  mean, "+/-", stddev}'`

cd ..
echo -e "$gene\t$Sample_size\t$GC_content_mean_sd\t$GC_Stretch_mean_sd\t$Seq_iden_mean_sd" > GC_content_seq_divergence/Gene_sequence_properties.tsv
echo $gene
done

#################################### HUMAN and CHICKEN deletion fst ########################
#https://asia.ensembl.org/info/genome/compara/analyses.html
#chicken 27 pairwise alignment
for i in `rsync rsync://ftp.ebi.ac.uk/ensemblorg/pub/current_maf/ensembl-compara/pairwise_alignments/ | grep 'ggal_bgalgal1.mat.broiler.grcg7b.v.' | awk '{print $NF}'`
do
wget https://ftp.ensembl.org/pub/current_maf/ensembl-compara/pairwise_alignments/$i
tar -xvzf $i
cd ${i%.tar.gz}
focal_species=`head "$(ls -1 ggal_bgalgal1.mat.broiler.grcg7b.v.*.maf | head -1)" | grep "^s" | cut -f2 -d' ' | grep -v "gallus_gallus" | cut -f1 -d'.'|sort -u`
for maf in `ls -1 ggal_bgalgal1.mat.broiler.grcg7b.v.*.maf`
do
sed -i 's/^a# /a\n# /' $maf
maffilter param=/media/lokdeep/sdf/BUDDHA/BMCBIO/Squamates_segmental_deletion/maffilter.bpp DATA=$maf focal_species=$focal_species
done
cat *coordinates.txt |grep -v "gallus_gallus.chr"|awk '{print $1"\t"$3"\t"$4}'|sort -k1,1 -k2,2n|bedtools merge -i "stdin" -d 5 > ../"$focal_species".chicken.bed
cd ..
rm -r ${i%.tar.gz} $i
echo $focal_species
done


wget https://ftp.ensembl.org/pub/release-113/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.113.gtf.gz
wget https://ftp.ensembl.org/pub/release-113/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-113/regulation/gallus_gallus/GRCg7b/annotation/Gallus_gallus.GRCg7b.regulatory_features.v113.gff3.gz
wget https://ftp.ensembl.org/pub/release-113/regulation/gallus_gallus/GRCg7b/annotation/Gallus_gallus.GRCg7b.EMARs.v113.gff.gz
#location of Epigenetically Modified Accessible Regions (EMARs)
#https://asia.ensembl.org/info/genome/funcgen/data/regulatory-features.html
samtools faidx Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa
cut -f1,2 Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.fai > Gallus_gallus.chrom.sizes
bedtools makewindows -g Gallus_gallus.chrom.sizes -w 100 > Gallus_gallus.100bp.windows.bed
bedtools nuc -fi Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -bed Gallus_gallus.100bp.windows.bed > Gallus_gallus.100bp.gc_content.txt
awk 'NR>1 {print $1, $2, $3, $5}' Gallus_gallus.100bp.gc_content.txt|sed 's/ /\t/g' > Gallus_gallus.100bp.gc_content.bed
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_feature_table.txt.gz
zcat GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_feature_table.txt.gz|awk -F '\t' '$1=="gene" && $2=="protein_coding" {print $6"\t"$8"\t"$9"\t"$15}' > Gallus_gallus.gene.bed
cut -f1-3,11 Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.repeatmasker.bed> ../Gallus_gallus.repeat.bed
gff2bed < Gallus_gallus.GRCg7b.EMARs.v113.gff > Gallus_gallus.GRCg7b.EMARs.v113.gff.bed
cut -f1-3,8 Gallus_gallus.GRCg7b.EMARs.v113.gff.bed> Gallus_gallus.EMAR.bed
gunzip Gallus_gallus.GRCg7b.regulatory_features.v113.gff3.gz
cut -f1-3,8 Gallus_gallus.GRCg7b.regulatory_features.v113.gff3.bed >  Gallus_gallus.regulatory_features.bed

for window in 1000 5000 10000 20000 50000 100000 
do
win_kb=$((window / 1000))kb
bedtools makewindows -g Gallus_gallus.chrom.sizes -w $window > Gallus_gallus.${win_kb}.windows.bed
for i in *.chicken.bed
do
sp=$(echo $i | sed 's/\.chicken\.bed//g' | cut -f1,2 -d'_')
bedtools intersect -a Gallus_gallus.${win_kb}.windows.bed -b $i -wao | awk -v sp=$sp -v win=${win_kb} '{print $0"\t"sp"\t"win}' >> ${win_kb}.overlap.bed
echo "Processed $sp with ${win_kb} windows"
done
bedtools groupby -i ${win_kb}.overlap.bed -g 1,2,3,8 -c 7 -o sum >> ${win_kb}.merged_overlap.bed
bedtools intersect -a ${win_kb}.merged_overlap.bed -b Gallus_gallus.gene.bed -wa -wb >> ${win_kb}.gene_annotated.overlap.bed
rm ${win_kb}.overlap.bed ${win_kb}.merged_overlap.bed
done



sp1="Homo sapiens"
sp2="Gallus gallus"
sp1sp2=`echo "${sp1}-${sp2}.gene_orthologs.tsv"|sed 's/ /_/g'`
echo -e "Gene\t${sp1}_LONGEST_acc\t${sp1}_LONGEST_length\tGene\t${sp2}_LONGEST_acc\t${sp2}_LONGEST_length" > $sp1sp2

# Ensure output FASTA files are empty before starting
> "${sp1}.cds.fasta"
> "${sp2}.cds.fasta"

while IFS=':' read -r id gene; do
    transcript_file="$id/ncbi_dataset/data/transcript_protein.tsv"
    
    if [[ ! -f "$transcript_file" ]]; then
        echo "Warning: File $transcript_file not found for gene $gene. Skipping..."
        continue
    fi
    
    # Extract longest transcript for species 1
    sp1_LONGEST_acc=$(grep "$sp1" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f4)
    sp1_LONGEST_length=$(grep "$sp1" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f5)
    sp1_gene=$(grep "$sp1" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f3)
    
    # Extract longest transcript for species 2
    sp2_LONGEST_acc=$(grep "$sp2" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f4)
    sp2_LONGEST_length=$(grep "$sp2" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f5)
    sp2_gene=$(grep "$sp2" "$transcript_file" | sort -t$'\t' -nr -k7 | head -n1 | cut -f3)    
   
    
    echo -e "$sp1_gene\t$sp1_LONGEST_acc\t$sp1_LONGEST_length\t$sp2_gene\t$sp2_LONGEST_acc\t$sp2_LONGEST_length" >> $sp1sp2
done < Complete_loss.lst

cut -f4 Homo_sapiens-Gallus_gallus.gene_orthologs.tsv|sort -u |tr '\n' ' '
#chicken orthologs
#C11orf74==IFTAP
#HRASLS==PLAAT1
#C10orf128==TMEM273
echo -e "Gene_name\tSpecies_name\tNumber_of_unaligned_1kb_windows\tGene_length" > Number_of_unaligned_1kb_windows.tsv
for gene in ADAP2 ADRA1B ARL11 C11orf74 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HRASLS HSD17B1 IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 C10orf128 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
for species in `cut -f4 1kb.gene_annotated.overlap.bed|sort -u`
do
wind_no_align=$(awk -v sp="$species" -v g="$gene" '$4 == sp && $9 == g && $5 == 0' 1kb.gene_annotated.overlap.bed | cut -f1-4 | sort -u | wc -l)
gene_length=$(awk -v sp="$species" -v g="$gene" '$4 == sp && $9 == g' 1kb.gene_annotated.overlap.bed|cut -f6-9|bedtools merge -d 300000|awk '{print ($3-$2)/1000}')
echo -e "$gene\t$species\t$wind_no_align\t$gene_length" >> Number_of_unaligned_1kb_windows.tsv
echo -e "$gene\t$species\t$wind_no_align\t$gene_length" 
done
done
cut -f4 1kb.gene_annotated.overlap.bed|sort -u > Species.lst
for i in `cat  Species.lst`; do species=`echo $i |sed 's/_/ /g'`; order=`esearch -db taxonomy -query "$species" | efetch -format xml | xmllint --xpath "//LineageEx/Taxon[Rank='order']/ScientificName/text()" -`; class=`esearch -db taxonomy -query "$species" | efetch -format xml | xmllint --xpath "//LineageEx/Taxon[Rank='class']/ScientificName/text()" -`; echo -e "$class\t$order\t$i" >> taxonomy_info.tsv; echo $i; echo $i $order $class; done

echo -e "Gene_name\tClass\tOrder\tSpecies_name\tNumber_of_unaligned_1kb_windows\tGene_length" > Final.Number_of_unaligned_1kb_windows.tsv
for gene in ADAP2 ADRA1B ARL11 C11orf74 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HRASLS HSD17B1 IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 C10orf128 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
for species in `cut -f2 Number_of_unaligned_1kb_windows.tsv|tail -n+2|sort -u`
do
class=`awk -v sp="$species" '$3 == sp {print $1}' taxonomy_info.tsv`
order=`awk -v sp="$species" '$3 == sp {print $2}' taxonomy_info.tsv`
wind_no_align=$(awk -v sp="$species" -v g="$gene" '$2 == sp && $1 == g {print $3}' Number_of_unaligned_1kb_windows.tsv)
Gene_length=$(awk -v sp="$species" -v g="$gene" '$2 == sp && $1 == g {print $4}' Number_of_unaligned_1kb_windows.tsv)
echo -e "$gene\t$class\t$order\t$species\t$wind_no_align\t$Gene_length" >> Final.Number_of_unaligned_1kb_windows.tsv
echo -e "$gene\t$class\t$order\t$species\t$wind_no_align\t$Gene_length"
done
done




# human 171 pairwise alignment GRCh38 https://ftp.ensembl.org/pub/current_maf/ensembl-compara/pairwise_alignments/

for i in `rsync rsync://ftp.ebi.ac.uk/ensemblorg/pub/current_maf/ensembl-compara/pairwise_alignments/ | grep 'hsap_grch38.v.' | awk '{print $NF}'`
do
wget https://ftp.ensembl.org/pub/current_maf/ensembl-compara/pairwise_alignments/$i
tar -xvzf $i
cd ${i%.tar.gz} 
focal_species=`head "$(ls -1 hsap_grch38.v.*.maf | head -1)" | grep "^s" | cut -f2 -d' ' | grep -v "homo_sapiens" | cut -f1 -d'.'|sort -u`
for maf in `ls -1 hsap_grch38.v.*.maf`
do
sed -i 's/^a# /a\n# /' $maf
maffilter param=/media/lokdeep/sdf/BUDDHA/BMCBIO/Squamates_segmental_deletion/maffilter_human.bpp DATA=$maf focal_species=$focal_species
done
cat *coordinates.txt |grep -v "homo_sapiens.chr"|awk '{print $1"\t"$3"\t"$4}'|sort -k1,1 -k2,2n|bedtools merge -i "stdin" -d 5 > ../"$focal_species".homo_sapiens.bed
cd /media/lokdeep/sdf/BUDDHA/BMCBIO/Squamates_segmental_deletion/human/
rm -r ${i%.tar.gz} $i
echo $focal_species
echo $focal_species >> species.lst
done

wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.*
samtools faidx Homo_sapiens.GRCh38.dna.toplevel.fa
cut -f1,2 Homo_sapiens.GRCh38.dna.toplevel.fa.fai > Homo_sapiens.chrom.sizes
bedtools makewindows -g Homo_sapiens.chrom.sizes -w 100 > Homo_sapiens.100bp.windows.bed
bedtools nuc -fi Homo_sapiens.GRCh38.dna.toplevel.fa -bed Homo_sapiens.100bp.windows.bed > Homo_sapiens.100bp.gc_content.txt
awk 'NR>1 {print $1, $2, $3, $5}' Homo_sapiens.100bp.gc_content.txt|sed 's/ /\t/g' > Homo_sapiens.100bp.gc_content.bed
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_feature_table.txt.gz
zcat GCF_000001405.40_GRCh38.p14_feature_table.txt.gz|awk -F '\t' '$1=="gene" && $2=="protein_coding" {print $6"\t"$8"\t"$9"\t"$15}' > Homo_sapiens.gene.bed

for window in 1000
do
win_kb=$((window / 1000))kb
bedtools makewindows -g Homo_sapiens.chrom.sizes -w $window > Homo_sapiens.${win_kb}.windows.bed
for i in *.homo_sapiens.bed
do
sp=$(echo $i | sed 's/\.homo_sapiens\.bed//g' | cut -f1,2 -d'_')
bedtools intersect -a Homo_sapiens.${win_kb}.windows.bed -b $i -wao | awk -v sp=$sp -v win=${win_kb} '{print $0"\t"sp"\t"win}' >> ${win_kb}.overlap.bed
echo "Processed $sp with ${win_kb} windows"
done
bedtools groupby -i ${win_kb}.overlap.bed -g 1,2,3,8 -c 7 -o sum >> ${win_kb}.merged_overlap.bed
bedtools intersect -a ${win_kb}.merged_overlap.bed -b Homo_sapiens.gene.bed -wa -wb >> ${win_kb}.gene_annotated.overlap.bed
rm ${win_kb}.overlap.bed ${win_kb}.merged_overlap.bed
done
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
grep -w "$gene" Homo_sapiens.gene.bed >> Candidate_genes.bed
echo $gene
done

# count the 1kb window which shows no alignment in chain
echo -e "Gene_name\tSpecies_name\tNumber_of_unaligned_1kb_windows" > Number_of_unaligned_1kb_windows.tsv

for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
for species in `cut -f4 1kb.gene_annotated.overlap.bed|sort -u`
do
wind_no_align=$(awk -v sp="$species" -v g="$gene" '$4 == sp && $9 == g && $5 == 0' 1kb.gene_annotated.overlap.bed | cut -f1-4 | sort -u | wc -l)
echo -e "$gene\t$species\t$wind_no_align" >> Number_of_unaligned_1kb_windows.tsv
echo -e "$gene\t$species\t$wind_no_align" 
done
done
cut -f4 1kb.gene_annotated.overlap.bed|sort -u > Species.lst
for i in `cat  Species.lst`; do species=`echo $i |sed 's/_/ /g'`; order=`esearch -db taxonomy -query "$species" | efetch -format xml | xmllint --xpath "//LineageEx/Taxon[Rank='order']/ScientificName/text()" -`; class=`esearch -db taxonomy -query "$species" | efetch -format xml | xmllint --xpath "//LineageEx/Taxon[Rank='class']/ScientificName/text()" -`; echo -e "$class\t$order\t$i" >> taxonomy_info.tsv; echo $i; echo $i $order $class; done

head -1 HSD17B1.Number_of_unaligned_1kb_windows.tsv > Number_of_unaligned_1kb_windows.tsv
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
tail -n +2 "${gene}.Number_of_unaligned_1kb_windows.tsv" | awk 'NF == 4' >> Number_of_unaligned_1kb_windows.tsv
done


echo -e "Gene_name\tClass\tOrder\tSpecies_name\tNumber_of_unaligned_1kb_windows\tGene_length" > Final.Number_of_unaligned_1kb_windows.tsv
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
for species in `cut -f2 Number_of_unaligned_1kb_windows.tsv|tail -n+2|sort -u`
do
class=`awk -v sp="$species" '$3 == sp {print $1}' taxonomy_info.tsv`
order=`awk -v sp="$species" '$3 == sp {print $2}' taxonomy_info.tsv`
wind_no_align=$(awk -v sp="$species" -v g="$gene" '$2 == sp && $1 == g {print $3}' Number_of_unaligned_1kb_windows.tsv)
Gene_length=$(awk -v sp="$species" -v g="$gene" '$2 == sp && $1 == g {print $4}' Number_of_unaligned_1kb_windows.tsv)
echo -e "$gene\t$class\t$order\t$species\t$wind_no_align\t$Gene_length" >> Final.Number_of_unaligned_1kb_windows.tsv
echo -e "$gene\t$class\t$order\t$species\t$wind_no_align\t$Gene_length"
done
done


for i in `cut -f1 Candidate_genes.bed|sort -u`; do  awk -v ch="$i" '$1==ch {print $0}' Gallus_gallus.chrom.sizes; done|sort -nk1,1 > Candidate_gene.Gallus_gallus.chrom.sizes


##gene cluseter within 2Mb region
sort -k1,1 -k2,2n genespace.candidate_gene.bed > sorted.bed
bedtools cluster -d 2000000 -i sorted.bed > clustered.bed
sort -nk5,5 clustered.bed

#REGION_1 1	167777273	171335742
#1	167777273	167813521	SIAH3	10
#1	167910600	167959913	LCP1	10
#1	167993956	168013779	RUBCNL	10
#1	168080147	168204194	LRCH1	10
#1	169146947	169152771	ARL11	10
#1	169704013	169713454	DLEU7	10
#1	169978669	170020436	INTS6	10
#1	171141631	171157413	STOML3	10
#1	171213674	171335742	FREM2	10

#REGION_2
#4	79231547	79261514	MAN2B2	30
#4	79387287	79452683	TBC1D14	30
#4	80670212	80677171	GPR78	30
#4	81719570	81754068	GRK4	30
#4	81949813	81957766	TNIP2	30

#REGION_3
#6	17564778	17569174	ELOVL3	33
#6	19374587	19394523	C10orf128	33

#REGION_4
#8	4165167	4244193	NOTCH2	35
#8	4530786	4697586	DNM3	35

#REGION_5
#9	12827756	12830575	HRASLS	36
#9	13401200	13408002	UTS2B	36

#REGION_6 
#11	1354457	1365252	SLC9A5	12
#11	1738296	1742625	IL34	12
#REGION_7
#12	18269804	18285446	IL5RA	13
#12	19458124	19461512	CAV3	13






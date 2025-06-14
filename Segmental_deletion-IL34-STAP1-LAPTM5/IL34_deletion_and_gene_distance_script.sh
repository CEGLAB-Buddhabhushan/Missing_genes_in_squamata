conda activate ncbi_datasets
rm -r IL34_deletion
mkdir  IL34_deletion
cd IL34_deletion
mkdir CDS_filtered
for i in FCSK-197258 COG4-25839 SF3B3-23450 MTSS2-92154 VAC14-55697 HYDIN-54768
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

for species_name in `cut -f9 gene_product.tsv|tail -n+2|sed 's/ /_/g'|sort -u`
do
j=`echo $species_name |sed 's/_/ /g'`
LONGEST=`grep "$j" transcript_protein.tsv|sort -t$'\t' -nr -k7|head -n1|cut -f4`
grep -w "$LONGEST" gene_product.tsv|awk -F'\t' '{print $18"\t"$21"\t"$22"\t"$9"\t"$19"\t"$24}'|sed 's/plus/\+/g;s/minus/\-/g;s/ /_/g'|sort -k1,1 -k2,2n|bedtools merge -i stdin -d 7500000 -s > CDS_sorted/"$species_name"."$gene".bed
[ -s CDS_sorted/"$species_name"."$gene".bed ] || rm -f CDS_sorted/"$species_name"."$gene".bed
echo $species_name
done

cd CDS_sorted
ls -1 *.bed|cut -f1 -d'.'|sort -u > /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/CDS_filtered/"$gene".species.lst
cd /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion
done

# get common species which have annoated in all genes
#FCSK-197258 COG4-25839 SF3B3-23450 MTSS2-92154 VAC14-55697 HYDIN-54768
cd CDS_filtered
comm -12 <(sort -u FCSK.species.lst) <(sort -u COG4.species.lst) | comm -12 - <(sort -u SF3B3.species.lst) | comm -12 - <(sort -u MTSS2.species.lst) | comm -12 - <(sort -u VAC14.species.lst)| comm -12 - <(sort -u HYDIN.species.lst) > common.species.lst 
for gene in FCSK COG4 SF3B3 MTSS2 VAC14 HYDIN
do
mkdir $gene
for species_name in `cat common.species.lst`
do
cp /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/$gene/ncbi_dataset/data/CDS_sorted/"$species_name"."$gene".bed $gene
done
echo $gene
done

##### Get gene distance table
genes=("FCSK" "COG4" "SF3B3" "MTSS2" "VAC14" "HYDIN")

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
conda activate klumpy_env
echo -e "Species\tNum_Gaps" > gap_summary.tsv  # Header for TSV file

for i in $(cat species_list.txt); do
    # Merge BED regions
    merged_region=$(cat FCSK/"$i".FCSK.bed HYDIN/"$i".HYDIN.bed | \
                    sort -k1,1 -k2,2n | \
                    bedtools merge -i - -d 100000000)

    # Extract chromosome ID, start, and stop positions
    id=$(echo "$merged_region" | cut -f1)
    chr_start=$(echo "$merged_region" | cut -f2)
    chr_stop=$(echo "$merged_region" | cut -f3)

    # Fetch FASTA and find gaps
    efetch -db nuccore -id $id -format fasta -chr_start $chr_start -chr_stop $chr_stop | \
        gzip > "$i".FCSK_HYDIN.fasta.gz

    klumpy find_gaps --fasta "$i".FCSK_HYDIN.fasta.gz 

    # Check if the gaps file exists, else set num_gaps = 0
    if [[ -f "$i".FCSK_HYDIN_gaps.tsv ]]; then
        num_gaps=$(tail -n +2 "$i".FCSK_HYDIN_gaps.tsv | wc -l)
    else
        num_gaps=0
    fi

    echo -e "$i\t$num_gaps" >> gap_summary.tsv
    echo -e "$i\t$num_gaps"
    rm "$i".FCSK_HYDIN.fasta.gz  "$i".FCSK_HYDIN_gaps.tsv
done


head -n 1 gene_distances.sorted.tsv | awk '{print "Group\tOrder\tSpecies_name\tSpecies\tNum_Gaps\t"$0}' > final_output.tsv
paste taxonomy_info.tsv <(tail -n +2 gap_summary.tsv) <(tail -n +2 gene_distances.sorted.tsv) >> final_output.tsv
head -1 final_output.tsv > final_output.no_gaps.tsv
awk '$5==0' final_output.tsv >> final_output.no_gaps.tsv


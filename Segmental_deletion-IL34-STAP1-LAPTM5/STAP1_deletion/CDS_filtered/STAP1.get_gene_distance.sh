conda activate ncbi_datasets
rm -r STAP1_deletion
mkdir  STAP1_deletion
cd STAP1_deletion
mkdir CDS_filtered
for i in EPHA5-2044 CENPC-1060 UBA6-55236 YTHDC1-91746
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

for species_name in `cut -f2 -d'-' /mnt/disk4/BUDDHA/BMC-Bio/Vertebrate-RefSeqAnno.444.lst|sort -u`
    do
        j=`echo $species_name | sed 's/_/ /g'`
        LONGEST=`grep "$j" transcript_protein.tsv | sort -t$'\t' -nr -k7 | head -n1 | cut -f4`

        # Proceed only if LONGEST is found (not empty)
        if [ -n "$LONGEST" ]; then
            accession=`grep "$species_name" /mnt/disk4/BUDDHA/BMC-Bio/Vertebrate-RefSeqAnno.444.lst | cut -f1 -d'-'`
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
ls -1 *.bed|cut -f1 -d'.'|sort -u > /mnt/disk4/BUDDHA/BMC-Bio/STAP1_deletion/CDS_filtered/"$gene".species.lst
cd /mnt/disk4/BUDDHA/BMC-Bio/STAP1_deletion
done

# get common species which have annoated in all genes

cd CDS_filtered
comm -12 <(sort -u EPHA5.species.lst) <(sort -u CENPC.species.lst) | comm -12 - <(sort -u UBA6.species.lst) | comm -12 - <(sort -u YTHDC1.species.lst) > common.species.lst 
for gene in EPHA5 CENPC UBA6 YTHDC1
do
mkdir $gene
for species_name in `cat common.species.lst`
do
cp /mnt/disk4/BUDDHA/BMC-Bio/STAP1_deletion/$gene/ncbi_dataset/data/CDS_sorted/"$species_name"."$gene".bed $gene
done
echo $gene
done

##### Get gene distance table
genes=("EPHA5" "CENPC" "UBA6" "YTHDC1")

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
    merged_region=$(cat EPHA5/"$i".EPHA5.bed YTHDC1/"$i".YTHDC1.bed | \
                    sort -k1,1 -k2,2n | \
                    bedtools merge -i - -d 100000000)

    # Extract chromosome ID, start, and stop positions
    id=$(echo "$merged_region" | cut -f1)
    chr_start=$(echo "$merged_region" | cut -f2)
    chr_stop=$(echo "$merged_region" | cut -f3)

    # Fetch FASTA and find gaps
    efetch -db nuccore -id $id -format fasta -chr_start $chr_start -chr_stop $chr_stop | \
        gzip > "$i".EPHA5_YTHDC1.fasta.gz

    klumpy find_gaps --fasta "$i".EPHA5_YTHDC1.fasta.gz 

    # Check if the gaps file exists, else set num_gaps = 0
    if [[ -f "$i".EPHA5_YTHDC1_gaps.tsv ]]; then
        num_gaps=$(tail -n +2 "$i".EPHA5_YTHDC1_gaps.tsv | wc -l)
    else
        num_gaps=0
    fi

    echo -e "$i\t$num_gaps" >> gap_summary.tsv
    echo -e "$i\t$num_gaps"
    rm "$i".EPHA5_YTHDC1.fasta.gz  "$i".EPHA5_YTHDC1_gaps.tsv
done


head -n 1 gene_distances.sorted.tsv | awk '{print "Group\tOrder\tSpecies\tSpecies\tNum_Gaps\t"$0}'| awk 'BEGIN{OFS="\t"} $3==$4 && $4==$6 {print $2, $3, $5, $7, $8, $9, $10, $11, $12}' |sed 's/Species/Species_name/g;s/-/_/g' > final_output.tsv
paste taxonomy_info.tsv <(tail -n +2 gap_summary.tsv) <(tail -n +2 gene_distances.sorted.tsv) | awk 'BEGIN{OFS="\t"} $3==$4 && $4==$6 {print $2, $3, $5, $7, $8, $9, $10, $11, $12}' >> final_output.tsv
awk 'BEGIN{OFS="\t"} NR==1 {print $0, "Gene_status"; next} $1=="Squamata" {print $0, "Loss"; next} {print $0, "Intact"}' final_output.tsv > final_output_with_status.tsv
awk 'BEGIN{OFS="\t"} $1 != ""' final_output_with_status.tsv > final_output_with_status.filtered.tsv

cut -f2 final_output_with_status.filtered.tsv|tail -n+2> species_tree.lst
#get species tree from time tree website



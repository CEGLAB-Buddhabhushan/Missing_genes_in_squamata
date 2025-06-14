#get the list of chromosome leve assebly of vertbrate genome
#https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=7742&reference_only=true&annotated_only=true&refseq_annotation=true&assembly_level=2:3

cut -f1,3 Vertebrate-RefSeqAnno.467_genomes.tsv|cut -f1,2 -d' '|sed 's/\t/-/g;s/ /_/g' > Vertebrate-RefSeqAnno.467_genomes.reformated.lst
cut -f2 -d'-' Vertebrate-RefSeqAnno.467_genomes.reformated.lst|tail -n+2 |sort -u > Species.lst
##get species tree from time tree website
# remove internodes
sed -e "s/'[^()]*'//g" Species.nwk > Species.no_internodes.nwk
# get the species found on timetree webside ##444
sed 's/(/\n/g' Species.no_internodes.nwk |sed 's/)/\n/g' |sed 's/;/\n/g' |sed 's/:/\n/g' |sed 's/,/\n/g' |grep "^[A-Z]"|sort -u > Species.no_internodes.nwk.txt
comm -12 Species.no_internodes.nwk.txt Species.lst > species.filtered.lst
## download species tree for filtered species
# fulter genome for cosideration in phytools
for i in `cat species.filtered.lst`; do grep -w "$i" Vertebrate-RefSeqAnno.467_genomes.reformated.lst; done > Vertebrate-RefSeqAnno.444.lst
# Download the genomic information calcualte the paiwise gene distance
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

for i in `cat Vertebrate-RefSeqAnno.444.lst|head -1`
do
accession=`echo $i |cut -f1 -d'-'`
species_name=`echo $i |cut -f2 -d'-'`
datasets download genome accession $accession --include genome,gff3 --no-progressbar --api-key c181c44f2dca87803a9c0f9a24a32f67a608 --filename "$i".zip
rm -rf $i
unzip -q -o "$i".zip -d $i
rm "$i".zip
cd $i/ncbi_dataset/data/$accession/
grep -P "\tCDS\t" genomic.gff | awk 'BEGIN {OFS="\t"} $3 == "CDS" { 
    chrom=$1; start=$4-1; end=$5; 
    gene = gensub(/.*gene=([^;]+).*/, "\\1", "g", $9); 
    key = chrom "_" start "_" end; 
    if (!seen[key]) { print chrom, start, end, gene, 0, $7; seen[key]=1 } 
}' > "$i"_cds.bed
cd /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/
echo $Species_name
done

grep -P "\tCDS\t" genomic.gff | awk 'BEGIN {OFS="\t"} $3 == "CDS" { 
    chrom=$1; start=$4-1; end=$5; 
    gene = gensub(/.*gene=([^;]+).*/, "\\1", "g", $9); 
    key = chrom "_" start "_" end; 
    if (!seen[key]) { print chrom, start, end, gene, 0, $7; seen[key]=1 } 
}' |grep -i -w "il34"|sort -t$'\t' -nr -k2|head -n1 

###############
mkdir filtered
for i in FCSK-197258 COG4-25839 SF3B3-23450 MTSS2-92154 VAC14-55697 HYDIN-54768
do
gene=`echo $i | cut -f1 -d'-'`
gene_id=`echo $i | cut -f2 -d'-'`
    
datasets download gene gene-id $gene_id --include product-report --ortholog all --api-key --no-progressbar --api-key c181c44f2dca87803a9c0f9a24a32f67a608 --filename "$gene".zip
rm -rf "$gene"
unzip -q -o "$gene".zip -d "$gene"
rm "$gene".zip
cd "$gene/ncbi_dataset/data"
dataformat tsv gene --inputfile data_report.jsonl --fields annotation-assembly-accession,annotation-genomic-range-accession,annotation-genomic-range-range-start,annotation-genomic-range-range-stop,symbol,protein-count,orientation > data_report.sorted.tsv
    
mkdir -p sorted
    
for i in `cat ../../../Vertebrate-RefSeqAnno.444.lst`
do
accession=`echo $i | cut -f1 -d'-'`
species_name=`echo $i | cut -f2 -d'-'`
        
grep -w "$accession" data_report.sorted.tsv | \
awk -v gene="$gene" '{print $2, $3, $4, gene, ".", $7}' | \
sed 's/plus/\+/g;s/minus/\-/g;s/ /\t/g' > sorted/"$species_name"."$gene".bed
[ -s sorted/"$species_name"."$gene".bed ] || rm -f sorted/"$species_name"."$gene".bed
echo "$accession $species_name"
done
cd sorted
ls -1 *.bed|cut -f1 -d'.'|sort -u > /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/filtered/"$gene".species.lst
cd /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion
done

# get common species which have annoated in all genes
cd filtered
comm -12 <(sort -u COG4.species.lst) <(sort -u FCSK.species.lst) | comm -12 - <(sort -u HYDIN.species.lst) | comm -12 - <(sort -u MTSS2.species.lst) | comm -12 - <(sort -u SF3B3.species.lst) | comm -12 - <(sort -u VAC14.species.lst) > common.species.lst 
for gene in FCSK COG4 SF3B3 MTSS2 VAC14 HYDIN
do
mkdir $gene
for species_name in `cat common.species.lst`
do
cp /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/$gene/ncbi_dataset/data/sorted/"$species_name"."$gene".bed $gene
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
cut -f1  gene_distances.sorted.tsv|tail -n+2 |sort -u > gene_distances.sorted.species.lst
## get species tree for gene_distances.sorted.species.lst from time tree website
sed -e "s/'[^()]*'//g" gene_distances.sorted.species.nwk > gene_distances.sorted.species.no_internodes.nwk






######### Get bed for all species

rm -r All_filtered
mkdir All_filtered
for i in FCSK-197258 COG4-25839 SF3B3-23450 MTSS2-92154 VAC14-55697 HYDIN-54768
do
gene=`echo $i | cut -f1 -d'-'`
gene_id=`echo $i | cut -f2 -d'-'`
cd "$gene/ncbi_dataset/data"
dataformat tsv gene --inputfile data_report.jsonl --fields annotation-assembly-accession,annotation-genomic-range-accession,annotation-genomic-range-range-start,annotation-genomic-range-range-stop,symbol,protein-count,orientation,tax-name|sed 's/ /_/g' > data_report.tsv
    

rm -r All_sorted
mkdir -p All_sorted
for species_name in `cut -f8 data_report.tsv|tail -n+2`
do
grep -w "$species_name" data_report.tsv | \
awk -v gene="$gene" '{print $2, $3, $4, gene, ".", $7}' | \
sed 's/plus/\+/g;s/minus/\-/g;s/ /\t/g' > All_sorted/"$species_name"."$gene".bed
[ -s All_sorted/"$species_name"."$gene".bed ] || rm -f All_sorted/"$species_name"."$gene".bed
echo "$gene $species_name"
done
cd All_sorted
ls -1 *.bed|cut -f1 -d'.'|sort -u > /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/All_filtered/"$gene".species.lst
cd /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion
done

# get common species which have annoated in all genes
cd All_filtered
comm -12 <(sort -u COG4.species.lst) <(sort -u FCSK.species.lst) | comm -12 - <(sort -u HYDIN.species.lst) | comm -12 - <(sort -u MTSS2.species.lst) | comm -12 - <(sort -u SF3B3.species.lst) | comm -12 - <(sort -u VAC14.species.lst) > common.species.lst 
for gene in FCSK COG4 SF3B3 MTSS2 VAC14 HYDIN
do
mkdir $gene
for species_name in `cat common.species.lst`
do
cp /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/$gene/ncbi_dataset/data/All_sorted/"$species_name"."$gene".bed $gene
done
echo $gene
done
##### Get gene distance table
genes=("FCSK" "COG4" "SF3B3" "MTSS2" "VAC14" "HYDIN")

# Create an output file
output_file="gene_distances.All.tsv"

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
grep -v "EBR" gene_distances.All.tsv > gene_distances.All.sorted.tsv
## add class and order info
cut -f1 gene_distances.All.sorted.tsv | tail -n +2 > species_list.txt
for i in `cat species_list.txt`
do
species=`echo $i |sed 's/_/ /g'`
order=`esearch -db taxonomy -query "$species" | efetch -format xml | xmllint --xpath "//LineageEx/Taxon[Rank='order']/ScientificName/text()" -`
class=`esearch -db taxonomy -query "$species" | efetch -format xml | xmllint --xpath "//LineageEx/Taxon[Rank='class']/ScientificName/text()" -`
echo -e "$class\t$order\t$i" >> taxonomy_info.tsv
echo $i
done

head -n 1 gene_distances.All.sorted.tsv | awk '{print "Class\tOrder\tSpecies_name\t"$0}' > final_output.tsv
paste taxonomy_info.tsv <(tail -n +2 gene_distances.All.sorted.tsv) >> final_output.tsv


cut -f1   gene_distances.All.sorted.tsv|tail -n+2 |sort -u > gene_distances.All.sorted.species.lst
## get species tree for gene_distances.sorted.species.lst from time tree website
sed -e "s/'[^()]*'//g" gene_distances.sorted.species.nwk > gene_distances.sorted.species.no_internodes.nwk





################ BASED ON CDS ####################
rm -r CDS_filtered
mkdir CDS_filtered
for i in FCSK-197258 COG4-25839 SF3B3-23450 MTSS2-92154 VAC14-55697 HYDIN-54768
do
gene=`echo $i | cut -f1 -d'-'`
gene_id=`echo $i | cut -f2 -d'-'`
    
#datasets download gene gene-id $gene_id --include product-report --ortholog all --api-key --no-progressbar --api-key c181c44f2dca87803a9c0f9a24a32f67a608 --filename "$gene".zip
#rm -rf "$gene"
#unzip -q -o "$gene".zip -d "$gene"
#rm "$gene".zip
cd "$gene/ncbi_dataset/data"
dataformat tsv gene-product --inputfile product_report.jsonl  > gene_product.tsv
rm -r CDS_sorted
mkdir -p CDS_sorted

for species_name in `cat /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/Common_species_n310_withVertebrate-RefSeqAnno.lst|sort -u`
do
j=`echo $species_name |sed 's/_/ /g'`

grep -w "$j" gene_product.tsv | sort -nk33 | awk -F'\t' '{print $9"\t"$18"\t"$14"\t"$15"\t"$9"\t"0"\t"$24"\t"$25"\t"$26}'|awk -F'\t' 'NF==9 && !($0 ~ /\t\t/)' |head -1> CDS_sorted/"$species_name".tmp.bed
awk -F'\t' '{
    if ($7 == "plus") {
        genomic_start = $8 + $3 - 1;
        genomic_end = $8 + $4 - 1;
    } else {
        genomic_start = $9 - $4 + 1;
        genomic_end = $9 - $3 + 1;
    }
    print $2, genomic_start, genomic_end, $1, 0, ($7 == "plus" ? "+" : "-");
}' OFS='\t' CDS_sorted/"$species_name".tmp.bed | sed 's/plus/\+/g;s/minus/\-/g;s/ /_/g' > CDS_sorted/"$species_name"."$gene".bed
rm CDS_sorted/"$species_name".tmp.bed
[ -s CDS_sorted/"$species_name"."$gene".bed ] || rm -f CDS_sorted/"$species_name"."$gene".bed
#rm tmp.bed
echo $species_name
done

cd CDS_sorted
ls -1 *.bed|cut -f1 -d'.'|sort -u > /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/CDS_filtered/"$gene".species.lst
cd /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion
done

# get common species which have annoated in all genes
cd CDS_filtered
comm -12 <(sort -u COG4.species.lst) <(sort -u FCSK.species.lst) | comm -12 - <(sort -u HYDIN.species.lst) | comm -12 - <(sort -u MTSS2.species.lst) | comm -12 - <(sort -u SF3B3.species.lst) | comm -12 - <(sort -u VAC14.species.lst) > common.species.lst 
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

head -n 1 gene_distances.sorted.tsv | awk '{print "Class\tOrder\tSpecies_name\t"$0}' > final_output.tsv
paste taxonomy_info.tsv <(tail -n +2 gene_distances.sorted.tsv) >> final_output.tsv


########################## ALL CDS ##################
rm -r CDS_filtered
mkdir CDS_filtered
for i in FCSK-197258 COG4-25839 SF3B3-23450 MTSS2-92154 VAC14-55697 HYDIN-54768
do
gene=`echo $i | cut -f1 -d'-'`
gene_id=`echo $i | cut -f2 -d'-'`
    
#datasets download gene gene-id $gene_id --include product-report --ortholog all --no-progressbar --api-key c181c44f2dca87803a9c0f9a24a32f67a608 --filename "$gene".zip
#rm -rf "$gene"
#unzip -q -o "$gene".zip -d "$gene"
#rm "$gene".zip
cd "$gene/ncbi_dataset/data"
rm  gene_product.tsv
dataformat tsv gene-product --inputfile product_report.jsonl  > gene_product.tsv

rm -r CDS_sorted
mkdir -p CDS_sorted

for species_name in `cut -f9 gene_product.tsv|sed 's/ /_/g'|sort -u`
do
j=`echo $species_name |sed 's/_/ /g'`
grep -w "$j" gene_product.tsv | sort -nk33 | awk -F'\t' '{print $18"\t"$25"\t"$26"\t"$9"\t"0"\t"$24}'|awk -F'\t' 'NF==6 && !($0 ~ /\t\t/)' |head -1| sed 's/plus/\+/g;s/minus/\-/g;s/ /_/g'> CDS_sorted/"$species_name"."$gene".bed
[ -s CDS_sorted/"$species_name"."$gene".bed ] || rm -f CDS_sorted/"$species_name"."$gene".bed
echo $species_name
done

cd CDS_sorted
ls -1 *.bed|cut -f1 -d'.'|sort -u > /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion/CDS_filtered/"$gene".species.lst
cd /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_deletion
done

# get common species which have annoated in all genes
cd CDS_filtered
comm -12 <(sort -u COG4.species.lst) <(sort -u FCSK.species.lst) | comm -12 - <(sort -u HYDIN.species.lst) | comm -12 - <(sort -u MTSS2.species.lst) | comm -12 - <(sort -u SF3B3.species.lst) | comm -12 - <(sort -u VAC14.species.lst) > common.species.lst 
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

head -n 1 gene_distances.sorted.tsv | awk '{print "Class\tOrder\tSpecies_name\t"$0}' > final_output.tsv
paste taxonomy_info.tsv <(tail -n +2 gene_distances.sorted.tsv) >> final_output.tsv


### find gaps

echo -e "Species\tNum_Gaps" > gap_summary.tsv  # Header for TSV file

for i in $(cat species_list.txt); do
    # Merge BED regions
    merged_region=$(cat COG4/"$i".COG4.bed VAC14/"$i".VAC14.bed | \
                    sort -k1,1 -k2,2n | \
                    bedtools merge -i - -d 100000000)

    # Extract chromosome ID, start, and stop positions
    id=$(echo "$merged_region" | cut -f1)
    chr_start=$(echo "$merged_region" | cut -f2)
    chr_stop=$(echo "$merged_region" | cut -f3)

    # Fetch FASTA and find gaps
    efetch -db nuccore -id $id -format fasta -chr_start $chr_start -chr_stop $chr_stop | \
        gzip > "$i".COG4_VAC14.fasta.gz

    klumpy find_gaps --fasta "$i".COG4_VAC14.fasta.gz 

    # Check if the gaps file exists, else set num_gaps = 0
    if [[ -f "$i".COG4_VAC14_gaps.tsv ]]; then
        num_gaps=$(tail -n +2 "$i".COG4_VAC14_gaps.tsv | wc -l)
    else
        num_gaps=0
    fi

    echo -e "$i\t$num_gaps" >> gap_summary.tsv
    echo -e "$i\t$num_gaps"
    rm "$i".COG4_VAC14.fasta.gz 
done

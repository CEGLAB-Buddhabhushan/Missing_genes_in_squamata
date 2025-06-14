*********pacbio **********************************************************************************************
naja naja take chromosome co-ordinates 
chr=$(cat co_ordinates.txt | cut -f1 -d':') && esearch -db nucleotide -query "$chr" | efetch -format fasta > "Naja_naja.fa"

######________downlaoding data__________ ########################################### 
conda activate pfd

for i in `cat SRA.lst`
do
parallel-fastq-dump --sra-id $i --threads 16 --outdir . --gzip
done

######_________mapping_____minimap______ ##############

species2="Naja_naja"  

for fastq in *.fastq.gz
do
  /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/minimap2-2.28_x64-linux/minimap2 -ax map-pb -t 36 ../"$species2".fa $fastq | samtools sort -@36 -O BAM - > ${fastq%.fastq.gz}.bam
  samtools index ${fastq%.fastq.gz}.bam                   
done

samtools merge -@36 packbio_merged.bam *.bam 
samtools sort -@36 packbio_merged.bam -o packbio_merged.sorted.bam
samtools index -@36 packbio_merged.sorted.bam

#######________Klumpy____________________###############
klumpy scan_alignments --alignment_map packbio_merged.sorted.bam  --threads 32  --annotation /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Naja_naja/GCA_009733165.1_Nana_v5_genomic.gtf

klumpy find_gaps --fasta  /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Naja_naja/Naja_naja.fa

klumpy alignment_plot --alignment_map packbio_merged.sorted.bam  --reference CM019161.1 --candidates  packbio_merged.sorted_Candidate_Regions.tsv --window_size 10000 --window_step 5000 --color red --vertical_line_gaps --vertical_line_klumps --format svg --leftbound 768622 --rightbound 781041 --annotation /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Naja_naja/GCA_009733165.1_Nana_v5_genomic.gtf  --gap_file Naja_naja_gaps.tsv  --height 1500 --width 4000 --min_len 14000   --vertical_line_exons --number


BLAST ######################################################################################################################
zcat *fastq.gz |sed  -n '1~4s/^@/>/p;2~4p' > "$species2".fa
makeblastdb -in "$species".fa -out "$species".fa -dbtype nucl  done 

######________blast and store outfmt 1 and 7_____####### 
species2="Naja_naja" 
genome=/media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Naja_naja/PacBio/Naja_naja.fa 
for query in IL34_query_multispecies.fa Sphenodon_punctatus_IL34_Exons.fa
do
blastn -task blastn -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/outfmt1/blastn."$query"."$species2".e-value0.05.outfmt1.out -num_threads 32 -outfmt 1 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/outfmt7/blastn."$query"."$species2".e-value0.05.outfmt7.out -num_threads 32 -outfmt 7 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -word_size 7 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/word_size7_outfmt1/blastn."$query"."$species2".word_size7.e-value0.05.outfmt1.out -num_threads 32 -outfmt 1 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -word_size 7 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/word_size7_outfmt7/blastn."$query"."$species2".word_size7.e-value0.05.outfmt7.out -num_threads 32 -outfmt 7 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

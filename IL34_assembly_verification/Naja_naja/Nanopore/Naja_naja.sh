********* Nanopore **************************************************************************************************************

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
  /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/minimap2-2.28_x64-linux/minimap2 -ax map-ont -t 36 ../"$species2".fa $fastq | samtools sort -@36 -O BAM - > ${fastq%.fastq.gz}.bam                  
done

samtools merge -@36 Nanopore_merged.bam *.bam 
samtools sort -@36 Nanopore_merged.bam -o packbio_merged.sorted.bam
samtools index -@36 Nanopore_merged.sorted.bam


#######________Klumpy____________________###############

conda activate klumpy_env
klumpy scan_alignments --alignment_map Nanopore_merged.sorted.bam  --threads 32  --annotation /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Naja_naja/GCA_009733165.1_Nana_v5_genomic.gtf 

klumpy alignment_plot --alignment_map Nanopore_merged.sorted.bam  --reference CM019161.1 --candidates  Nanopore_merged.sorted_Candidate_Regions.tsv--window_size 10000 --window_step 5000 --color red --vertical_line_gaps --vertical_line_klumps --format svg --leftbound 768622 --rightbound 781041 --annotation /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Naja_naja/GCA_009733165.1_Nana_v5_genomic.gtf  --gap_file Naja_naja_gaps.tsv  --height 1500 --width 4000 --min_len 10000   --vertical_line_exons --numberBLAST


 ######################################################################################################################
species="Naja_naja"  
zcat *fastq.gz |sed  -n '1~4s/^@/>/p;2~4p' > "$species".fa
makeblastdb -in "$species".fa -out "$species".fa -dbtype nucl

species3="Naja_naja" 
genome=/media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Naja_naja/Nanopore/Naja_naja.fa 
for query in IL34_query_multispecies.fa Sphenodon_punctatus_IL34_Exons.fa
do
blastn -task blastn -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/nano_blast/outfmt1/blastn.nanopore."$query"."$species3".e-value0.05.outfmt1.out -num_threads 32 -outfmt 1 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/nano_blast/outfmt7/blastn.nanopore."$query"."$species3".e-value0.05.outfmt7.out -num_threads 32 -outfmt 7 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -word_size 7 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/nano_blast/word_size7_outfmt1/blastn.nanopore."$query"."$species3".word_size7.e-value0.05.outfmt1.out -num_threads 32 -outfmt 1 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -word_size 7 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/nano_blast/word_size7_outfmt7/blastn.nanopore."$query"."$species3".word_size7.e-value0.05.outfmt7.out -num_threads 32 -outfmt 7 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

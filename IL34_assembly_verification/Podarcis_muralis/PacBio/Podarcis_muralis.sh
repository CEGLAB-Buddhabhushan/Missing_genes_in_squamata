************** pacbio hifi ***************************************************************************************************
chromosome downlaod 
chr=$(cat co_ordinates.txt | cut -f1 -d':') && esearch -db nucleotide -query "$chr" | efetch -format fasta > "Podarcis_muralis.fa"


####________downlaoding data__________####################

conda activate pfd

for i in `cat SRA.lst`
do
parallel-fastq-dump --sra-id $i --threads 24 --outdir . --gzip
done

####_________mapping____minimap________####################
species3="Podarcis_muralis"  

for fastq in *.fastq.gz
do
  /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/minimap2-2.28_x64-linux/minimap2 -ax map-hifi -t 36 ../"$species3".fa $fastq | samtools sort -@36 -O BAM - > ${fastq%.fastq.gz}.bam
  samtools index ${fastq%.fastq.gz}.bam                   
done


####__________klumpy___________________####################

conda activate klumpy_env
klumpy scan_alignments --alignment_map SRR8468524_subreads.bam   --threads 32  --annotation /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Podarcis_muralis/GCF_004329235.1_PodMur_1.0_genomic.gtf

klumpy find_gaps --fasta  /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Podarcis_muralis/Podarcis_muralis.fa

klumpy alignment_plot --alignment_map SRR8468524_subreads.bam --reference NC_041319.1 --candidates  SRR8468524_subreads_Candidate_Regions.tsv  --window_size 10000 --window_step 5000 --color red --vertical_line_gaps --vertical_line_klumps --format svg --leftbound 39114790 --rightbound 39128268 --annotation /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Podarcis_muralis/GCF_004329235.1_PodMur_1.0_genomic.gtf --gap_file Podarcis_muralis_gaps.tsv  --height 1500 --width 3500 --min_len 10000 --vertical_line_exons   --number


BLAST ######################################################################################################################
species="Podarcis_muralis"  
zcat *fastq.gz |sed  -n '1~4s/^@/>/p;2~4p' > "$species".fa
makeblastdb -in "$species".fa -out "$species".fa -dbtype nucl


######________blast and store outfmt 1 and 7_____####### 
species3="Podarcis_muralis" 
genome=/media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/Podarcis_muralis/PacBio/Podarcis_muralis.fa 
for query in IL34_query_multispecies.fa Sphenodon_punctatus_IL34_Exons.fa
do
blastn -task blastn -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/outfmt1/blastn."$query"."$species3".e-value0.05.outfmt1.out -num_threads 32 -outfmt 1 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/outfmt7/blastn."$query"."$species3".e-value0.05.outfmt7.out -num_threads 32 -outfmt 7 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -word_size 7 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/word_size7_outfmt1/blastn."$query"."$species3".word_size7.e-value0.05.outfmt1.out -num_threads 32 -outfmt 1 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

blastn -task blastn -evalue 0.05 -word_size 7 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/BLASTn_output/word_size7_outfmt7/blastn."$query"."$species3".word_size7.e-value0.05.outfmt7.out -num_threads 32 -outfmt 7 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/$query

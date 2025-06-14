#map-pb map-ont map-hifi
species=$1
data_type=$2
minimap2_flag=$3
cd $species
chr=`cat co_ordinates.txt|cut -f1 -d':'`
esearch -db nucleotide -query "$chr"|efetch -format fasta > "$species".fa

cd $data_type

conda activate pfd

for i in `cat SRA.lst`
do
#parallel-fastq-dump --sra-id $i --threads 16 --outdir . --split-files --gzip
parallel-fastq-dump --sra-id $i --threads 16 --outdir . --gzip
done


for fastq in *.fastq.gz
do
/media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_long_read_verification/minimap2-2.28_x64-linux/minimap2 -ax $minimap2_flag -t 36 ../"$species".fa $fastq | samtools sort -O BAM - > "$species".bam
samtools index "$species".bam
done

## makeblastdb
for i in *_1.fastq.gz
do
f=${i%_1.fastq.gz}_2.fastq.gz
zcat $i $f |sed  -n '1~4s/^@/>/p;2~4p' > "$species".fa

makeblastdb -in "$species".fa -out "$species".fa -dbtype nucl
done

zcat *fastq.gz |sed  -n '1~4s/^@/>/p;2~4p' > "$species".fa
makeblastdb -in "$species".fa -out "$species".fa -dbtype nucl


samtools merge packbio_merged.bam *.bam
samtools sort packbio_merged.bam -o packbio_merged.sorted.bam
samtools index packbio_merged.sorted.bam

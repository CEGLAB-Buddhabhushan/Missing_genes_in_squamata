# download  orthologs of vertebrate from NCBI
#get accessions and download longets isoform using cdskit
for gene in LAPTM5 STAP1 TNIP2
do
cd $gene
rm *.fa *.aggregate.fasta tmp.fa *.accession2fasta.fasta *.acc.txt
for i in *_refseq_transcript.fasta
do
clade=`echo $i|sed 's/_refseq_transcript.fasta//g'`
grep ">" $i|cut -f1 -d' '|sed 's/>//g' > "$clade".acc.txt
grep -v "XM" "$clade".acc.txt ## mannualy download and add this accession id
cdskit accession2fasta --accession_file "$clade".acc.txt -o "$clade".accession2fasta.fasta
cdskit aggregate --seqfile "$clade".accession2fasta.fasta --outfile "$clade".aggregate.fasta
cut -f1,2 -d'_' "$clade".aggregate.fasta|seqtk seq > "$clade".fa
for i in `grep "NM_" "$clade".acc.txt`
do
accession=$i
species_name=`grep "$accession" "$clade"_refseq_transcript.fasta|cut -f2,3 -d' '|sed 's/ /_/g'`
efetch -db nuccore -id "$accession" -format fasta_cds_na | sed "s/^>\(.*\) \[gene=.*$/>\1 $species_name/"|sed 's/^>\S*_cds_\S*_\S*_\S* />/' >> "$clade".fa
done
cdskit aggregate --seqfile "$clade".fa --outfile "$clade".fa.aggregate.fasta
cut -f1,2 -d'_' "$clade".fa.aggregate.fasta|seqtk seq > tmp.fa
mv tmp.fa "$clade".fa
done
cd ..
done

for i in *.fa
do
perl ~/Documents/Gc_content/ORFvalidator.pl $i
done


mkdir output
for gene in LAPTM5 STAP1 TNIP2
do
cp $gene/*.fa output
done
cd output

echo -e "Gene_name\tGroup\tSpecies_name\tGC_content\tGC_Stretch" > GC_content-GC_Stretch.tsv
for i in *.fa
do
Gene=`echo $i|cut -f2 -d'.'`
Group=`echo $i|cut -f1 -d'.'`
seqkit fx2tab --name --gc $i > "$i".gc.tsv
perl ../GC_Stretch_finder.pl $i > "$i".stretch.tsv
for Species_name in `cut -f1 "$i".gc.tsv`
do
GC_content=`grep "$Species_name" "$i".gc.tsv|cut -f2`
GC_Stretch=`grep "$Species_name" "$i".stretch.tsv|cut -f2 -d','|sed 's/Average Length: //g'`
echo -e "$Gene\t$Group\t$Species_name\t$GC_content\t$GC_Stretch" >> GC_content-GC_Stretch.tsv
done
rm "$i".gc.tsv "$i".stretch.tsv
done
## Generate plot using R script... this can be found in respective folder

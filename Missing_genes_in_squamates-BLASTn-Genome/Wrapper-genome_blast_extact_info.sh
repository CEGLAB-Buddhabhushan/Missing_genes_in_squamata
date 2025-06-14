curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat

cp ../HybPiper/*.fa .
#keep human gene sequence only
for i in `ls -1 *.fa`; do seqtk seq $i |head -2 >tmp; mv tmp $i; rm tmp; done
mkdir -p BLASTn_output
for i in `cat Squamata_Reference_scaffold_chromososme.sorted.lst`
do
accession=`echo $i|cut -f1 -d'-'`
Species_name=`echo $i|cut -f2 -d'-'`

./datasets download genome accession $accession --include genome --filename "$i".zip --api-key 61131f45442b718e60396eb45ce73b65bc08 --no-progressbar
rm -rf $i
unzip -q -o "$i".zip -d $i
rm "$i".zip
cd $i/ncbi_dataset/data/$accession/
mkdir -p /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/BLASTn_output/$i
genome=`ls -1 *.fna|head -1`
makeblastdb -in $genome -out $genome -dbtype nucl
for query in ADAP2.fa ADRA1B.fa ARL11.fa CAV3.fa CYYR1.fa DIP2A.fa DLEU7.fa DNM3.fa ELOVL3.fa FREM2.fa GPR78.fa GPR82.fa GRK4.fa HHLA2.fa HPGDS.fa HSD17B1.fa IFTAP.fa IL13RA2.fa IL26.fa IL34.fa IL5RA.fa INTS6.fa KDM3A.fa KREMEN1.fa LAPTM5.fa LCP1.fa LRCH1.fa MAN2B2.fa MATK.fa MOB1A.fa NOTCH2.fa OLAH.fa PLAAT1.fa RBBP7.fa RIPPLY3.fa RUBCNL.fa SH2D1A.fa SH2D2A.fa SIAH3.fa SLC24A1.fa SLC9A5.fa SSTR4.fa STAP1.fa STOML3.fa SUV39H2.fa TBC1D14.fa TMEM273.fa TNIP2.fa UNKL.fa UTS2B.fa WNT2.fa YBX3.fa ZNF438.fa
do
gene=`echo $query|sed 's/\.fa//g'`
blastn -task blastn -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/BLASTn_output/$i/blastn."$gene"."$i".e-value0.05.outfmt1.out -num_threads 32 -outfmt 1 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/$query
blastn -task blastn -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/BLASTn_output/$i/blastn."$gene"."$i".e-value0.05.outfmt7.bls -num_threads 32 -outfmt 7 -query /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/$query
blastn -task blastn  -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/BLASTn_output/$i/blastn."$gene"."$i".e-value0.05.outfmt6.tsv -num_threads 32 -outfmt '6 qseqid sseqid evalue qlen qcovs' -query /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/$query
blastn -task blastn  -evalue 0.05 -db $genome -out /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/BLASTn_output/$i/blastn."$gene"."$i".e-value0.05.outfmt17.sam -num_threads 32 -outfmt '17 SQ' -query /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/$query
done

cd /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/
rm -r $i
echo $Species_name
done

 
cd BLASTn_output
# remove empty directory
find . -type d -empty -delete
## genome blast done on 265
# copy file for each gene

mkdir -p /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
mkdir -p /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise/$gene
for i in `ls -1`
do
cp $i/*"$gene"* /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise/$gene/
echo $gene $i
done
done
cd ../GENE_wise/
# create IGV report
mkdir Genome_BLASTn_IGV_report
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
cd $gene
cp /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/"$gene".fa .
header=`grep ">" "$gene".fa |sed 's/>//g'`

    empty_sam_files=$(find . -maxdepth 1 -name "*.sam" -size 0 | wc -l)
    total_sam_files=$(find . -maxdepth 1 -name "*.sam" | wc -l)

    if [ "$total_sam_files" -eq 0 ]; then
        echo "$gene: No SAM files found" >> /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise/Genome_BLASTn_IGV_report.output.log
        cd /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise/
        continue
    fi

    if [ "$empty_sam_files" -eq "$total_sam_files" ]; then
        echo -e "$gene\t$total_sam_files\t$empty_sam_files\tNo BLASTn hits found" >> /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise/Genome_BLASTn_IGV_report.output.log
        cd /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise/
        continue
    fi

for sam in `ls -1 *.sam`
do
species_name=`echo $sam|cut -f2 -d'-'|sed 's/\.e-value0.05.outfmt17.sam//g'`
sed -i "s/Query_1/$header/g" $sam
samtools view -bhS $sam > "$species_name"-"$gene".bam
samtools sort "$species_name"-"$gene".bam -o "$species_name"-"$gene".sorted.bam
samtools index "$species_name"-"$gene".sorted.bam
done
faidx --transform bed "$gene".fa > "$gene".bed
bam_list=`ls -1 *.sorted.bam|tr '\n' ' '`
create_report "$gene".bed --standalone --fasta "$gene".fa --tracks $bam_list --output "$gene".Genome_BLASTn.IGV_report.html --translate-sequence-track
cp "$gene".Genome_BLASTn.IGV_report.html /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise/Genome_BLASTn_IGV_report/
cd /media/lokdeep/sdf/BUDDHA/BMCBIO/other_gene_loss/Squamata_genone_blast/GENE_wise/
done

echo -e "Gene_name Query_Seq-id Subject_Seq-id Expect_value Query_sequence_length Query_Coverage_Per_Subject"|sed 's/ /\t/g' > Maximum_query_coverage.tsv
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
max_cov=`cat $gene/*.tsv|sort -nrk5,5|head -1`
echo -e "$gene\t$max_cov" >> Maximum_query_coverage.tsv
echo $gene
done
#gene + species wise
echo -e "Gene_name\tSpecies_name\tMax_cov" > Gene_species_Max_coverage.tsv
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
cd $gene
for tsv in `ls -1 *.outfmt6.tsv`
do
species_name=`echo $tsv|cut -f2 -d'-'|sed 's/\.e-value0.05.outfmt6.tsv//g;s/\.e//g'`
max_cov=`cat $tsv|sort -nrk5,5|head -1|cut -f5`
echo -e "$gene\t$species_name\t$max_cov" >> ../Gene_species_Max_coverage.tsv
done
cd ..
echo $gene
done





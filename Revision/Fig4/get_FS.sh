for gene in ADAP2 ADRA1B ARL11 C11orf74 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HRASLS HSD17B1 IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 C10orf128 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
grep -w "$gene" -A10 -B10 Gallus_gallus.gene.bed|cut -f4|sort -u >> focal_and_syntenic_gene.lst
done


for gene in `cat focal_and_syntenic_gene.lst`
do
id=`grep -w "$gene" longest_isoform.lst|cut -f1 -d':'`
awk '$8 == "exon"' GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gtf2bed.bed| grep -w "$id" |cut -f1-4 > FS_"$gene".exon.bed
awk '$8 == "CDS"' GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gtf2bed.bed| grep -w "$id" |cut -f1-4 >  FS_"$gene".CDS.bed

bedtools intersect -a FS_"$gene".CDS.bed -b FS_"$gene".exon.bed >> FS_Exonic.Gallus_gallus.GRCg7b.longest_isoform.All_genes.bed
rm FS_"$gene".CDS.bed FS_"$gene".exon.bed
echo $id $gene
done

for i in `cut -f2 -d'-' chr_name.lst`
do
rename_chr=`grep "$i"  chr_name.lst|cut -f1 -d'-'`
sed -i "s/$i/$rename_chr/g" FS_Exonic.Gallus_gallus.GRCg7b.longest_isoform.All_genes.bed
done
grep -v "scaffold" FS_Exonic.Gallus_gallus.GRCg7b.longest_isoform.All_genes.bed | grep -v "NC_"|grep -v "NW_"> FS_Exonic.Gallus_gallus.GRCg7b.longest_isoform.All_genes.reformated.bed


for i in *.chicken.bed
do
sp=$(echo $i | sed 's/\.chicken.bed//g')
bedtools intersect -a FS_Exonic.Gallus_gallus.GRCg7b.longest_isoform.All_genes.reformated.bed -b $i -wao | awk -v sp=$sp '{print $0"\t"sp}' >> FS_Exonic.Gallus_gallus.GRCg7b.longest_isoform.all_genes.overlap.bed
echo "Processed: $sp"
done



bedtools groupby -i FS_Exonic.Gallus_gallus.GRCg7b.longest_isoform.all_genes.overlap.bed -g 1,2,3,4,9 -c 8 -o sum > FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed
sort -u FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed > tmp.bed
mv tmp.bed FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed 
sed -i 's/Anas_platyrhynchos_platyrhynchos/Anas_platyrhynchos/g' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed
sed -i 's/Aquila_chrysaetos_chrysaetos/Aquila_chrysaetos/g' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed
sed -i 's/Chrysemys_picta_bellii/Chrysemys_picta/g' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed
sed -i 's/Struthio_camelus_australis/Struthio_camelus/g' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed
sed -i 's/Terrapene_carolina_triunguis/Terrapene_carolina/g' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed


### For focal gene
echo -e "Gene_name\tClass\tOrder\tSpecies_name\tAligned_region\tGene_length\tNormalized_value" > FS_Exonic.Final.Gallus_gallus.53_genes.focal.tsv

for gene in ADAP2 ADRA1B ARL11 C11orf74 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HRASLS HSD17B1 IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 C10orf128 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
for species in `cut -f5 FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed|sort -u`
do
class=`awk -v sp="$species" '$3 == sp {print $1}' taxonomy_info.tsv`
order=`awk -v sp="$species" '$3 == sp {print $2}' taxonomy_info.tsv`
aligned_region=$(awk -v sp="$species" -v g="$gene" '$5 == sp && $4 == g {sum+=$6} END{print sum}' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed)
Gene_length=$(awk -v sp="$species" -v g="$gene" '
    $5 == sp && $4 == g { sum += ($3 - $2) }
    END { print sum }
' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed)
Normalized_value=$(awk -v ur="$aligned_region" -v gl="$Gene_length" 'BEGIN{if(gl>0) print ur/gl; else print "NA"}')
echo -e "$gene\t$class\t$order\t$species\t$aligned_region\t$Gene_length\t$Normalized_value" >> FS_Exonic.Final.Gallus_gallus.53_genes.focal.tsv
echo -e "$gene\t$class\t$order\t$species\t$aligned_region\t$Gene_length\t$Normalized_value" 
done
done

### focal syntenic genes
for focal_gene in ADAP2 ADRA1B ARL11 C11orf74 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HRASLS HSD17B1 IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 C10orf128 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
echo -e "Gene_name\tClass\tOrder\tSpecies_name\tAligned_region\tGene_length\tNormalized_value" > "$focal_gene".FS_Exonic.Final.Gallus_gallus.21_genes.tsv
for gene in `awk -v FG="$focal_gene" '$11 == FG {print $0}' syntenic_genes.lst`
do
for species in `cut -f5 FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed |sort -u`
do
class=`awk -v sp="$species" '$3 == sp {print $1}' taxonomy_info.tsv`
order=`awk -v sp="$species" '$3 == sp {print $2}' taxonomy_info.tsv`
aligned_region=$(awk -v sp="$species" -v g="$gene" '$5 == sp && $4 == g {sum+=$6} END{print sum}' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed )
Gene_length=$(awk -v sp="$species" -v g="$gene" '
    $5 == sp && $4 == g { sum += ($3 - $2) }
    END { print sum }
' FS_Exonic.Gallus_gallus.All_genes.merged_overlap.bed)
Normalized_value=$(awk -v ur="$aligned_region" -v gl="$Gene_length" 'BEGIN{if(gl>0) print ur/gl; else print "NA"}')
echo -e "$gene\t$class\t$order\t$species\t$aligned_region\t$Gene_length\t$Normalized_value" >> "$focal_gene".FS_Exonic.Final.Gallus_gallus.21_genes.tsv
done
done
echo $focal_gene 
done



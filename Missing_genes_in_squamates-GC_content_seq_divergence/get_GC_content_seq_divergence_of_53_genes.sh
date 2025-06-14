### get gc content, gc strenches and sequence divergence
rm -r GC_content_seq_divergence
mkdir GC_content_seq_divergence
echo -e "Gene_name\tSample_size\tGC_content\tGC_Stretch\tSeq_iden" > GC_content_seq_divergence/Gene_sequence_properties.tsv
for i in `cat Complete_loss.lst`
do
gene_id=`echo $i|cut -f1 -d':'`
gene=`echo $i|cut -f2 -d':'`
cut -f2,4,5 $gene_id/ncbi_dataset/data/transcript_protein.tsv|sed 's/ /_/g;s/\t/-/g'|tail -n+2 > GC_content_seq_divergence/"$gene".transcript_info.lst
cd GC_content_seq_divergence/
for i in `cat "$gene".transcript_info.lst|sort -u`
do
species_name=`echo $i|cut -f1 -d'-'`
accession=`echo $i|cut -f2 -d'-'`
length=`echo $i|cut -f3 -d'-'`
efetch -db nuccore -id "$accession" -format fasta_cds_na | sed "s/^>\(.*\) \[gene=.*$/>\1 $species_name:$accession:$length/"|sed 's/^>\S*_cds_\S*_\S*_\S* />/' >> "$gene".fa
done
cdskit aggregate --seqfile "$gene".fa --outfile "$gene".longest.fa --expression ":.*"
perl ../mulfifasta.pl "$gene".longest.fa |cut -f1 -d':'|seqtk seq > "$gene".valid.longest.fa
Sample_size=`grep ">" "$gene".valid.longest.fa |wc -l`
GC_content_mean_sd=`seqkit fx2tab --name --gc "$gene".valid.longest.fa| awk '{sum+=$2; sumsq+=$2*$2; n+=1} END {mean=sum/n; stddev=sqrt((sumsq - (sum^2)/n)/n); print  mean, "+/-", stddev}'`
GC_Stretch_mean_sd=`perl ../GC_Stretch_finder.pl "$gene".valid.longest.fa | cut -f5 -d' ' | awk '{sum+=$1; sumsq+=$1*$1; n+=1} 
     END {
         mean=sum/n; 
         stddev=sqrt((sumsq - (sum^2)/n)/n); 
         print mean, "+/-", stddev
     }'`
mkdir "$gene"_query/
faSplit byname "$gene".valid.longest.fa "$gene"_query/
echo "species ident"|sed 's/ /\t/g' > "$gene".identity_wrt_human.out
cd  "$gene"_query/
for i in *.fa
do
n=`echo $i|sed 's/\.fa//g'`
needle -asequence Homo_sapiens.fa -bsequence $i -outfile Homo_sapiens-"$n".needle.out -gapopen 10.0 -gapextend 0.5
sim=`grep "Similarity" Homo_sapiens-"$n".needle.out |awk '{print $4}'|sed 's/(//g'|sed 's/)//g'|sed 's/%//g'`
ident=`grep "Identity" Homo_sapiens-"$n".needle.out |awk '{print $4}'|sed 's/(//g'|sed 's/)//g'|sed 's/%//g'`
gaps=`grep "Gaps" Homo_sapiens-"$n".needle.out |awk '{print $5}'|sed 's/(//g'|sed 's/)//g'|sed 's/%//g'`
echo $n $ident|sed 's/ /\t/g' >> ../"$gene".identity_wrt_human.out
done
cd ..
Seq_iden_mean_sd=`cat "$gene".identity_wrt_human.out| awk '{sum+=$2; sumsq+=$2*$2; n+=1} END {mean=sum/n; stddev=sqrt((sumsq - (sum^2)/n)/n); print  mean, "+/-", stddev}'`
cd ..
echo -e "$gene\t$Sample_size\t$GC_content_mean_sd\t$GC_Stretch_mean_sd\t$Seq_iden_mean_sd" >> GC_content_seq_divergence/Gene_sequence_properties.tsv
echo $gene
done

echo -e "Gene_name\tSpecies_name\tSeq.Identity" > Identity_wrt_human.tsv
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
awk -v gene="$gene" 'BEGIN {OFS="\t"} NR==1 {print "Gene_name", $0} NR>1 && NF==2 {print gene, $0}' "${gene}.identity_wrt_human.out" |grep -v "species">> Identity_wrt_human.tsv
done
#########

echo -e "Gene_name\tSpecies_name\tGC_content" > GC_content.tsv
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
seqkit fx2tab --name --gc "$gene".valid.longest.fa \
    | sed 's/:/\t/g' \
    | awk -v gene="$gene" 'NR>1 {print gene, $0}' OFS="\t" \
    >> GC_content.tsv
done

####
echo -e "Gene_name\tSpecies_name\tGC_Stretch" > GC_Stretch.tsv
for gene in ADAP2 ADRA1B ARL11 CAV3 CYYR1 DIP2A DLEU7 DNM3 ELOVL3 FREM2 GPR78 GPR82 GRK4 HHLA2 HPGDS HSD17B1 IFTAP IL13RA2 IL26 IL34 IL5RA INTS6 KDM3A KREMEN1 LAPTM5 LCP1 LRCH1 MAN2B2 MATK MOB1A NOTCH2 OLAH PLAAT1 RBBP7 RIPPLY3 RUBCNL SH2D1A SH2D2A SIAH3 SLC24A1 SLC9A5 SSTR4 STAP1 STOML3 SUV39H2 TBC1D14 TMEM273 TNIP2 UNKL UTS2B WNT2 YBX3 ZNF438
do
perl GC_Stretch_finder.pl "$gene".valid.longest.fa|cut -f2,5 -d' '|sed 's/,//g;s/ /\t/g' \
    | awk -v gene="$gene" 'NR>1 {print gene, $0}' OFS="\t" \
    >> GC_Stretch.tsv
done

